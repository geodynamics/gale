import os, sys
import SCons.Script
import SConfig
from Package import Package
from Installation import Installation
from Configuration import Configuration

class Project(SConfig.Node):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Node.__init__(self, scons_env, scons_opts, required)
        self.checks += [self.check_libs, self.select_config, self.print_results]

    def setup_options(self):
        self.opts.AddOptions(
            SCons.Script.BoolOption('with_debug',
                                    'Generate debugging symbols', 0),
            SCons.Script.BoolOption('static_libraries',
                                    'Build static libraries', 1),
            SCons.Script.BoolOption('shared_libraries',
                                    'Build shared libraries', 1),
            SCons.Script.BoolOption('with_optimise',
                                    'optimising with -03 flag', 1),
            SCons.Script.BoolOption('with_eptiming',
                                    'enable parallel EP execution timing', 0),
            SCons.Script.BoolOption('with_warnings',
                                    'enable compile warnings', 0),	    
            SCons.Script.BoolOption('with_tau',
                                    'use tau instrumentation', 0),
            ('tau_cc', 'tau compiler to use if tau is enabled', 'tau_cc.sh'),
            ('build_dir', 'temporary build directory', 'build')
            )

    def check_libs(self):
        if not self.env['static_libraries'] and not self.env['shared_libraries']:
            self.ctx.Display("      Both static and shared libraries disabled!\n")
            return False
        return True

    def select_config(self):
        """Decide which set of dependencies to select as the desired
        configuration."""

        # Don't know why, but I need to do the following lengthy thing because
        # this doesn't work: deps = [d for d in self.reduce_dependencies()]
        deps = []
        clashes = []
        for d in self.reduce_dependencies():
            if isinstance(d, tuple):
                clashes += [d]
            else:
                cur = [s for s in d]
                if cur:
                    deps += [cur]

        # If there were no valid configurations found, then we're pretty much
        # screwed. Display the clashes so the user knows why.
        if not len(deps):
            self.ctx.Display('  There are no valid combinations of the required packages.\n')
            self.ctx.Display('  This means that one or more packages have conflicting\n')
            self.ctx.Display('  dependencies.\n\n')
            if not len(clashes):
                raise 'Error: Should be clashes!'

            # Need to process the list of clashes and construct some useful information.
            clash_dict = {}
            for clash in clashes:
                if clash[2].inst.pkg != clash[3].inst.pkg:
                    raise 'Error: Inconsistent clash report!'
                cause = clash[2].inst.pkg
                top_pkgs = (clash[0].inst.pkg, clash[1].inst.pkg)
                if cause not in clash_dict:
                    clash_dict[cause] = {}
                top_pkg_dict = clash_dict[cause]
                for top_pkg, useable in zip(top_pkgs, clash[2:]):
                    if top_pkg == useable.inst.pkg:
                        continue
                    if top_pkg not in top_pkg_dict:
                        top_pkg_dict[top_pkg] = []
                    if useable not in top_pkg_dict[top_pkg]:
                        top_pkg_dict[top_pkg] += [useable]

            for cause_pkg, useable_dict in clash_dict.iteritems():
                self.ctx.Display('  There were conflicts with package \'%s\':\n' % cause_pkg.name)
                for useable_pkg, useable_list in useable_dict.iteritems():
                    self.ctx.Display('    %s can use:\n' % useable_pkg.name)
                    for useable in useable_list:
                        self.ctx.Display('      %s\n' % useable.inst.base_dir)
            return False

        # Decide which selection of dependencies to use. We use the 'value'
        # method for each configuration.
        max_val = 0
        selected = deps[0]
        for d in deps:
            cur_val = 0
            for cfg in d:
                if isinstance(cfg, Configuration):
                    cur_val += cfg.value()
            if cur_val > max_val:
                max_val = cur_val
                selected = d

        # Set the 'selected' member of each package with configurations.
        done = []
        rem = list(selected)
        while len(rem):
            cur = rem.pop()
            if isinstance(cur, Configuration):
                if cur.inst.pkg not in done:
                    cur.inst.pkg.selected = cur
                    rem += cur.deps
                    done += [cur.inst.pkg]
        return True

    def reduce_dependencies(self, deps=[], pkgs={}, cur_index=0):
        """Each combination of dependent installations represents a
        unique installation of this package. This method generates sets
        of unique dependency combinations."""

        # The dictionary 'pkgs' is to map packages to installations. This is
        # needed to prevent multiple installations being used for the same
        # packages.

        # Once we're at the end of the list of dependencies we can
        # yield a result.
        if cur_index == len(self.deps):
            yield deps

        else:
            cur_dep, required = self.deps[cur_index]

            # If the dependency isn't actually required by this package, include
            # a combination that doesn't use it.
            if not required:
                for d in self.reduce_dependencies(deps, pkgs, cur_index + 1):
                    yield d

            # 'Package's are the only kind of object to have installations.
            if isinstance(cur_dep, Package):

                # If we've already included an installation of this package
                # in this set of dependencies we must use the same one in
                # every other dependency.
                if cur_dep in pkgs:
                    deps += [pkgs[cur_dep]]
                    for d in self.reduce_dependencies(deps, pkgs, cur_index + 1):
                        yield d
                    del deps[-1]

                # If we havn't already foudn this dependency anywhere, try
                # out all of it's installations, if it has any.
                elif len(cur_dep.configurations):
                    for dep_cfg in cur_dep.configurations:

                        # We have to keep the current state of 'pkgs' clean,
                        # so copy it to a new one.
                        new_pkgs = dict(pkgs)

                        # Traverse the dependency tree for this installation,
                        # adding each dependent installation to the 'new_pkgs'
                        # mapping.
                        rem = [dep_cfg]
                        okay = True
                        while len(rem):
                            cur = rem.pop()

                            # If the installation doesn't exist in the mapping, just
                            # add it in and continue.
                            if cur.inst.pkg not in new_pkgs:
                                new_pkgs[cur.inst.pkg] = (cur, dep_cfg)

                            # If we've already got a version of this dependency
                            # in the mapping, it means the version we're currently
                            # looking at must be the same or we can't use this
                            # combination.
                            else:
                                entry = new_pkgs[cur.inst.pkg]
                                if entry[0] != cur:

                                    # Report conflicts using a tuple of both the
                                    # conflicting installations.
                                    yield (dep_cfg, entry[1], cur, entry[0])
                                    okay = False

                                    # Get out of this loop, we can't use this
                                    # dependency combination now. Move on to the
                                    # next installation.
                                    break

                            # Add the dependencies of the current installation to be
                            # checked.
                            rem += cur.deps

                        # If we were able to match up all previous installations with
                        # this one, continue recursing for every other dependency.
                        if okay:
                            deps += [dep_cfg]
                            for d in self.reduce_dependencies(deps, new_pkgs, cur_index + 1):
                                yield d
                            del deps[-1]

                # There are no installations for this dependency. If it's
                # a required dependency then something has gone very wrong.
                elif cur_dep.required:

                    # Throw an error here.
                    # TODO
                    raise 'Error: Should be a valid configuration.'

                else:
                    # There are no installations for this package but it's also not required,
                    # so we can safely ignore it and continue on our merry way.
                    for d in self.reduce_dependencies(deps, pkgs, cur_index + 1):
                        yield d

            else:
                deps += [cur_dep]
                for d in self.reduce_dependencies(deps, pkgs, cur_index + 1):
                    yield d
                del deps[-1]

    def print_results(self):
        self.ctx.Display("  Static libraries: %s\n" % str(bool(self.env['static_libraries'])))
        self.ctx.Display("  Shared libraries: %s\n" % str(bool(self.env['shared_libraries'])))
        self.ctx.Display("  Using build directory: %s\n" % self.env['build_dir'])
        self.ctx.Display("  Debugging symbols: %s\n" % str(bool(self.env['with_debug'])))
        self.ctx.Display("  Optimisation symbols: %s\n" % str(bool(self.env['with_optimise'])))
        if self.env['with_tau']:
            self.ctx.Display('  Instrumenting with tau compiler: %s\n' % self.env['tau_cc'])
        return True

    def enable(self, scons_env, old_state=None):
        SConfig.Node.enable(self, scons_env, old_state)

        # Setup tau.
        if self.env['with_tau']:
            self.backup_variable(scons_env, 'tau_old_cc', old_state)
            self.backup_variable(scons_env, 'CC', old_state)
            scons_env['tau_old_cc'] = self.env['CC']
            scons_env['CC'] = self.env['tau_cc']
            scons_env.AppendUnique(CONFIGVARS=['tau_old_cc'])

        # Setup debugging flags.
        if self.env['with_debug']:
            d = scons_env.ParseFlags('-g')
            self.backup_variable(scons_env, d.keys(), old_state)
            scons_env.MergeFlags(d)
            d = scons_env.ParseFlags('-DDEBUG')
            self.backup_variable(scons_env, d.keys(), old_state)
            scons_env.MergeFlags(d)
        else:
            d = scons_env.ParseFlags('-DNDEBUG')
            self.backup_variable(scons_env, d.keys(), old_state)
            scons_env.MergeFlags(d)

        # Setup the optimisation flags
        if self.env['with_optimise'] and not self.env['with_debug']:
            d = scons_env.ParseFlags('-O3')
            self.backup_variable(scons_env, d.keys(), old_state)
            scons_env.MergeFlags(d)

	# Setup the warnings flags
        if self.env['with_warnings']:
            d = scons_env.ParseFlags('-Wall')
            self.backup_variable(scons_env, d.keys(), old_state)
            scons_env.MergeFlags(d)

        # Setup the ep timing flags
        if self.env['with_eptiming']:
            d = scons_env.ParseFlags('-DENABLE_STGERMAIN_LOG')
            self.backup_variable(scons_env, d.keys(), old_state)
            scons_env.MergeFlags(d)

        # Setup the include paths.
        inc_dir = self.env.get_build_path('include')
        self.backup_variable(scons_env, 'CPPPATH', old_state)
        scons_env.PrependUnique(CPPPATH=[inc_dir])

        # Setup LIB_DIR.
        lib_dir = self.env.get_build_path('lib')
        self.backup_variable(scons_env, 'LIBPATH', old_state)
        scons_env.PrependUnique(LIBPATH=[lib_dir])

        # Setup the RPATH.
        self.backup_variable(scons_env, 'RPATH', old_state)
        scons_env.PrependUnique(RPATH=[scons_env.Dir(lib_dir).abspath])

	# TODO print compiler flags without quotes and commas so they can be copied clenly
	self.ctx.Display( "Compiler flags: " )
	print scons_env['CCFLAGS'][:]
	self.ctx.Display( "Precompiler defines" )
	print scons_env['CPPDEFINES'][:]
