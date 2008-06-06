import os, sys
import SCons.Script
import SConfig
from Package import Package
from Package import Installation

class Project(SConfig.Node):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Node.__init__(self, scons_env, scons_opts, required)
        self.checks += [self.check_libs, self.select_config, self.print_results]

    def setup_options(self):
        self.opts.AddOptions(
            SCons.Script.BoolOption('with_debug',
                                    'Generate debugging symbols', 1),
            SCons.Script.BoolOption('static_libraries',
                                    'Build static libraries', 1),
            SCons.Script.BoolOption('shared_libraries',
                                    'Build shared libraries', 1),
            ('build_dir', 'Temporary build directory', 'build')
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
        for d in self.reduce_dependencies():
            cur = [s for s in d]
            if cur:
                deps += [cur]

        # Need to determine the best set of dependencies to use. For now
        # just take the first one we can find.
        if len(deps):
            found = False
            for selected in deps:
                if isinstance(selected, tuple):
                    # A tuple represents a conflicting set of packages.
                    continue
                found = True
                done = []
                rem = list(selected)
                while len(rem):
                    cur = rem.pop()
                    if isinstance(cur, Installation):
                        if cur.pkg not in done:
                            cur.pkg.selected = cur
                            rem += cur.deps
                            done += [cur.pkg]
            if not found:
                return False

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
            cur_dep = self.deps[cur_index][0]

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
                elif len(cur_dep.installations):
                    for dep_inst in cur_dep.installations:

                        # We have to keep the current state of 'pkgs' clean,
                        # so copy it to a new one.
                        new_pkgs = dict(pkgs)

                        # Traverse the dependency tree for this installation,
                        # adding each dependent installation to the 'new_pkgs'
                        # mapping.
                        rem = [dep_inst]
                        okay = True
                        while len(rem):
                            cur = rem.pop()

                            # We can't do this unless the dependency is an installation.
                            if isinstance(cur, Installation):

                                # If the installation doesn't exist in the mapping, just
                                # add it in and continue.
                                if cur.pkg not in new_pkgs:
                                    new_pkgs[cur.pkg] = cur

                                # If we've already got a version of this dependency
                                # in the mapping, it means the version we're currently
                                # looking at must be the same or we can't use this
                                # combination.
                                else:
                                    if new_pkgs[cur.pkg] != cur:

                                        # Report conflicts using a tuple of both the
                                        # conflicting installations.
                                        yield (new_pkgs[cur.pkg], cur)
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
                            deps += [dep_inst]
                            for d in self.reduce_dependencies(deps, new_pkgs, cur_index + 1):
                                yield d
                            del deps[-1]

                # There are no installations for this dependency. If it's
                # a required dependency then something has gone very wrong.
                elif cur_dep.required:

                    # Throw an error here.
                    # TODO
                    sys.exit()

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
        return True

#    def setup(self):
#        SConfig.Node.setup(self)
#        modified = []
#        if self.env['shared_libraries']:
#            for pkg in self.env.package_list:
#                if isinstance(pkg, SConfig.Package):
#                    pkg.require_shared = True
#                    modified += [pkg]
#        return modified

    def enable(self, scons_env, old_state=None):
        SConfig.Node.enable(self, scons_env, old_state)

        # Setup debugging flags.
        if self.env['with_debug']:
            scons_env.MergeFlags('-g')

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
