import os, platform, glob, copy
import SCons.Script
import SConfig

class Package(SConfig.Node):
    """This is a big one. Describes how to search for a package
    that has been installed on a system. There are a lot of
    options that can be modified to refine how the search
    proceeds and also how results are generated."""

    def __init__(self, scons_env, scons_opts, required=False, **kw):
        SConfig.Node.__init__(self, scons_env, scons_opts, required)

        # This will be set in the preprocessor.
        self.have_define         = kw.get('have_define', '')

        # Search options.
        self.base_dirs           = [] #['']
        self.base_patterns       = [] #['']
        self.sub_dirs            = [] #[[[''], ['']]]
        self.header_sub_dir      = [] #['']
        self.system_header_dirs  = []
        self.system_library_dirs = []

        # Which headers do we require?
        self.headers             = [] #[['']]

        # Library options.
        self.libraries           = [] #[['']]
        self.shared_libraries    = [] # Only libraries listed here will be considered
                                      # when checking for shared libraries.
        self.extra_libraries     = [] # These libraries are not considered locatable.
        self.require_shared      = kw.get('require_shared', False)
        self.symbols             = [([], '')] #[([''], '')]
        self.symbol_setup        = ''
        self.symbol_teardown     = ''
        self.symbol_prototypes   = [] #['']
        self.symbol_calls        = [] #['']

        # Framework options.
        self.frameworks          = [] #[['']]

        # Will be set after configuration.
        self.candidates_built = False
        self.candidates = []
        self.installations = []
        self.selected = None

        # Private stuff.
        self.platform = self.dependency(SConfig.Platform, True) # Need this so we can get
                                                                # access to information about
                                                                # the platform we're running on.
        self.checks = [self.check_candidates]

        # Set everything up.
        self.setup_search_defaults()

    def setup_options(self):
        """Two things need to happen here. The first is to tell SCons about
        any options we're going to use. Do that by calling
        'self.opts.AddOptions'. The second is to add entries to
        'self.option_map', which maps a command-line option to it's
        envrionment equivalent."""

        SConfig.Node.setup_options(self)
        self.opts.AddOptions(
            SCons.Script.BoolOption('with_' + self.command_name,
                                    'Turn on/off package %s' % self.name, 1),
            SCons.Script.PathOption(self.command_name + '_dir',
                                    '%s installation path' % self.name,
                                    None, SCons.Script.PathOption.PathIsDir),
            SCons.Script.PathOption(self.command_name + '_inc_dir',
                                    '%s header installation path' % self.name,
                                    None, SCons.Script.PathOption.PathIsDir),
            SCons.Script.PathOption(self.command_name + '_lib_dir',
                                    '%s library installation path' % self.name,
                                    None, SCons.Script.PathOption.PathIsDir),
            (self.command_name + '_lib',
             '%s libraries' % self.name,
             None, None),
            (self.command_name + '_framework',
             '%s framework' % self.name,
             None, None))
        self.option_map = {'with_' + self.command_name: None,
                           self.command_name + '_dir': self.environ_name + '_DIR'}

    def setup_search_defaults(self):
        """Setup the usual search paths for packages depending on the kind of system
        we're on."""

        if self.platform.system in ['Darwin', '*ix']:
            self.base_dirs = ['/usr', '/usr/local']
            self.sub_dirs = [[['include'], ['lib']]]
            if self.platform.bits == 64:
                self.sub_dirs = [[['include'], ['lib64']],
                                 [['include'], [os.path.join('lib', '64')]]] + self.sub_dirs

        # Set Darwin specific defaults.
        if self.platform.system == 'Darwin':
            self.base_dirs += ['/sw']

        # Set Window specific defaults.
        if self.platform.system == 'Windows':
            pass # TODO

        # Combine these guys to build default system paths. We need these to ensure specific
        # include paths are used before generic ones.
        for base_dir in ['/usr', '/sw']:
            for hdr_dirs, lib_dirs in self.combine_sub_dirs(base_dir):
                hdr_dirs = [os.path.join(base_dir, h) for h in hdr_dirs]
                self.system_header_dirs += [h for h in hdr_dirs if h not in self.system_header_dirs]
                lib_dirs = [os.path.join(base_dir, l) for l in lib_dirs]
                self.system_library_dirs += [l for l in lib_dirs if l not in self.system_library_dirs]

    def setup(self):
        """Finalise everything before running the configuration checks."""

        SConfig.Node.setup(self)
        modified = []
        if self.require_shared: # Depends on 'dl' if shared.
            if self.name != 'dl': # Make sure we're not the 'dl' package.
                modified = [self.dependency(SConfig.packages.dl)]

        # If we havn't already done so, build the list of candidates.
        if not self.candidates_built:
            self.candidates_built = True
            self.build_candidates()

        # Return a list of modified packages.
        return modified

    def build_candidates(self):
        """For package configuration, here is where we need to begin filling
        out our list of installations. We build candidates, filling them with
        as much information as possible (usually only the location at which
        they can be found and perhaps some supporting libraries that are
        required by particular installations)."""

        for inst in self.generate_installations(): # Generate each candidate
            # Scan the installation for as much information as we can. If successful,
            # add the installation to our list of candidates.
            if self.process_installation(inst):
                self.add_candidate(inst)

    def add_candidate(self, inst):
        """Add a candidate installation if it is not already there."""

        if inst not in self.candidates:
            self.candidates += [inst]

    def process_installation(self, inst):
        """This method gives us a chance to modify any of the details of this
        installation before moving on to checking it's validity. Here we can
        also determine which, if any, other dependant installations are given
        by this installation's package config (or any other means)."""

        return True

    def check_candidates(self):
        """Runs a sequence of tests to confirm the validity of candidate
        installations. All that pass are moved to the 'installations'
        member."""

        # Expand the candidates to include dependent packages.
        self.add_dependencies()

        # If we have no candidates, report the problem.
        if not self.candidates:
            self.ctx.Display('  No candidates found.\n')
            return False

        # Try out all candidates.
        cur = 1 # Keep track of which one we're currently trying.
        for inst in self.candidates:
            # Print current status.
            self.ctx.Log('  Trialing candidate %d of %d ...\n' % (cur, len(self.candidates)))
            self.ctx.Log('    ' + inst.__str__(False)[:-1].replace('\n', '\n    ') + '\n')
            print '  Trialing candidate %d of %d ...\r' % (cur, len(self.candidates)),
            cur = cur + 1

            # Check for the headers and libraries.
            if self.check_headers(inst) and self.check_libraries(inst):
                
                # If the checks passed, include the 'have_define' if there and add
                # this installation to the list.
                if self.have_define:
                    inst.add_cpp_defines(self.have_define)
                self.installations += [inst]

        # Print results.
        if len(self.installations) == 1:
            self.ctx.Display('\n  Found %d valid configuration.\n' % len(self.installations))
        else:
            self.ctx.Display('\n  Found %d valid configurations.\n' % len(self.installations))

        # Log complete results to file.
        for inst in self.installations:
            self.ctx.Log(inst.__str__(False))

        # If we couldn't find a valid installation return negative.
        if len(self.installations) == 0:
            return False
        return True

    def add_dependencies(self):
        """Copies each of the candidate installations and includes all
        permutations of dependent packages."""

        new_candidates = []
        for deps in self.combine_dependencies():
            for inst in self.candidates:
                new_inst = copy.copy(inst)
                new_inst.deps = list(deps)
                new_candidates += [new_inst]
        self.candidates = new_candidates

    def combine_dependencies(self, deps=[], pkgs={}, cur_index=0):
        """Each combination of dependent installations represents a
        unique installation of this package. This method generates sets
        of unique dependency combinations."""

        # The dictionary 'pkgs' is to map packages to installations. This is
        # needed to prevent multiple installations being used for the same
        # packages.

        if cur_index == len(self.deps):
            yield deps # Complete permutation.
        else:
            cur_dep, required = self.deps[cur_index]

            # If the dependency isn't actually required by this package, include
            # a combination that doesn't use it.
            if not required:
                for d in self.combine_dependencies(deps, pkgs, cur_index + 1):
                    yield d

            # We can only iterate over installations if we're dealing with a package.
            if isinstance(cur_dep, Package):

                # Check if we already have this package selected.
                if cur_dep in pkgs:
                    deps += [pkgs[cur_dep]]
                    for d in self.combine_dependencies(deps, pkgs, cur_index + 1):
                        yield d
                    del deps[-1]

                elif len(cur_dep.installations):
                    # Try each installation.
                    for dep_inst in cur_dep.installations:
                        # Traverse the dependency and collect any sub-dependencies
                        # into a copied package list.
                        new_pkgs = dict(pkgs)
                        rem = [dep_inst]
                        while len(rem):
                            cur = rem.pop()
                            if isinstance(cur, Installation):
                                if cur.pkg not in new_pkgs:
                                    new_pkgs[cur.pkg] = cur
                                rem += cur.deps

                        # Set dependency and recurse.
                        deps += [dep_inst]
                        for d in self.combine_dependencies(deps, new_pkgs, cur_index + 1):
                            yield d
                        del deps[-1]

                elif required:
                    # There are no installations for this dependency. If it's
                    # a required dependency then something has gone very wrong.
                    # Throw an error here.
                    # TODO
                    sys.exit()
            else:
                deps += [cur_dep]
                for d in self.combine_dependencies(deps, pkgs, cur_index + 1):
                    yield d
                del deps[-1]

    def get_check_headers_fail_reason(self, fail_logs):
        return ''

    def get_check_symbols_fail_reason(self, fail_logs):
        return ''

    def configure(self, scons_ctx):
        # If this package is deslected just return now.
        if not self.required and not self.env['with_' + self.command_name]:
            return

        # Run the configuration.
        SConfig.Node.configure(self, scons_ctx)

        # If this was a critical fail, try and help the user.
        if self.required and not self.result:
            self.ctx.Display('\nThe required package ' + self.name + ' could not be found.\n')
            if len(self.option_map.keys()):
                self.ctx.Display('You can directly specify search parameters for this package via\n')
                self.ctx.Display('the following command line options:\n\n')
                for opt in self.option_map.iterkeys():
                    self.ctx.Display('  ' + opt + '\n')
                self.ctx.Display('\nRun \'scons help\' for details on available options.\n')
            self.ctx.Display('\n')
            self.env.Exit()

    def check_headers(self, inst):
        """Determine if the required headers are available with the current construction
        environment settings."""

        # If there are no headers to check, automatically pass.
        if not self.headers:
            inst.hdrs = []
            return True

        fail_logs = []
        result = [False, '']
        hdr_list = self.make_list(self.headers)
        for inst.hdrs in hdr_list: # Check each selection of headers.
            # See if we can locate the files themselves.
            found = False
            for hdr in inst.hdrs:
                for hdr_dir in inst.hdr_dirs:
                    path = os.path.join(inst.base_dir, hdr_dir, hdr)
                    if os.path.exists(path):
                        found = True
                        break
                if not found: break
            if not found: continue

            # Try to compile the source.
            src = self.get_header_source(inst)
            old_state = {}
            inst.enable(self.env, old_state)
            result = self.compile_source(src)
            self.pop_state(old_state)

            # Check the results.
            if result[0]: return True # If we found 'em, return.
            fail_logs += [result[1]]

        # If we failed to find the headers, try and determine what happened.
        inst.hdrs = None # Clear headers so we know what failed.
        msg = self.get_check_headers_fail_reason(fail_logs)
        if not msg:
            msg = 'Headers not found.'
        inst.fail_reason = msg
        return False

    def check_libraries(self, inst):
        """Check if the currently selected location is a valid installation of the
        required package. At this stage we know that the paths given in the location
        actually exist and we need to confirm that the libraries in 'libs' exist."""

        fail_logs = []
        result = [False, '']
        for inst.libs in self.generate_libraries(inst):

            # Only continue in to this nested section if there are
            # actually some libraries to look for.
            if inst.libs:

                # See if we can locate the files themselves.
                found = False
                for lib in inst.libs:

                    # Don't check for the existance of a library if it's considered
                    # an extra library. This is because extra libraries can often exist
                    # in only compiler known locations.
                    if lib in self.extra_libraries + inst.extra_libs:
                        found = True
                        continue

                    # Check each library directory extension.
                    for lib_dir in inst.lib_dirs:
                        static_name = self.env.subst('${LIBPREFIX}' + lib + '${LIBSUFFIX}')
                        shared_name = self.env.subst('${SHLIBPREFIX}' + lib + '${SHLIBSUFFIX}')
                        static_path = os.path.join(inst.base_dir, lib_dir, static_name)
                        shared_path = os.path.join(inst.base_dir, lib_dir, shared_name)
                        if os.path.exists(static_path) or os.path.exists(shared_path):
                            found = True
                            break
                        else:
                            found = False
                    if not found: break
                if not found: continue

            # If there are no libraries and no frameworks to check for,
            # we can just return success now.
            elif not inst.fworks:
                return True

            # Check that we can link against the libraries by trying to link against
            # a particular set of symbols.
            old_state = {}
            inst.enable(self.env, old_state)
            result = self.check_symbols(inst)
            self.pop_state(old_state)

            # Handle resutls.
            if result[0]:
                if not self.require_shared:
                    return True
                break
            fail_logs += [result[1]]

        # Check if we have shared libraries.
        if result[0]:
            inst.enable(self.env, old_state)
            inst.have_shared = self.check_shared(inst)
            self.pop_state(old_state)
            if inst.have_shared:
                return True

        # Figure out what to report.
        if self.require_shared:
            inst.fail_reason = 'No shared libraries.'
        else:
            inst.libs = None # No libraries found.
            for inst.fail_reason in fail_logs:
                if inst.fail_reason: break
            if not inst.fail_reason:
                inst.fail_reason = 'Libraries not found.'
        return False

    def check_symbols(self, inst):
        """We have our paths and libraries setup, now we need to see if we can find
        one of the set of required symbols in the libraries."""

        # Keep a list of the reasons why each set of symbols failed to compile.
        fail_logs = []

        # Try out each set of symbols.
        for syms in self.symbols:

            # Try to link the source and get the results.
            result = self.link_source(self.get_check_symbols_source(inst, syms[0]))
            if result[0]:

                # In addition to the results indicating success, we need to check
                # the output for any warnings that may indicate failure.
                if result[1].find('skipping incompatible') != -1:

                    # The library path specified wasn't actually used, instead the
                    # libraries in there were incompatible and the compiler was
                    # able to make use of libraries in default search paths.
                    fail_logs += ['Incompatible libraries.']
                    result[0] = 0
                    continue

                # Succeeded, if there are any preprocessor defines required by
                # this set of symbols, add them to the instance now.
                if syms[1]:
                    inst.add_cpp_defines(syms[1]) # Add the CPP defines.

                # We're done, break out of the symbols search.
                break

            # If we failed we need to log the reason why.
            fail_logs += [result[1]]

        if not result[0]:
            reason = self.get_check_symbols_fail_reason(fail_logs)
        else:
            reason = ''
        return [result[0], reason]

    def check_shared(self, inst):
        """Confirm that there are shared versions of this package's libraries available.
        At this point we know we can link against the libraries."""

        # If we don't have any libraries to check, return True.
        if not len(inst.libs):
            return True

        # Build a binary to try and dynamically open the libraries in order.
        result = [1, '', '']
        src = self.get_header_source(inst)
        src += """
int main(int argc, char* argv[]) {
  void* lib[%d];
""" % len(inst.libs)
        for l in inst.libs:
            if self.shared_libraries and l not in self.shared_libraries:
                continue
            if l in inst.extra_libs:
                continue
            offs = ''
            for p in self.generate_library_paths(inst, l):
                offs += '  '
                if len(offs) > 2:
                    src += '{\n'
                src += '%slib[%d] = dlopen("%s", RTLD_NOW);\n' % (offs, inst.libs.index(l), p)
                src += '%sif( !lib[%d] ) ' % (offs, inst.libs.index(l))
            src += 'return 1;\n'
            while len(offs) > 2:
                offs = offs[:-2]
                src += offs + '}\n'
        src += '  return 0;\n}\n'
        if not self.run_source(src)[0]:
            self.ctx.Log('  Failed to open shared libraries.\n')
            return False
        return True

    def generate_library_paths(self, inst, library):
        lib_name = self.env.subst('${SHLIBPREFIX}' + library + '${SHLIBSUFFIX}')
        if inst.lib_dirs:
            for d in inst.lib_dirs:
                path = os.path.join(inst.base_dir, d, lib_name)
                yield os.path.abspath(path)
        else:
            yield lib_name

    def generate_installations(self):
        """Generate a set of potential package locations. Locations are of the form
        ['base_dir', ['header_dirs'], ['lib_dirs'], ['frameworks']]."""

        # If we've been given options directly specifying the location of this
        # package we need to use those in place of searching for locations.
        base_dir = self.env.get(self.command_name + '_dir', '')
        inc_dir = self.env.get(self.command_name + '_inc_dir', '')
        lib_dir = self.env.get(self.command_name + '_lib_dir', '')
        fwork = self.env.get(self.command_name + '_framework', '')
        if inc_dir or lib_dir:
            if not (inc_dir and lib_dir):
                print '   Error: must specify both of'
                print '      ' + self.command_name + '_inc_dir'
                print '      ' + self.command_name + '_lib_dir'
                env.Exit()
            yield Installation(self, '', [inc_dir], [lib_dir], [fwork])
            return
        have_loc = base_dir or inc_dir or lib_dir

        # Produce an empty location to see if the package exists in a default
        # location.
#        if not have_loc:
#            yield Installation(self)

        # Combine all possible base directories.
        if not base_dir:
            base_dirs = list(self.base_dirs)
            for dir in self.base_dirs:
                for ptrn in self.base_patterns:
                    base_dirs += glob.glob(os.path.join(dir, ptrn))
        else:
            base_dirs = [base_dir]

        # Make sure there are no symbolic links in the base directories. If there are
        # any in there, expand them and make sure there are no duplicates.
        new_base_dirs = []
        for d in base_dirs:
            d = os.path.realpath(d)
            if d not in new_base_dirs:
                new_base_dirs += [d]
        base_dirs = new_base_dirs

        # Add an empty set to the beginning of the frameworks.
        if not fwork:
            frameworks = [[]] + self.frameworks
        else:
            frameworks = [fwork]

        for fw in frameworks:
            # If we have some frameworks, try them alone first.
            if not have_loc and fw:
                yield Installation(self, '', [], [], fw)

        # Traverse the list of base directories and form each
        # installation.
        for dir in base_dirs:
            for inst in self.combine_base_dir(dir):
                yield inst

    def combine_base_dir(self, base_dir):
        """From a given base directory, combine possible header and library
        extensions to yield a set of installations."""

        if os.path.exists(base_dir) and os.path.isdir(base_dir):
            for hdr, lib in self.combine_sub_dirs(base_dir):
                yield Installation(self, base_dir, hdr, lib)
                for sub in self.combine_header_sub_dir(base_dir, hdr):
                    yield Installation(self, base_dir, sub, lib)

    def combine_sub_dirs(self, base_dir):
        """Take a base directory and combine it with the set of header and library
        subdirectories. Yields (['header_dirs'], ['lib_dirs'])."""

        for hdr, lib in self.sub_dirs:
            loc_okay = True
            hdr_dirs = []
            lib_dirs = []

            # Combine header subdirectories.
            for h in hdr:
                dir = self.join_sub_dir(base_dir, h)
                if not (os.path.exists(dir) and os.path.isdir(dir)):
                    loc_okay = False
                    break
                hdr_dirs += [h]
            if not loc_okay:
                continue

            # Combine library subdirectories.
            for l in lib:
                dir = self.join_sub_dir(base_dir, l)
                if not (os.path.exists(dir) and os.path.isdir(dir)):
                    loc_okay = False
                    break
                lib_dirs += [l]
            if not loc_okay:
                continue

            yield (hdr_dirs, lib_dirs)

    def combine_header_sub_dir(self, base_dir, hdr_dirs):
        if not self.header_sub_dir or not hdr_dirs:
            return
        for sub_dir in self.header_sub_dir:
            cand = [os.path.join(h, sub_dir) for h in hdr_dirs if h]
            for d in cand:
                path = os.path.join(base_dir, d)
                if not (os.path.exists(path) and os.path.isdir(path)):
                    return
            yield cand

    def join_sub_dir(self, base_dir, sub_dir):
        if os.path.isabs(sub_dir):
            return sub_dir
        return os.path.join(base_dir, sub_dir)

    def generate_libraries(self, inst):
        if inst.fworks: # Try any frameworks by themselves first.
            yield inst.extra_libs
        for libs in self.libraries:
            yield libs + inst.extra_libs

    def enable(self, scons_env, old_state=None):
        SConfig.Node.enable(self, scons_env, old_state)
        if self.selected:
            self.selected.enable(scons_env, old_state)

    def get_all_headers(self, inst, headers):
        if not self.result:
            return
        for d in inst.deps:
            if isinstance(d, Installation):
                d.pkg.get_all_headers(d, headers)
        headers += [h for h in inst.hdrs if h not in headers]

    def get_header_source(self, inst):
        src = '#include<stdlib.h>\n#include<stdio.h>\n#include<string.h>\n'
        hdrs = []
        for d in inst.deps:
            if isinstance(d, Installation):
                d.pkg.get_all_headers(d, hdrs)
        hdrs += inst.hdrs
        for h in hdrs:
            src += '#include<' + h + '>\n'
        return src

    def get_check_version_source(self, inst):
        """Generates the source code required to check the package version."""

        src = self.get_header_source(inst)
        return src

    def get_check_symbols_source(self, inst, symbols):
        src = self.get_header_source(inst)
        for sym, proto in zip(symbols, self.symbol_prototypes):
            src += (proto % sym) + '\n'
        src += 'int main(int argc, char* argv[]) {\n'
        if self.symbol_setup:
            src += self.symbol_setup + '\n'
        for sym, call in zip(symbols, self.symbol_calls):
            src += (call % sym) + '\n'
        if self.symbol_teardown:
            src += self.symbol_teardown + '\n'
        src += 'return 0;\n}\n'
        return src

    def make_list(self, var):
        """Convert anything into a list. Handles things that are already lists,
        tuples and strings."""

        if isinstance(var, str):
            return [var]
        elif isinstance(var, (list, tuple)) and not var:
            return []
        else:
            return list(var)

class Installation:
    def __init__(self, package, base_dir='', hdr_dirs=[], lib_dirs=[], fworks=[]):
        self.pkg = package
        self.deps = []
        self.base_dir = base_dir
        self.hdr_dirs = list(hdr_dirs)
        self.lib_dirs = list(lib_dirs)
        self.hdrs = []
        self.libs = []
        self.extra_libs = []
        self.have_shared = None
        self.fworks = list(fworks)
        self.cpp_defines = []
        self.fail_reason = ''

    def __eq__(self, inst):
        return (self.pkg is inst.pkg and
                self.deps == inst.deps and
                self.base_dir == inst.base_dir and
                self.hdr_dirs == inst.hdr_dirs and self.lib_dirs == inst.lib_dirs and
                self.fworks == inst.fworks)

    def add_hdr_dirs(self, hdr_dirs):
        dir_list = self.pkg.make_list(hdr_dirs)
        self.hdr_dirs += [d for d in dir_list if d not in self.hdr_dirs]

    def add_lib_dirs(self, lib_dirs):
        dir_list = self.pkg.make_list(lib_dirs)
        for d in dir_list:
            d = os.path.normpath(d)
            if d not in self.lib_dirs:
                self.lib_dirs += [d]

    def add_libs(self, libs):
        lib_list = self.pkg.make_list(libs)
        self.libs += [l for l in lib_list if l not in self.libs]

    def add_extra_libs(self, libs):
        lib_list = self.pkg.make_list(libs)
        self.extra_libs += [l for l in lib_list if l not in self.extra_libs]

    def add_cpp_defines(self, defs):
        def_list = self.pkg.make_list(defs)
        self.cpp_defines += [d for d in def_list if d not in self.cpp_defines]

    def enable(self, scons_env, old_state=None):
        """Inserts this installation's information into the provided SCons
        environment."""

        # Insert dependencies.
        self.enable_dependencies(scons_env, old_state)

        if self.cpp_defines:
            self.pkg.backup_variable(scons_env, 'CPPDEFINES', old_state)
            scons_env.AppendUnique(CPPDEFINES=self.cpp_defines)

        if self.hdr_dirs:
            self.pkg.backup_variable(scons_env, 'CPPPATH', old_state)
            rev_hdr_dirs = list(self.hdr_dirs)
            rev_hdr_dirs.reverse()
            for d in rev_hdr_dirs:
                full_dir = self.pkg.join_sub_dir(self.base_dir, d)
                if full_dir in self.pkg.system_header_dirs:
                    if not os.path.isabs(full_dir): full_dir = '#' + full_dir
                    scons_env.AppendUnique(CPPPATH=[full_dir])
                else:
                    if not os.path.isabs(full_dir): full_dir = '#' + full_dir
                    scons_env.PrependUnique(CPPPATH=[full_dir])

        if self.fworks:
            self.pkg.backup_variable(scons_env, 'FRAMEWORKS', old_state)
            self.pkg.backup_variable(scons_env, 'CPPPATH', old_state)
            scons_env.PrependUnique(FRAMEWORKS=self.fworks)
            rev_fworks = list(self.fworks)
            rev_fworks.reverse()
            for fw in rev_fworks: # Sort of a hack for Mac OS X.
                path = '/System/Library/Frameworks/' + fw + '.framework/Headers'
                scons_env.PrependUnique(CPPPATH=[path])

        if self.lib_dirs:
            self.pkg.backup_variable(scons_env, 'LIBPATH', old_state)
            self.pkg.backup_variable(scons_env, 'RPATH', old_state)
            rev_lib_dirs = list(self.lib_dirs)
            rev_lib_dirs.reverse()
            for d in rev_lib_dirs:
                full_dir = self.pkg.join_sub_dir(self.base_dir, d)
                abs_dir = os.path.abspath(full_dir)
                if full_dir in self.pkg.system_library_dirs:
                    if not os.path.isabs(full_dir): full_dir = '#' + full_dir
                    scons_env.AppendUnique(LIBPATH=[full_dir])
                    scons_env.AppendUnique(RPATH=[abs_dir])
                else:
                    if not os.path.isabs(full_dir): full_dir = '#' + full_dir
                    scons_env.PrependUnique(LIBPATH=[full_dir])
                    scons_env.PrependUnique(RPATH=[abs_dir])

        if self.libs:
            self.pkg.backup_variable(scons_env, 'LIBS', old_state)
            scons_env.PrependUnique(LIBS=self.libs)

    def enable_dependencies(self, scons_env, old_state={}):
        """Enables all available dependencies."""

        for d in self.deps:
            d.enable(scons_env, old_state)

    def __str__(self, brief=True):
        """Convert to printable string."""

        txt = 'Package: %s\n' % self.pkg.name
        dep_txt = ''
        for dep in self.deps:
            if isinstance(dep, Installation):
                if brief:
                    dep_txt += '    %s\n' % dep.pkg.name
                else:
                    dep_txt += '    ' + str(dep)[:-1].replace('\n', '\n    ') + '\n'
        if dep_txt:
            txt += '  Dependencies:\n'
            txt += dep_txt
        txt += '  Base directory: %s\n' % self.base_dir
        txt += '  Header extensions: %s\n' % self.hdr_dirs
        txt += '  Library extensions: %s\n' % self.lib_dirs
        if self.libs:
            txt += '  Libraries: %s\n' % self.libs
        if self.fworks:
            txt += '  Frameworks: %s\n' % self.fworks
        if self.cpp_defines:
            txt += '  Exporting: %s\n' % self.cpp_defines
        return txt
