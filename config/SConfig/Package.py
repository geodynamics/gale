import os, platform, glob
import SCons.Script
import SConfig

class Package(SConfig.Node):
    def __init__(self, scons_env, scons_opts, required=False, **kw):
        SConfig.Node.__init__(self, scons_env, scons_opts, required)

        # This will be set in the preprocessor.
        self.have_define         = kw.get('have_define', '')

        # Search options.
        self.base_dirs           = [] #['']
        self.base_patterns       = [] #['']
        self.sub_dirs            = [] #[[[''], ['']]]
        self.header_sub_dir      = ''
        self.system_header_dirs  = []
        self.system_library_dirs = []

        # Which headers do we require?
        self.headers             = [[]] #[['']]

        # Library options.
        self.libraries           = [] #[['']]
        self.extra_libraries     = [] #['']
        self.shared_libraries    = [] # Only libraries listed here will be considered
                                      # when checking for shared libraries.
        self.require_shared      = False
        self.symbols             = [([], '')] #[([''], '')]
        self.symbol_setup        = ''
        self.symbol_teardown     = ''
        self.symbol_prototypes   = [] #['']
        self.symbol_calls        = [] #['']

        # Framework options.
        self.frameworks          = [] #[['']]

        # Will be set after configuration.
        self.base_dir = ''
        self.hdr_dirs = []
        self.lib_dirs = []
        self.hdrs = []
        self.libs = []
        self.have_shared = None
        self.fworks = []
        self.cpp_defines = []

        # Private stuff.
        self.platform = self.dependency(SConfig.Platform, True) # Need this so we can get
                                                                # access to information about
                                                                # the platform we're running on.
        self.checks = [self.find_package] # The basic package location check.

        # Set everything up.
        self.setup_search_defaults()

    def setup_options(self):
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

    def get_check_headers_fail_reason(self, fail_logs):
        return ''

    def get_check_symbols_fail_reason(self, fail_logs):
        return ''

    def configure(self, scons_ctx):
        if not self.required and not self.env['with_' + self.command_name]:
            return
        SConfig.Node.configure(self, scons_ctx)

    def setup(self):
        SConfig.Node.setup(self)
        if self.require_shared:
            self.dependency(SConfig.packages.dl)

    def find_package(self):
        """Basic check routine for locating the package."""
        result = False
        self.ctx.Display('  Searching locations:\n')
        for loc in self.generate_locations():
            self.ctx.Display('    %s\n' % str(loc))
            result = self.check_location(loc)
            if result:
                if self.have_define:
                    self.cpp_defines += [self.have_define]
                self.base_dir = loc[0]
                self.hdr_dirs = loc[1]
                self.lib_dirs = loc[2]
                self.fworks = loc[3]
                break
        return result

    def check_location(self, location):
        """Check if the currently selected location is a valid installation of the
        required package. At this stage we know that the paths given in the location
        actually exist."""
        old_state = self.enable_location_state(location)

        # Check for the headers.
        if not self.check_headers(location):
            self.pop_state(old_state)
            return False

        # Scan each set of libraries in turn.
        if not self.check_libraries(location):
            self.pop_state(old_state)
            return False

        self.pop_state(old_state)
        return True

    def check_headers(self, location):
        """Determine if the required headers are available with the current construction
        environment settings."""
        fail_logs = []
        for hdrs in self.headers:
            src = self.get_header_source(hdrs)
            result = self.compile_source(src)
            if result[0]:
                self.hdrs = list(hdrs)
                break
            fail_logs += [result[1]]
        if not result[0]:
            msg = self.get_check_headers_fail_reason(fail_logs)
            if not msg:
                msg = 'Headers not found.'
            self.ctx.Display('      ' + msg + '\n')
        return result[0]

    def check_libraries(self, location):
        """Check if the currently selected location is a valid installation of the
        required package. At this stage we know that the paths given in the location
        actually exist and we need to confirm that the libraries in 'libs' exist."""
        fail_reasons = []
        no_shared = False
        for libs in self.generate_libraries(location):
            old_state = self.enable_library_state(location, libs)

            # Check that we can link against the libraries by trying to link against
            # a particular set of symbols.
            result = self.check_symbols(location, libs)
            if not result[0]:
                fail_reasons += [result[1]]
                self.pop_state(old_state)
                continue

            # Check if we have shared libraries.
            if not self.require_shared or not libs:
                self.libs = list(libs)
                self.pop_state(old_state)
                return True
            elif self.check_shared(location, libs):
                self.libs = list(libs)
                self.pop_state(old_state)
                return True
            else:
                no_shared = True

        # Figure out what to report.
        if no_shared:
            reason = 'No shared libraries.'
        else:
            reason = ''
            for reason in fail_reasons:
                if reason:
                    break
            if not reason:
                reason = 'Libraries not found.'
        self.ctx.Display('      ' + reason + '\n')

        self.pop_state(old_state)
        return False

    def check_symbols(self, location, libraries):
        """We have our paths and libraries setup, now we need to see if we can find
        one of the set of required symbols in the libraries."""
        fail_logs = []
        for syms in self.symbols:
            result = self.link_source(self.get_check_symbols_source(syms[0]))
            if result[0]:
                if syms[1]:
                    self.cpp_defines += [syms[1]] # Add the CPP defines.
                break
            fail_logs += [result[1]]
        if not result[0]:
            reason = self.get_check_symbols_fail_reason(fail_logs)
        else:
            reason = ''
        return [result[0], reason]

    def check_shared(self, location, libraries):
        """Confirm that there are shared versions of this package's libraries available.
        At this point we know we can link against the libraries."""
        # Build a binary to try and dynamically open the libraries in order.
        result = [1, '', '']
        src = self.get_header_source()
        src += """
int main(int argc, char* argv[]) {
  void* lib[%d];
""" % len(libraries)
        for l in libraries:
            if self.shared_libraries and l not in self.shared_libraries:
                continue
            if l in self.extra_libraries:
                continue
            offs = ''
            for p in self.generate_library_paths(location, l):
                offs += '  '
                if len(offs) > 2:
                    src += '{\n'
                src += '%slib[%d] = dlopen("%s", RTLD_NOW);\n' % (offs, libraries.index(l), p)
                src += '%sif( !lib[%d] ) ' % (offs, libraries.index(l))
            src += 'return 1;\n'
            while len(offs) > 2:
                offs = offs[:-2]
                src += offs + '}\n'
        src += '  return 0;\n}\n'
        if not self.run_source(src)[0]:
            return False
        return True

    def generate_library_paths(self, location, library):
        lib_name = self.env.subst('${SHLIBPREFIX}' + library + '${SHLIBSUFFIX}')
        if location[2]:
            for d in location[2]:
                path = os.path.join(location[0], d, lib_name)
                yield os.path.abspath(path)
        else:
            yield lib_name

    def generate_locations(self):
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
            yield ['', [inc_dir], [lib_dir], [fwork]]
            return
        have_loc = base_dir or inc_dir or lib_dir

        # Produce an empty location to see if the package exists in a default
        # location.
        if not have_loc:
            yield ['', [], [], []]

        # Combine all possible base directories.
        if not base_dir:
            base_dirs = list(self.base_dirs)
            for dir in self.base_dirs:
                for ptrn in self.base_patterns:
                    base_dirs += glob.glob(os.path.join(dir, ptrn))
        else:
            base_dirs = [base_dir]

        # Add an empty set to the beginning of the frameworks.
        if not fwork:
            frameworks = [[]] + self.frameworks
        else:
            frameworks = [fwork]

        for fw in frameworks:
            # If we have some frameworks, try them alone first.
            if not have_loc and fw:
                yield ['', [], [], list(fw)]

            # Traverse the list of base directories, using them only if they exist.
            for dir in base_dirs:
                if os.path.exists(dir) and os.path.isdir(dir):
                    for hdr, lib in self.combine_sub_dirs(dir):
                        yield [dir, list(hdr), list(lib), list(fw)]
                        for sub in self.combine_header_sub_dir(dir, hdr):
                            yield [dir, list(sub), list(lib), list(fw)]

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
        cand = [os.path.join(h, self.header_sub_dir) for h in hdr_dirs if h]
        for d in cand:
            path = os.path.join(base_dir, d)
            if not (os.path.exists(path) and os.path.isdir(path)):
                return
        yield cand

    def join_sub_dir(self, base_dir, sub_dir):
        if os.path.isabs(sub_dir):
            return sub_dir
        return os.path.join(base_dir, sub_dir)

    def generate_libraries(self, location):
        if location[3]: # Try any frameworks by themselves first.
            yield self.extra_libraries
        for libs in self.libraries:
            yield libs + self.extra_libraries

    def enable_location_state(self, location):
        """Modify our environment to include search paths for the current location."""
        old_state = {}
        if location[1]:
            old_state['CPPPATH'] = self.env.get('CPPPATH', [])
            self.env.PrependUnique(CPPPATH=[self.join_sub_dir(location[0], l) for l in location[1]])
        if location[2]:
            old_state['LIBPATH'] = self.env.get('LIBPATH', [])
            old_state['RPATH'] = self.env.get('RPATH', [])
            lib_paths = [self.join_sub_dir(location[0], l) for l in location[2]]
            self.env.PrependUnique(LIBPATH=lib_paths)
            self.env.PrependUnique(RPATH=[os.path.abspath(p) for p in lib_paths])
        if location[3]:
            old_state['FRAMEWORKS'] = self.env.get('FRAMEWORKS', [])
            self.env.PrependUnique(FRAMEWORKS=location[3])
            if 'CPPPATH' not in old_state:
                old_state['CPPPATH'] = self.env.get('CPPPATH', [])
            for fw in location[3]: # Sort of a hack for Mac OS X.
                path = '/System/Library/Frameworks/' + fw + '.framework/Headers'
                self.env.PrependUnique(CPPPATH=[path])
        return old_state

    def enable_library_state(self, location, libs):
        """Take the current location and libraries and convert them into an SCons
        construction environment state dictionary."""
        old_state = {}
        if libs:
            old_state['LIBS'] = self.env.get('LIBS', [])
            self.env.PrependUnique(LIBS=libs)
        return old_state

    def enable(self, scons_env, old_state=None):
        SConfig.Node.enable(self, scons_env, old_state)
        if self.cpp_defines:
            self.backup_variable(scons_env, 'CPPDEFINES', old_state)
            scons_env.AppendUnique(CPPDEFINES=self.cpp_defines)

        if self.hdr_dirs:
            self.backup_variable(scons_env, 'CPPPATH', old_state)
            for d in self.hdr_dirs:
                full_dir = self.join_sub_dir(self.base_dir, d)
                if full_dir in self.system_header_dirs:
                    if not os.path.isabs(full_dir): full_dir = '#' + full_dir
                    scons_env.AppendUnique(CPPPATH=[full_dir])
                else:
                    if not os.path.isabs(full_dir): full_dir = '#' + full_dir
                    scons_env.PrependUnique(CPPPATH=[full_dir])

        if self.fworks:
            self.backup_variable(scons_env, 'FRAMEWORKS', old_state)
            self.backup_variable(scons_env, 'CPPPATH', old_state)
            scons_env.PrependUnique(FRAMEWORKS=self.fworks)
            for fw in self.fworks: # Sort of a hack for Mac OS X.
                path = '/System/Library/Frameworks/' + fw + '.framework/Headers'
                self.env.PrependUnique(CPPPATH=[path])

        if self.lib_dirs:
            self.backup_variable(scons_env, 'LIBPATH', old_state)
            self.backup_variable(scons_env, 'RPATH', old_state)
            for d in self.lib_dirs:
                full_dir = self.join_sub_dir(self.base_dir, d)
                abs_dir = os.path.abspath(full_dir)
                if full_dir in self.system_library_dirs:
                    if not os.path.isabs(full_dir): full_dir = '#' + full_dir
                    scons_env.AppendUnique(LIBPATH=[full_dir])
                    scons_env.AppendUnique(RPATH=[abs_dir])
                else:
                    if not os.path.isabs(full_dir): full_dir = '#' + full_dir
                    scons_env.PrependUnique(LIBPATH=[full_dir])
                    scons_env.PrependUnique(RPATH=[abs_dir])

        if self.libs:
            self.backup_variable(scons_env, 'LIBS', old_state)
            scons_env.PrependUnique(LIBS=self.libs)

    def get_all_headers(self, headers):
        if not self.result:
            return
        for d, r in self.deps:
            if hasattr(d, 'get_all_headers'):
                d.get_all_headers(headers)
        headers += [h for h in self.hdrs if h not in headers]

    def get_header_source(self, headers=None):
        src = '#include<stdlib.h>\n#include<stdio.h>\n#include<string.h>\n'
        hdrs = []
        for d, r in self.deps:
            if hasattr(d, 'get_all_headers'):
                d.get_all_headers(hdrs)
        if headers is None:
            hdrs += list(self.hdrs)
        else:
            hdrs += list(headers)
        for h in hdrs:
            src += '#include<' + h + '>\n'
        return src

    def get_check_symbols_source(self, symbols):
        src = self.get_header_source()
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
