import os, platform, glob
import SCons.Script
import SConfig

class Package(object):
    def __init__(self, env, options=None, required=False):
        # Construct this package's name.
        self.name = self.__module__.split('.')[-1]
        self.option_name = self.name.lower()
        self.environ_name = self.name.upper()
        self.command_options = {}
        self.environ_options = {}

        # Setup some system specific information.
        self.system = platform.system()
        if platform.platform().find('x86_64') != -1 or \
               platform.platform().find('ppc64') != -1 or \
               platform.architecture()[0].find('64') != -1:
            self.bits = 64
        else:
            self.bits = 32

        # Setup the options immediately, that way we can access them
        # in sub-classes.
        self.opts = options
        self.setup_options()
        if self.opts:
            self.opts.Update(env)

        # Language options.
        self.language          = 'C' # 'C' | 'CXX' | 'F77'
        self.have_define       = ''

        # Search options.
        self.base_dirs         = [] #['']
        self.base_patterns     = [] #['']
        self.sub_dirs          = [] #[[[''], ['']]]
        self.header_sub_dir    = ''

        # Header options.
        self.headers            = [[]] #[['']]

        # Library options.
        self.libraries         = [] #[['']]
        self.require_shared    = False
        self.use_rpath         = True
        self.symbols           = [([], '')] #[([''], '')]
        self.symbol_setup      = ''
        self.symbol_teardown   = ''
        self.symbol_prototypes = [] #['']
        self.symbol_calls      = [] #['']

        # Framework options.
        self.frameworks        = [] #[['']]

        # Version checking options.
        self.version           = [] #[['', 0]]

        # A list of checks to perform for this package.
        self.checks            = [self.find_package]

        # Once configured, these values will be set.
        self.configured = False
        self.result = [0, '', '']
        self.cpp_defines = []
        self.base_dir = ''
        self.hdr_dirs = []
        self.hdrs = []
        self.lib_dirs = []
        self.libs = []
        self.have_shared = None
        self.fworks = []
        self.state = {}

        # Private stuff.
        self.env = env
        self.deps = []
        self.required = required

        self.setup_search_defaults()

    def setup_search_defaults(self):
        # Set common defaults of Darwin and Linux.
        if self.system in ['Darwin', 'Linux']:
            self.base_dirs = ['/usr', '/usr/local']
            self.sub_dirs = [[['include'], ['lib']]]
            if self.bits == 64:
                self.sub_dirs = [[['include'], ['lib64']],
                                 [['include'], [os.path.join('lib', '64')]]] + self.sub_dirs

        # Set Darwin specific defaults.
        if self.system == 'Darwin':
            self.base_dirs += ['/sw']

        # Set Window specific defaults.
        if self.system == 'Windows':
            pass # TODO

    def setup_options(self):
        if not self.opts:
            return
        self.command_options = ['with_' + self.option_name,
                                self.option_name + 'Dir',
                                self.option_name + 'IncDir',
                                self.option_name + 'LibDir',
                                self.option_name + 'Lib',
                                self.option_name + 'Framework']
        self.environ_options = {self.option_name + 'Dir': self.environ_name + '_DIR',
                                self.option_name + 'IncDir': self.environ_name + '_INC_DIR',
                                self.option_name + 'LibDir': self.environ_name + '_LIB_DIR',
                                self.option_name + 'Lib': self.environ_name + '_LIB',
                                self.option_name + 'Framework': self.environ_name + '_FRAMEWORK'}
        self.opts.AddOptions(
            SCons.Script.BoolOption('with_' + self.option_name, 'Turn on/off %s' % self.name, 1),
            SCons.Script.PathOption(self.option_name + 'Dir',
                                    '%s installation path' % self.name,
                                    None, SCons.Script.PathOption.PathIsDir),
            SCons.Script.PathOption(self.option_name + 'IncDir',
                                    '%s header installation path' % self.name,
                                    None, SCons.Script.PathOption.PathIsDir),
            SCons.Script.PathOption(self.option_name + 'LibDir',
                                    '%s library installation path' % self.name,
                                    None, SCons.Script.PathOption.PathIsDir),
            (self.option_name + 'Lib',
             '%s libraries' % self.name,
             None, None),
            (self.option_name + 'Framework',
             '%s framework' % self.name,
             None, None))

    def get_headers_error_message(self, console):
        return ''

    def get_run_error_message(self, console):
        return ''

    def configure(self, configure_context):
        """Perform the configuration of this package. Override this method to perform
        custom configuration checks."""
        # If we've already configured this package, or it is deselected, just return.
        if self.configured or not (self.required or self.env['with_' + self.option_name]):
            return

        # Setup our configure context and environment.
        self.ctx = configure_context
        self.configured = True

        # Process dependencies first.
        result = self.process_dependencies()
        if not result[0]:
            self.ctx.Display('  ' + result[2])
            self.result = result
            return result

        # Process options.
        self.process_options()

        # Perfrom actual configuration.
        self.ctx.Message('Checking for package %s ... ' % self.name)
        self.ctx.Display('\n')
        for check in self.checks:
            result = check()
            self.result = result
            if not self.result[0]:
                break

        # Required?
        if self.required and not self.result[0]:
            self.ctx.Display('\nThe required package ' + self.name + ' could not be found.\n')
            self.ctx.Display('The printouts above should provide some information on what went wrong,\n')
            self.ctx.Display('To see further details, please read the \'config.log\' file.\n')
            if len(self.command_options):
                self.ctx.Display('You can directly specify search parameters for this package via\n')
                self.ctx.Display('the following command line options:\n\n')
                for opt in self.command_options:
                    self.ctx.Display('  ' + opt + '\n')
                self.ctx.Display('\nRun \'scons help\' for more details on these options.\n\n')
            self.env.Exit()

        # If we succeeded, store the resulting environment.
        if self.result[0]:
            if self.have_define:
                self.cpp_defines += [self.have_define]
            self.build_state()
            self.push_state(self.state)

        # Display results.
        if not self.result[0] and self.result[2]:
            self.ctx.Display('   ' + self.result[2])
        self.ctx.Display('   ')
        self.ctx.Result(self.result[0])

        return self.result

    def dependency(self, package_module, required=True):
        if self.configured:
            print 'Error: Cannot add a dependency during configuration.'
            self.env.Exit()
        pkg = self.env.Package(package_module, self.opts)
        if pkg not in [d[0] for d in self.deps]:
            self.deps += [(pkg, required)]
        return pkg

    def process_dependencies(self):
        """Ensure all dependencies have been configured before this package."""
        for pkg, req in self.deps:
            pkg.configure(self.ctx)
            if req and not pkg.result[0]:
                return [0, '', 'Missing dependency: ' + pkg.name]
        return [1, '', '']

    def process_options(self):
        """Do any initial option processing, including importing any values from
        the environment and validating that all options are consistent."""
        cmd_opts = False
        for opt in self.command_options:
            if opt in self.opts.args:
                cmd_opts = True
                break
        if cmd_opts:
            return
        for cmd, env in self.environ_options.iteritems():
            if cmd not in self.opts.args and env in self.env['ENV']:
                self.env[cmd] = self.env['ENV'][env]

    def find_package(self):
        # Search for package locations.
        self.ctx.Display('   Searching locations:\n')
        for loc in self.generate_locations():
            result = self.check_location(loc)
            self.ctx.Display('      %s\n' % str(loc))
            if result[0]: # If we succeeded, back out here.
                break
            if result[2]: # Display an error message.
                self.ctx.Display('         %s\n' % result[2])

        result = [result[0], '', '']
        return result

    def check_location(self, location):
        """Check if the currently selected location is a valid installation of the
        required package. At this stage we know that the paths given in the location
        actually exist."""
        # Validate the headers, first.
        result = self.validate_location(location)
        if not result[0]:
            return result

        # Construct our path state.
        path_state = self.build_header_state(location)
        old = self.push_state(path_state)

        # Check for the headers.
        result = self.check_headers(location)
        if not result[0]:
            self.pop_state(old)
            return result

        # Scan each set of libraries in turn.
        libs = []
        for libs in self.generate_libraries(location):
            result = self.check_libs(location, libs)
            if result[0]:
                break

        # Store last known configuration.
        self.store_result(result, location, libs)

        # Roll-back on state.
        self.pop_state(old)
        return result

    def check_headers(self, location):
        """Determine if the required headers are available with the current construction
        environment settings."""
        for hdrs in self.headers:
            src = self.get_header_source(hdrs)
            result = self.run_scons_cmd(self.ctx.TryCompile, src, '.c')
            if result[0]:
                self.hdrs = list(hdrs)
                break
        if not result[0]:
            msg = self.get_headers_error_message(result[1])
            if not msg:
                msg = 'Failed to locate headers.'
        else:
            msg = ''
        return [result[0], '', msg]

    def check_libs(self, location, libs):
        """Check if the currently selected location is a valid installation of the
        required package. At this stage we know that the paths given in the location
        actually exist and we need to confirm that the libraries in 'libs' exist."""
        # Validate the libraries.
        result = self.validate_libraries(location, libs)
        if not result[0]:
            return result

        # Construct the library state.
        lib_state = self.build_lib_state(location, libs)
        old = self.push_state(lib_state)

        # Check that we can link against the libraries by trying to link against
        # a particular set of symbols.
        result = self.check_symbols()
        if not result[0]:
            if not result[2]:
                result[2] = 'Failed to link against library(s).'
            self.pop_state(old)
            return result

        # Check if we have shared libraries.
        if self.require_shared:
            result = self.check_shared(location, libs)
            if not result[0] and not result[2]:
                result[2] = 'No shared library(s) available.'

        # Roll-back on our state.
        self.pop_state(old)
        return result

    def check_symbols(self):
        """We have our paths and libraries setup, now we need to see if we can find
        one of the set of required symbols in the libraries."""
        for syms in self.symbols:
            result = self.run_source(self.get_check_symbols_source(syms[0]))
            if result[0]:
                if syms[1]:
                    self.cpp_defines += [syms[1]] # Add the CPP defines.
                break
        return result

    def check_shared(self, location, libs):
        """Confirm that there are shared versions of this package's libraries available."""
        if not self.pkg_dl.result[0]:
            return [0, '', 'No dynamic loader found (libdl).']

        # Build a binary to try and dynamically open the library.
        result = [1, '', '']
        for l in libs:
            src = self.get_header_source()
            src += """
int main(int argc, char* argv[]) {
  void* lib;
  lib = dlopen("%s", RTLD_LAZY);
  return lib ? 0 : 1;
}
""" % self.env.subst('${SHLIBPREFIX}' + l + '${SHLIBSUFFIX}')
            result = self.run_source(src)
            if not result[0]:
                break
        return result

    def run_source(self, source):
        """At this point we know all our construction environment has been set up,
        so we should be able to build and run the application."""
        result = self.run_scons_cmd(self.ctx.TryRun, source, '.c')
        msg = self.get_run_error_message(result[1])
        return [result[0][0], result[0][1], msg]

    def validate_location(self, location):
        """Confirm that the location is okay, possibly modifying it in place if
        there are additional locations to search."""
        return [1, '', '']

    def validate_libraries(self, location, libs):
        """Confirm that the specified libraries are okay, possibly modifying in
        place the list of libraries."""
        return [1, '', '']

    def generate_locations(self):
        """Generate a set of potential package locations. Locations are of the form
        ['base_dir', ['header_dirs'], ['lib_dirs'], ['frameworks']]."""
        # If we've been given options directly specifying the location of this
        # package we need to use those in place of searching for locations.
        base_dir = self.env.get(self.option_name + 'Dir', '')
        inc_dir = self.env.get(self.option_name + 'IncDir', '')
        lib_dir = self.env.get(self.option_name + 'LibDir', '')
        fwork = self.env.get(self.option_name + 'Framework', '')
        if inc_dir or lib_dir:
            if not (inc_dir and lib_dir):
                print '   Error: must specify both of'
                print '      ' + self.option_name + 'IncDir'
                print '      ' + self.option_name + 'LibDir'
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

    def generate_libraries(self, location):
        for libs in self.libraries:
            yield libs

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

    def build_header_state(self, location):
        """Build a construction state for including headers."""
        state = {}
        if location[1]:
            state['CPPPATH'] = [self.join_sub_dir(location[0], l) for l in location[1]]
        if location[3]:
            state['FRAMEWORKS'] = location[3]
        return state

    def build_lib_state(self, location, libs):
        """Take the current location and libraries and convert them into an SCons
        construction environment state dictionary."""
        state = {}
        if location[2]:
            state['LIBPATH'] = [self.join_sub_dir(location[0], l) for l in location[2]]
            if self.use_rpath:
                state['RPATH'] = [os.path.abspath(p) for p in state['LIBPATH']]
        if location[3]:
            state['FRAMEWORKS'] = location[3]
        if libs:
            state['LIBS'] = libs
        return state

    def build_state(self):
        self.state = {}
        if self.cpp_defines:
            self.state['CPPDEFINES'] = self.cpp_defines
        if self.hdr_dirs:
            self.state['CPPPATH'] = [self.join_sub_dir(self.base_dir, d) \
                                         for d in self.hdr_dirs]
        if self.fworks:
            self.state['FRAMEWORKS'] = self.fworks
        if self.lib_dirs:
            self.state['LIBPATH'] = [self.join_sub_dir(self.base_dir, d) \
                                         for d in self.lib_dirs]
            if self.use_rpath:
                self.state['RPATH'] = [os.path.abspath(p) for p in self.state['LIBPATH']]
        if self.libs:
            self.state['LIBS'] = self.libs

    def store_result(self, result, location, libs):
        self.result = result
        self.base_dir = location[0]
        self.hdr_dirs = location[1]
        self.lib_dirs = location[2]
        self.libs = libs
        if self.require_shared:
            self.have_shared = True
        else:
            self.have_shared = None
        self.fworks = location[3]

    def push_state(self, state, append=False):
        old = {}
        copy = dict(state)
        for k, v in copy.iteritems():
            if not v:
                continue
            if not isinstance(v, list):
                copy[k] = [v]
            else:
                copy[k] = v
            old[k] = self.env.get(k, [])
        if append:
            self.env.AppendUnique(**copy)
        else:
            self.env.PrependUnique(**copy)
        return old

    def pop_state(self, old):
        self.env.Replace(**old)

    def get_all_headers(self, headers):
        if not self.result[0]:
            return
        headers += [h for h in self.hdrs if h not in headers]
        for d, r in self.deps:
            d.get_all_headers(headers)

    def get_header_source(self, headers=None):
        src = '#include<stdlib.h>\n#include<stdio.h>\n#include<string.h>\n'
        if headers is None:
            hdrs = list(self.hdrs)
        else:
            hdrs = list(headers)
        for d, r in self.deps:
            d.get_all_headers(hdrs)
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

    def run_scons_cmd(self, cmd, *args, **kw):
        # Capture the log.
        old_log = self.ctx.sconf.logstream
        self.ctx.sconf.logstream = open('sconfig.log', 'w')

        # Execute the command.
        res = cmd(*args, **kw)

        # Make sure the file is closed.
        try:
            self.ctx.sconf.logstream.close()
        finally:
            pass

        # Replace the old log.
        self.ctx.sconf.logstream = old_log

        # Return results.
        log_file = open('sconfig.log', 'r')
        log = log_file.read()
        log_file.close()
        os.remove('sconfig.log')
        old_log.write(log)
        return [res, log]

    def __str__(self):
        str =  'Package name:  %s\n' % self.name
        str += 'Found:         %s\n' % self.result[0]
        str += 'Base path:     %s\n' % self.base_dir
        str += 'Header paths:  %s\n' % self.hdr_dirs
        str += 'Library paths: %s\n' % self.lib_dirs
        str += 'Libraries:     %s\n' % self.libs
        str += 'Have shared:   %s\n' % self.have_shared
        str += 'Frameworks:    %s\n' % self.fworks
        return str
