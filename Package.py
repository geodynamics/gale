import os, SCons.SConf, checks
import utils, platform

class Package:

    def __init__(self, name, env, **kw):
        self.name = name
        self.env = env
        self.required = kw.get('required', True)
        self.result = None
        self.setup_dependencies()
        self.setup_options()

    def setup_dependencies(self):
        # self.add_dependency(...)
        pass

    def setup_options(self):
        from SCons.Script.Main import AddOption
        n = self.name
        ln = n.lower()
        AddOption('--%s-dir'%ln, dest='%s_dir'%ln, nargs=1, type='string',
                  action='store', help='Base directory for %s installation.'%n)
        AddOption('--%s-inc-dir'%ln, dest='%s_inc_dir'%ln, nargs=1, type='string',
                  action='store', help='Directory for %s headers.'%n)
        AddOption('--%s-lib-dir'%ln, dest='%s_lib_dir'%ln, nargs=1, type='string',
                  action='store', help='Directory for %s libraries.'%n)
        AddOption('--%s-libs'%ln, dest='%s_libs'%ln, nargs=1, type='string',
                  action='store', help='%s libraries.'%n)

    def gen_locations(self):
	yield

    def gen_envs(self, loc):
        env = self.env.Clone()
        env.AppendUnique(CPPPATH=loc[1])
        env.AppendUnique(LIBPATH=loc[2])
        env.AppendUnique(RPATH=loc[2])
        yield env

    def check(self, conf, env):
        return conf.CheckLibWithHeader(None, env.get('pkg_headers', []), 'c',
                                       autoadd=0)

    def __call__(self, **kw):

        # Don't run twice.
        if self.result is not None:
            return self.env

        # Process dependencies.
        if hasattr(self, 'deps'):
            for d in self.deps:
                d()

        # Check for overrides.
        if 'required' in kw:
            self.required = kw['required']

        # Handle dependencies.
        self.setup_dependencies()

        # Show an initial message.
        SCons.SConf.progress_display('Checking for package %s... '%self.name,
                                     append_newline=0)

        # Switch off intermediate display.
        SCons.SConf.progress_display.set_mode(0)

        # Initialise the result to False.
        self.result = False

        for loc in self._gen_locations():

            # Stash on the package.
            self.location = loc

            for env in self.gen_envs(loc):

                # Setup the  environment.
                conf = env.Configure(
                    custom_tests={
                        'CheckLibs': checks.CheckLibs,
                        'CheckLibsWithHeader': checks.CheckLibsWithHeader
                        }
                    )

                # Run the checks and remove temporary header variable.
                self.result = self.check(conf, env)
                if 'pkg_headers' in env:
                    del env['pkg_headers']

                # If successful, export a definition.
                if self.result:
                    env.AppendUnique(CPPDEFINES=['HAVE_' + self.name.upper()])

                # Finish the check.
                conf.Finish()
                if self.result:
                    break

            if self.result:
                break

        # Reset the display.
        SCons.SConf.progress_display.set_mode(1)
        res = self.result and 'yes' or 'no'
        SCons.SConf.progress_display(res)

        # If this package was required and it failed print out a message.
        if not self.result and self.required:
            print
            print 'Failed to locate required package %s.'%self.name
            print 'Details in \'config/config2/config.log\''
            print
            self.env.Exit()

        return self.result and env or self.env

    def add_dependency(self, mod, **kw):
        if not hasattr(self, 'deps'):
            self.deps = []
        self.deps.append(self.env.UsePackage(mod, **kw))
        return self.deps[-1]

    def get_option(self, name):
        # First check command line.
        val = self.env.GetOption(name)
        if not val:
            # Now check environment.
            return os.environ.get(name.upper(), None)
        return val

    def find_libraries(self, dirs, libs):
        # If no directories were given, return okay.
        if not dirs:
            return True

        prefixes = ['lib']
        suffixes = ['.a']
        if platform.system() == 'Linux':
            suffixes += ['.so']
        if platform.system() == 'Darwin':
            suffixes += ['.dylib']
        if platform.system() == 'Windows':
            prefixes += ['']
            suffixes += ['.dll']

        libs = utils.conv.to_list(libs)
        for l in libs:
            if not [p for p in utils.path.find(l, dirs, prefixes, suffixes)]:
                return False
        return True

    def _gen_locations(self):
        ln = self.name.lower()
        base = self.get_option(ln + '_dir')
        if base:
            yield (base, [os.path.join(base, 'include')],
                   [os.path.join(base, 'lib')])

        else:
            inc_dir = self.get_option(ln + '_inc_dir')
            lib_dir = self.get_option(ln + '_lib_dir')
            if inc_dir and lib_dir:
                yield ('', [inc_dir], [lib_dir])

            else:
                yield ('', [], [])
                for loc in self.gen_locations():
                    yield loc
