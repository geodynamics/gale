import os, SCons.SConf, checks

class package_check:

    def __init__(self, env, name):
        self.env = env
        self.name = name
        self.result = False
        self.locations = setup_locations(self.env, name.lower())

    def __enter__(self):
        # Show an initial message.
        SCons.SConf.progress_display('Checking for package %s... '%self.name,
                                     append_newline=0)

        # Switch off intermediate display.
        SCons.SConf.progress_display.set_mode(0)

        return self

    def __exit__(self, type, value, traceback):
        SCons.SConf.progress_display.set_mode(1)
        res = self.result and 'yes' or 'no'
        SCons.SConf.progress_display(res)
        
        # Add an entry into the environment's found_packages value.
        self.env.AppendUnique(found_packages=[self.name])

def setup_conf(env):
    return env.Configure(
        custom_tests={
            'CheckLibs': checks.CheckLibs,
            'CheckLibsWithHeader': checks.CheckLibsWithHeader
            },
        )

def get_option(env, name):
    # First check command line.
    val = env.GetOption(name)
    if not val:
        # Now check environment.
        return os.environ.get(name.upper(), None)
    return val

def setup_locations(env, name):
    base = get_option(env, name + '_dir')
    if base:
        return [(base, [os.path.join(base, 'include')],
                 [os.path.join(base, 'lib')])]

    inc_dir = get_option(env, name + '_inc_dir')
    lib_dir = get_option(env, name + '_lib_dir')
    if inc_dir and lib_dir:
        return [('', [inc_dir], [lib_dir])]

    return None

def use_location(env, loc):
    env.AppendUnique(CPPPATH=loc[1])
    env.AppendUnique(LIBPATH=loc[2])

def print_failed_message(name):
    print "Failed to locate required package %s."%name
