from SCons.Script import *

def ConfigurePackage(env, mod, **kw):
    name = mod.__module__[mod.__module__.rfind('.') + 1:]
    pkg = env['packages'].get(mod, None)
    if not pkg:
        pkg = mod(name, env, **kw)
        env['packages'][mod] = pkg
        # Don't configure if we're cleaning or helping.
        if not (GetOption('clean') or GetOption('help')):
            env._dict.update(pkg()._dict)
    return pkg

def generate(env, options=[]):
    import platform

    # Add the options to SCons.
    for o in options:
        if len(o) > 4:
            type = o[4]
        else:
            type = 'string'
        AddOption(o[0], dest=o[1], nargs=1, type=type,
                  action='store', help=o[2], default=o[3])

    # Store option values on the environment.
    for o in options:
        if GetOption(o[1]) is not None:
            env[o[1]] = GetOption(o[1])
        elif o[1].upper() in os.environ:
            env[o[1]] = os.environ[o[1].upper()]

    # Need to handle Darwin shared libraries.
    if platform.system() == 'Darwin':
        env['_RPATH'] = ''
        env.Append(SHLINKFLAGS=['-dynamiclib', '-flat_namespace',
                                '-single_module', '-undefined', 'suppress',
                                '-install_name', '${_abspath(TARGET)}'])
        env['_abspath']=lambda x: File(x).abspath

    env.AddMethod(ConfigurePackage)
    env['packages'] = {}

def exists(env):
    return True
