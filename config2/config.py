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

    env.AddMethod(ConfigurePackage)
    env['packages'] = {}

def exists(env):
    return True
