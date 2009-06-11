from SCons.Script import *

def UsePackage(env, mod, **kw):
    name = mod.__module__[mod.__module__.rfind('.') + 1:]
    pkg = env['packages'].get(mod, None)
    if not pkg:
        pkg = mod(name, env, **kw)
        env['packages'][mod] = pkg
    return pkg

def ConfigurePackage(env, mod, **kw):
    # Don't configure if we're cleaning or helping.
    pkg = env['packages'][mod]
    if not (GetOption('clean') or GetOption('help')):
        env._dict.update(pkg(**kw)._dict)
    return pkg

def SaveConfig(env, filename='config.cfg'):
    if not (GetOption('help') or GetOption('clean')):

        out = open(filename, 'w')
        opts = [o[1] for o in env['cfg_options']] + [
            'CPPPATH', 'LIBPATH', 'RPATH', 'LIBS', 'CPPDEFINES'
            ]
        for o in opts:
            v = env.get(o, None)
            if v is not None:
                out.write('%s = %s\n'%(o, repr(env[o])))
        out.close()

def generate(env, options=[]):
    import platform

    # Add the options to SCons.
    env['cfg_options'] = options
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

    env.AddMethod(UsePackage)
    env.AddMethod(ConfigurePackage)
    env.AddMethod(SaveConfig)
    env['packages'] = {}

def exists(env):
    return True
