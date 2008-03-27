import os, sys, platform

# Source the configuration scripts.
sys.path.insert(0, os.path.abspath(os.path.join(os.getcwd(), 'config')))
SConscript('config/SConfig/SConscript')

# Create our base construction environment.
if platform.platform().find('ia64') != -1:
    # Hack needed for APAC, damn it all!
    env = Environment(ENV=os.environ, tools=['gcc', 'gnulink'])
else:
    env = Environment(ENV=os.environ)
if 'CC' in env['ENV']:
    env['CC'] = env['ENV']['CC']

# Need to modify building shared libraries when on Mac OS X.
if platform.system() == 'Darwin':
    env.AppendUnique(SHLINKFLAGS=['-flat_namespace',
                                  '-single_module',
                                  '-undefined', 'suppress'])
    import SCons.Util # And fix RPATHs.
    env['LINKFLAGS'] = SCons.Util.CLVar('')
    env['RPATHPREFIX'] = ''
    env['RPATHSUFFIX'] = ''
    env['_RPATH'] = ''

# Configuring or building? Or helping?
if 'config' in COMMAND_LINE_TARGETS or 'help' in COMMAND_LINE_TARGETS:
    import packages
    opts = Options() # Create our options database.

    # Setup all the packages we want configured.
    env.Package(packages.stgUnderworld, opts, True)

    # Displaying help?
    if 'help' in COMMAND_LINE_TARGETS:
        env.Alias('help', '.')
        print opts.GenerateHelpText(env)

    # Or configuring?
    else:
        env.configure_packages(opts)

else:
    # Prepare our construction environment.
    env.load_config('config.cfg') # Load configuration.
    SConscript('StgSCons', exports='env') # Setup our StG specific utils.
    env.Default(env['buildPath']) # Needed for different build paths.
    env.Append(CPPDEFINES=['HAVE_SCONS'])

    # Another Mac OS X hack.
    if platform.system() == 'Darwin':
	env['_abspath'] = lambda x: File(x).abspath
        env.Append(SHLINKFLAGS=['-install_name', '${_abspath(TARGET)}'])

    # Prepare library builders.
    env.library_builders = []
    if env['static_libraries']:
        env.library_builders += [env.StaticLibrary]
    if env['shared_libraries']:
        env.library_builders += [env.SharedLibrary]

    # Specify targets.
    SConscript('StGermain/SConscript', exports='env')
    env.Prepend(LIBS='StGermain')
    SConscript('StgDomain/SConscript', exports='env')
    env.Prepend(LIBS='StgDomain')
    SConscript('StgFEM/SConscript', exports='env')
    env.Prepend(LIBS='StgFEM')
    SConscript('PICellerator/SConscript', exports='env')
    env.Prepend(LIBS='PICellerator')
    SConscript('Underworld/SConscript', exports='env')
    if env['with_glucifer']:
        SConscript('gLucifer/SConscript', exports='env')
