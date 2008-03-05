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

# Configuring or building? Or helping?
if 'config' in COMMAND_LINE_TARGETS or 'help' in COMMAND_LINE_TARGETS:
    import packages
    opts = Options('config.cache') # Create our options database.

    # Setup all the packages we want configured.
    env.Package(packages.stgUnderworld, opts)
    env.setup_packages(opts)

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
    if env['staticLibraries']: # Static libs or shared?
        env.library_builder = env.StaticLibrary
    else:
        env.library_builder = env.SharedLibrary
    env.Default(env['buildPath']) # Needed for different build paths.
    if platform.system() == 'Darwin': # Need to modify building shared libraries when on Mac OS X.
        env.AppendUnique(SHLINKFLAGS=['-flat_namespace',
                                      '-single_module',
                                      '-undefined', 'suppress'])
        import SCons.Util # And fix RPATHs.
        env['LINKFLAGS'] = SCons.Util.CLVar('')
        env['RPATHPREFIX'] = ''
        env['RPATHSUFFIX'] = ''
        env['_RPATH'] = ''

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
