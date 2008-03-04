import os, sys

# Source the configuration scripts.
sys.path.insert(0, os.path.abspath(os.path.join(os.getcwd(), 'config')))
SConscript('config/SConfig/SConscript')

# Create our base construction environment.
env = Environment()

# Configuring or building? Or helping?
if 'config' in COMMAND_LINE_TARGETS or 'help' in COMMAND_LINE_TARGETS:
    import packages
    opts = Options('config.cache') # Create our options database.

    # Setup all the packages we want configured.
    env.Package(packages.StgDomain, opts)
    env.setup_packages(opts)

    if 'help' in COMMAND_LINE_TARGETS:
        env.Alias('help', '.')
        print opts.GenerateHelpText(env)

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

    SConscript('SConscript', exports='env')
