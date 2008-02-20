import os

# Check versions of some things.
EnsurePythonVersion(2, 5)
EnsureSConsVersion(0, 97)

#
# Create our substitution environment.
#

env = Environment(ENV=os.environ)

# Include our SCons utils (hopefully they'll add stuff to do 
# these things we need).
SConscript('SConsUtils', exports='env')

# Are we configuring or building?
if ('config' in COMMAND_LINE_TARGETS or 'help' in COMMAND_LINE_TARGETS) \
        and not env.GetOption('clean'):
    SConscript('SConfigure', exports='env')
    env.Alias('config', '.')
    env.Alias('help', '.')

else:
    # Setup all our StGermain specific utilities.
    SConscript('StgSCons', exports='env')

    # Load the configuration.
    env.load_cfg('config.state')

    # Add any variables that get used throughout the whole build.
    if env['debug']:
        env.Append(CCFLAGS='-g')
    env.Prepend(CPPPATH=['#build/include'])
    env.Prepend(LIBPATH=['#build/lib'])

    # Setup any additional build targets.
    env.Alias('install', env['prefix'])

    #
    # Call target SConscripts.
    #

    env.BuildDir('build', '.', duplicate=0)

    env.clear_all()
    SConscript('build/StGermain/SConscript', exports='env')
    env.Prepend(LIBS=['StGermain'])

    env.clear_all()
    SConscript('build/StgDomain/SConscript', exports='env')
    env.Prepend(LIBS=['StgDomain'])

    env.clear_all()
    SConscript('build/StgFEM/SConscript', exports='env')
    env.Prepend(LIBS=['StgFEM'])

    env.clear_all()
    SConscript('build/PICellerator/SConscript', exports='env')
    env.Prepend(LIBS=['PICellerator'])

    env.clear_all()
    SConscript('build/Underworld/SConscript', exports='env')
    env.Prepend(LIBS=['Underworld'])

#    env.clear_all()
#    SConscript('build/Experimental/SConscript', exports='env')
#    env.Prepend(LIBS=['ExperimentalUnderworld'])

    env.clear_all()
    SConscript('build/gLucifer/SConscript', exports='env')
