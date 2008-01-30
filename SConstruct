import os, glob

# Check versions of some things.
EnsurePythonVersion(2, 5)
EnsureSConsVersion(0, 97)

#
# Setup our option database.
#

opts = Options('config.cache')
opts.AddOptions(
    BoolOption('debug', 'Enable debugging version', 1),
    BoolOption('useMpiRecord', 'Don''t know what this does...', 0),
    PathOption('prefix', 'Installation path',
               '/usr/local', PathOption.PathIsDirCreate),
    PathOption('mpichDir', 'MPI installation path',
               None, PathOption.PathIsDir),
    PathOption('mpichIncDir', 'MPICH header installation path',
               None, PathOption.PathIsDir),
    PathOption('mpichLibDir', 'MPICH library installation path',
               None, PathOption.PathIsDir),
    PathOption('libxml2Dir', 'libXML2 installation path',
               None, PathOption.PathIsDir),
    PathOption('libxml2IncDir', 'libXML2 header installation path',
               None, PathOption.PathIsDir),
    PackageOption('csoap', 'Enable use of the package CSoap', 'no'))

#
# Create our substitution environment.
#

env = Environment(CC='cc', ENV=os.environ, options=opts)
env['CPPDEFINES'] = env['CPPDEFINES'] if 'CPPDEFINES' in env._dict else []
env['LIBS'] = env['LIBS'] if 'LIBS' in env._dict else []

# Add any variables that get used throughout the whole build.
env.proj = 'StGermain'
if env['debug']:
    env.Append(CCFLAGS='-g')
env.Append(CPPPATH=['#build/include'])
env.Append(CPPPATH=['#build/include/' + env.proj])
env.Append(LIBPATH=['#build/lib'])

# Add any helper functions we may need.
SConscript('StgSCons', exports='env')

# Setup any additional builds.
env.Alias('install', env['prefix'])

#
# Configuration section.
#

SConscript('SConfigure', exports='env opts')

#
# Target specification section.
#

env.BuildDir('build', '.', duplicate=0)

# TODO: Build sources if csoap present.
env.build_headers(glob.glob('Regresstor/libRegresstor/src/*.h'),
                  'StGermain/Regresstor/libRegresstor')

# Some special cases we need to copy first.
env.build_headers(['Base/Foundation/src/ClassSetup.h',
                   'Base/Foundation/src/ClassEmpty.h'],
                  'StGermain/Base/Foundation', force_copy=True)

# We have this funny 'mpirecord' business to handle now.
if env['useMpiRecord']:
    # TODO
    pass
else:
    env.build_headers('Base/IO/src/mpirecord/none/mpimessaging.h',
                      'StGermain/Base/IO/mpirecord')

env.build_directory('Base/Foundation')
env.build_directory('Base/IO',
                    test_libs='StGermainBaseFoundation')
env.build_directory('Base/Container',
                    test_libs=['StGermainBaseFoundation',
                               'StGermainBaseIO'])
env.build_directory('Base/Automation',
                    test_libs=['StGermainBaseFoundation',
                               'StGermainBaseIO',
                               'StGermainBaseContainer'])
env.build_directory('Base/Extensibility',
                    test_libs=['StGermainBaseFoundation',
                               'StGermainBaseIO',
                               'StGermainBaseContainer',
                               'StGermainBaseAutomation'])
env.build_directory('Base/Context',
                    test_libs=['StGermainBaseFoundation',
                               'StGermainBaseIO',
                               'StGermainBaseContainer',
                               'StGermainBaseAutomation',
                               'StGermainBaseExtensibility'])
env.build_directory('Base',
                    extra_objs=env.obj_nodes['StGermainBaseFoundation'] + \
                        env.obj_nodes['StGermainBaseIO'] + \
                        env.obj_nodes['StGermainBaseContainer'] + \
                        env.obj_nodes['StGermainBaseAutomation'] + \
                        env.obj_nodes['StGermainBaseExtensibility'] + \
                        env.obj_nodes['StGermainBaseContext'])
env.build_directory('Utils', test_libs='StGermainBase')

env.build_headers(glob.glob('libStGermain/src/*.h'), 'StGermain')
env.build_objects(glob.glob('libStGermain/src/*.c'), 'libStGermain')
env.build_metadata(glob.glob('libStGermain/src/*.meta'), 'libStGermain')

allNodes = env.obj_nodes['StGermainBaseFoundation'] + \
    env.obj_nodes['StGermainBaseIO'] + \
    env.obj_nodes['StGermainBaseContainer'] + \
    env.obj_nodes['StGermainBaseAutomation'] + \
    env.obj_nodes['StGermainBaseExtensibility'] + \
    env.obj_nodes['StGermainBaseContext'] + \
    env.obj_nodes['StGermainUtils'] + \
    env.obj_nodes['libStGermain']
env.build_library(allNodes, 'StGermain')

env.build_tests(glob.glob('libStGermain/tests/test*.c'), 'StGermain',
                libs='StGermain')

env.Program('build/bin/FlattenXML', 'build/Base/FlattenXML/src/main.c',
            LIBS=env['LIBS'] + ['StGermainBase'])
env.Program('build/bin/StGermain', 'build/src/main.c',
            LIBS=env['LIBS'] + ['StGermain'])
