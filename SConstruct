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
    PathOption('petscDir', 'PETSc installation path',
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
env.proj = 'StgDomain'
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

env.build_directory('Geometry')
env.build_directory('Shape', test_libs=['StgDomainGeometry'])
env.build_directory('Mesh', test_libs=['StgDomainGeometry',
                                       'StgDomainShape'])
env.build_directory('Utils', test_libs=['StgDomainGeometry',
                                        'StgDomainShape',
                                        'StgDomainMesh'])
env.build_directory('Swarm', test_libs=['StgDomainGeometry',
                                        'StgDomainShape',
                                        'StgDomainMesh',
                                        'StgDomainUtils'])

env.build_headers(glob.glob('libStgDomain/src/*.h'), 'StgDomain')
env.build_objects(glob.glob('libStgDomain/src/*.c'), 'libStgDomain')
env.build_objects(glob.glob('libStgDomain/Toolbox/*.c'),
                  'StgDomain_Toolboxmodule')
env.build_metadata(glob.glob('libStgDomain/src/*.meta'), 'libStgDomain')
env.build_metadata(glob.glob('libStgDomain/Toolbox/*.meta'),
                   'StgDomain_Toolboxmodule')

allNodes = env.obj_nodes['StgDomainGeometry'] + \
    env.obj_nodes['StgDomainShape'] + \
    env.obj_nodes['StgDomainMesh'] + \
    env.obj_nodes['StgDomainUtils'] + \
    env.obj_nodes['StgDomainSwarm'] + \
    env.obj_nodes['libStgDomain']
env.build_library(allNodes, 'StgDomain')

env.build_library(env.obj_nodes['StgDomain_Toolboxmodule'],
                  'StgDomain_Toolboxmodule',
                  True)

env.build_tests(glob.glob('libStGermain/test/test*.c'),
                'StgDomain', 'StgDomain')
