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
env['CPPPATH'] = env['CPPPATH'] if 'CPPPATH' in env._dict else []
env['LIBPATH'] = env['LIBPATH'] if 'LIBPATH' in env._dict else []
env['CPPDEFINES'] = env['CPPDEFINES'] if 'CPPDEFINES' in env._dict else []
env['LIBS'] = env['LIBS'] if 'LIBS' in env._dict else []
env['RPATH'] = env['RPATH'] if 'RPATH' in env._dict else []

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
env.build_directory('Shape')
env.build_directory('Mesh')
env.build_directory('Utils')
env.build_directory('Swarm')

env.build_headers(glob.glob('libStgDomain/src/*.h'), 'StgDomain')
env.build_objects(glob.glob('libStgDomain/src/*.c'), 'libStgDomain')
env.build_objects(glob.glob('libStgDomain/Toolbox/*.c'), 'Toolbox')
env.build_metadata(glob.glob('libStgDomain/src/*.meta'), 'libStgDomain')
env.build_metadata(glob.glob('libStgDomain/Toolbox/*.meta'), 'Toolbox')

env.build_library(env.obj_nodes['.'], 'StgDomain')

env.build_library(env.obj_nodes['Toolbox']['.'],
                  'StgDomain_Toolboxmodule',
                  True)

env.build_tests(glob.glob('libStgDomain/tests/test*.c'),
                'StgDomain', libs='StgDomain')
