Import('env')

env = env.Copy()
env.project_name = 'StgDomain' # Need a project name.
env.clear_all() # ... so that our structures are ready.

env.build_directory('Geometry')
env.build_directory('Shape')
env.build_directory('Mesh')
env.build_directory('Utils')
env.build_directory('Swarm')

env.build_headers(env.glob('libStgDomain/src/*.h'), 'StgDomain')
env.build_objects(env.glob('libStgDomain/src/*.c'), 'libStgDomain')
env.build_objects(env.glob('libStgDomain/Toolbox/*.c'), 'Toolbox')
env.build_metadata(env.glob('libStgDomain/src/*.meta'), 'libStgDomain')
env.build_metadata(env.glob('libStgDomain/Toolbox/*.meta'), 'Toolbox')

env.build_library(env.get_hnodes(env.SharedObject), 'StgDomain')

env.build_library(env.get_hnodes(env.SharedObject, 'Toolbox'),
                  'StgDomain_Toolboxmodule', ['StgDomain'],
                  True)

env.build_tests(env.glob('libStgDomain/tests/test*.c'),
                'StgDomain', libs='StgDomain')
