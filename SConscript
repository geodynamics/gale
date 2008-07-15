import os
Import('env')

#
# Prepare the construction environment by copying the one we
# were given.
env = env.Copy()
env.project_name = 'StgDomain'
env.AppendUnique(CPPPATH=[env.get_build_path('include/' + env.project_name)])
env.src_objs = []
env.suite_hdrs = []
env.suite_objs = []

#
# Build standard stg directories.
env.build_directory('Geometry')
env.build_directory('Shape')
env.build_directory('Mesh')
env.build_directory('Utils')
env.build_directory('Swarm')

#
# Need to handle libStGermain differently.
env.build_headers(env.glob('libStgDomain/src/*.h'), 'include/StgDomain')
env.src_objs += env.build_sources(env.glob('libStgDomain/src/*.c'), 'StgDomain/libStgDomain')
env.src_objs += env.build_metas(env.glob('libStgDomain/src/*.meta'), 'StgDomain/libStgDomain')

# Build shared library.
if env['shared_libraries']:
    env.SharedLibrary(env.get_build_path('lib/StgDomain'), env.src_objs)

# Build toolbox.
objs = env.build_sources(env.glob('libStgDomain/Toolbox/*.c'), 'StgDomain/libStgDomain/Toolbox')
objs += env.build_metas(env.glob('libStgDomain/Toolbox/*.meta'), 'StgDomain/libStgDomain/Toolbox')
if env['shared_libraries']:
    env.SharedLibrary(env.get_target_name('lib/StgDomain_Toolboxmodule'), objs,
                      SHLIBPREFIX='',
                      LIBPREFIXES=env.make_list(env['LIBPREFIXES']) + [''],
                      LIBS=['StgDomain'] + env.get('LIBS', []))
if env['static_libraries']:
    env.src_objs += objs

#
# Build static library.
if env['static_libraries']:
    env.Library(env.get_build_path('lib/StgDomain'), env.src_objs)

# Build unit test runner.
env['PCURUNNERINIT'] = ''
env['PCURUNNERSETUP'] = """StGermain_Init( &argc, &argv );
   StgDomain_Init( &argc, &argv );"""
env['PCURUNNERTEARDOWN'] = """StgDomain_Finalise();
   StGermain_Finalise();"""
runner_src = env.PCUSuiteRunner(env.get_build_path('StgDomain/testStgDomain.c'), env.suite_hdrs)
runner_obj = env.SharedObject(runner_src)
env.Program(env.get_build_path('bin/testStgDomain'),
            runner_obj + env.suite_objs,
            LIBS=['StgDomain', 'pcu'] + env.get('LIBS', []))
