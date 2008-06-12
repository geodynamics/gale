import os
Import('env')

env = env.Copy()
env.AppendUnique(CPPPATH=[env.get_build_path('include/StgDomain')])

# Collect our inputs from the directory structure.
bases = ['Geometry', 'Shape', 'Mesh', 'Utils', 'Swarm']
src_objs = []
suite_hdrs = []
suite_objs = []
for base in bases:
    env.build_files(env.glob(base + '/src/*.def'), 'include/StgDomain/' + base)
    env.build_headers(env.glob(base + '/src/*.h'), 'include/StgDomain/' + base)
    src_objs += env.build_sources(env.glob(base + '/src/*.c'), 'StgDomain/' + base)
    src_objs += env.build_metas(env.glob(base + '/src/*.meta'), 'StgDomain/' + base)
    suite_hdrs += env.glob(base + '/tests/*Suite.h')
    suite_objs += env.build_sources(env.glob(base + '/tests/*Suite.c'), 'StgDomain/' + base)

env.build_headers(env.glob('libStgDomain/src/*.h'), 'include/StgDomain')
src_objs += env.build_sources(env.glob('libStgDomain/src/*.c'), 'StgDomain/libStgDomain')
src_objs += env.build_sources(env.glob('libStgDomain/src/*.meta'), 'StgDomain/libStgDomain')

# Build library.
if env['static_libraries']:
    env.Library(env.get_build_path('lib/StgDomain'), src_objs)
if env['shared_libraries']:
    env.SharedLibrary(env.get_build_path('lib/StgDomain'), src_objs)

# Build toolbox.
if env['shared_libraries']:
    objs = env.build_sources(env.glob('libStgDomain/Toolbox/*.c'),
                             'StgDomain/libStgDomain/Toolbox')
    objs += env.build_metas(env.glob('libStgDomain/Toolbox/*.meta'),
                            'StgDomain/libStgDomain/Toolbox')
    env.SharedLibrary(env.get_target_name('lib/StgDomain_Toolboxmodule'), objs,
                      SHLIBPREFIX='',
                      LIBPREFIXES=env.make_list(env['LIBPREFIXES']) + [''],
                      LIBS=['StgDomain'] + env.get('LIBS', []))

# Build unit test runner.
env['PCURUNNERINIT'] = ''
env['PCURUNNERSETUP'] = """StGermain_Init( &argc, &argv );
   StgDomain_Init( &argc, &argv );"""
env['PCURUNNERTEARDOWN'] = """StgDomain_Finalise();
   StGermain_Finalise();"""
runner_src = env.PCUSuiteRunner(env.get_build_path('StgDomain/testStgDomain.c'), suite_hdrs)
runner_obj = env.SharedObject(runner_src)
env.Program(env.get_build_path('bin/testStgDomain'),
            runner_obj + suite_objs,
            LIBS=['StgDomain', 'pcu'] + env.get('LIBS', []))
