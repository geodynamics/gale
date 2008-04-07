import os
Import('env')

env = env.Copy()
env.AppendUnique(CPPPATH=[env.get_build_path('include/StGermain')])

# Forced copy.
env.copy_file(env.get_build_path('include/StGermain/Base/IO/mpirecord/mpimessaging.h'),
              'Base/IO/src/mpirecord/none/mpimessaging.h')
env.build_files(['Base/Foundation/src/ClassSetup.h', 'Base/Foundation/src/ClassEmpty.h'],
                'include/StGermain/Base/Foundation')

bases = ['Base/Foundation',
         'Base/IO',
         'Base/Container',
         'Base/Automation',
         'Base/Extensibility',
         'Base/Context',
         'Base',
         'Utils']
src_objs = []
suite_hdrs = []
suite_objs = []
for base in bases:
    env.build_files(env.glob(base + '/src/*.def'), 'include/StGermain/' + base)
    env.build_headers(env.glob(base + '/src/*.h'), 'include/StGermain/' + base)
    src_objs += env.build_sources(env.glob(base + '/src/*.c'), 'StGermain/' + base)
    src_objs += env.build_metas(env.glob(base + '/src/*.meta'), 'StGermain/' + base)
    suite_hdrs += env.glob(base + '/tests/*Suite.h')
    suite_objs += env.build_sources(env.glob(base + '/tests/*Suite.c'), 'StGermain/' + base)

env.build_headers(env.glob('libStGermain/src/*.h'), 'include/StGermain')
src_objs += env.build_sources(env.glob('libStGermain/src/*.c'), 'StGermain/libStGermain')
src_objs += env.build_sources(env.glob('libStGermain/src/*.meta'), 'StGermain/libStGermain')

# Build library.
if env['static_libraries']:
    env.Library(env.get_build_path('lib/StGermain'), src_objs)
if env['shared_libraries']:
    env.SharedLibrary(env.get_build_path('lib/StGermain'), src_objs)

# FlattenXML program.
src = 'Base/FlattenXML/src/main.c'
obj = env.SharedObject(env.get_build_path('StGermain/' + src[:-2]), src)
env.Program(env.get_build_path('bin/FlattenXML'), obj,
            LIBS=['StGermain'] + env.get('LIBS', []))

# StGermain program.
src = 'src/main.c'
obj = env.SharedObject(env.get_build_path('StGermain/' + src[:-2]), src)
env.Program(env.get_build_path('bin/StGermain'), obj,
            LIBS=['StGermain'] + env.get('LIBS', []))

# Build pcu.
SConscript('pcu/SConscript', exports='env')

# Build unit test runner.
env['PCURUNNERINIT'] = ''
env['PCURUNNERSETUP'] = 'StGermain_Init( &argc, &argv );'
env['PCURUNNERTEARDOWN'] = 'StGermain_Finalise();'
runner_src = env.PCUSuiteRunner(env.get_build_path('StGermain/testStGermain.c'), suite_hdrs)
runner_obj = env.SharedObject(runner_src)
env.Program(env.get_build_path('bin/testStGermain'),
            runner_obj + suite_objs,
            LIBS=['StGermain', 'pcu'] + env.get('LIBS', []))

# Copy scripts to correct destinations.
env.Install(env.get_build_path('script/StGermain'), 'script/scons.py')
