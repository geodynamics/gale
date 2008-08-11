import os
Import('env')

#
# Prepare the construction environment by copying the one we
# were given.
env = env.Copy()
env.project_name = 'StGermain'
env.AppendUnique(CPPPATH=[env.get_build_path('include/' + env.project_name)])
SConscript('pcu/script/scons.py', exports='env')
SConscript('script/scons.py', exports='env')
env.src_objs = []
env.suite_hdrs = []
env.suite_objs = []

#
# These couple of files need to be explicitly copied immediately.
env.build_files(['Base/IO/src/mpirecord/none/mpimessaging.h'],
                'include/StGermain/Base/IO/mpirecord')
env.build_files(['Base/Foundation/src/ClassSetup.h', 'Base/Foundation/src/ClassEmpty.h'],
                'include/StGermain/Base/Foundation')

#
# Build standard stg directories.
env.build_directory('Base/Foundation')
env.build_directory('Base/IO')
env.build_directory('Base/Container')
env.build_directory('Base/Automation')
env.build_directory('Base/Extensibility')
env.build_directory('Base/Context')
env.build_directory('Base')
env.build_directory('Utils')

#
# Need to handle libStGermain differently.
env.build_headers(env.glob('libStGermain/src/*.h'), 'include/StGermain')
env.src_objs += env.build_sources(env.glob('libStGermain/src/*.c'), 'StGermain/libStGermain')
env.src_objs += env.build_metas(env.glob('libStGermain/src/*.meta'), 'StGermain/libStGermain')

#
# Build libraries.
if env['static_libraries']:
    env.Library(env.get_build_path('lib/StGermain'), env.src_objs)
if env['shared_libraries']:
    env.SharedLibrary(env.get_build_path('lib/StGermain'), env.src_objs)

#
# FlattenXML program.
if env['shared_libraries']:
    src = 'Base/FlattenXML/src/main.c'
    if env.check_dir_target(src):
        obj = env.SharedObject(env.get_build_path('StGermain/' + src[:-2]), src)
        env.Program(env.get_build_path('bin/FlattenXML'), obj,
                    LIBS=['StGermain'] + env.get('LIBS', []))

#
# StGermain program. We can only build this guy here if we're using
# shared libraries.
if env['shared_libraries']:
    src = 'src/main.c'
    if env.check_dir_target(src):
        obj = env.SharedObject(env.get_build_path('StGermain/' + src[:-2]), src)
        env.Program(env.get_build_path('bin/StGermain'), obj,
                    LIBS=['StGermain'] + env.get('LIBS', []))

#
# Build unit testing framework.
SConscript('pcu/SConscript', exports='env')

#
# Build unit test runner.
if not env.get('dir_target', ''):
    env['PCURUNNERINIT'] = ''
    env['PCURUNNERSETUP'] = 'StGermain_Init( &argc, &argv );'
    env['PCURUNNERTEARDOWN'] = 'StGermain_Finalise();'
    runner_src = env.PCUSuiteRunner(env.get_build_path('StGermain/testStGermain.c'),
                                    env.suite_hdrs)
    runner_obj = env.SharedObject(runner_src)
    env.Program(env.get_build_path('bin/testStGermain'),
                runner_obj + env.suite_objs,
                LIBS=['StGermain', 'pcu'] + env.get('LIBS', []))

#
# Copy scripts to correct destinations.
if env.check_dir_target('script/scons.py'):
    env.Install(env.get_build_path('script/StGermain'), 'script/scons.py')

#
# Copy XML validation file to correct destination.
if env.check_dir_target('Base/IO/src/StGermain.xsd'):
    env.Install(env.get_build_path('lib'), 'Base/IO/src/StGermain.xsd')

#
# Return any module code we need to build into a static binary.
module = (env.get('STGMODULEPROTO', ''), env.get('STGMODULECODE', ''))
Return('module')
