import os
Import('env')

#
# Prepare the construction environment by copying the one we
# were given.
env = env.Copy()
env.project_name = 'glucifer'
env.AppendUnique(CPPPATH=[env.get_build_path('include/' + env.project_name)])
env.src_objs = []
env.suite_hdrs = []
env.suite_objs = []

#
# Build standard stg directories.
env.build_directory('Base')
env.build_directory('Windowing')
env.build_directory('RenderingEngines')
env.build_directory('OutputFormats')
env.build_directory('InputFormats')
env.build_directory('DrawingObjects')
env.build_directory('WindowInteractions')

#
# Need to handle libglucifer differently.
env.build_headers(env.glob('libglucifer/src/*.h'), 'include/glucifer')
env.src_objs += env.build_sources(env.glob('libglucifer/src/*.c'), 'glucifer/libglucifer')
env.src_objs += env.build_metas(env.glob('libglucifer/src/*.meta'), 'glucifer/libglucifer')

#
# Build shared library.
if env['shared_libraries']:
    env.SharedLibrary(env.get_build_path('lib/glucifer'), env.src_objs)

#
# Build plugins. Note that this must happen after the shared library
# has been built.
env.build_plugin('plugins/lucPlugin')

#
# Build static library.
if env['static_libraries']:
    env.Library(env.get_build_path('lib/glucifer'), env.src_objs)

#
# Build unit test runner.
env['PCURUNNERINIT'] = ''
env['PCURUNNERSETUP'] = """StGermain_Init( &argc, &argv );
   StgDomain_Init( &argc, &argv );
   StgFEM_Init( &argc, &argv );"""
env['PCURUNNERTEARDOWN'] = """StgFEM_Finalise();
   StgDomain_Finalise();
   StGermain_Finalise();"""
runner_src = env.PCUSuiteRunner(env.get_build_path('glucifer/testglucifer.c'), env.suite_hdrs)
runner_obj = env.SharedObject(runner_src)
env.Program(env.get_build_path('bin/testglucifer'),
            runner_obj + env.suite_objs,
            LIBS=['glucifer', 'pcu'] + env.get('LIBS', []))

#
# Copy over XML files.
xml_bases = ['']
for base in xml_bases:
    dst = env.get_build_path('lib/StGermain/glucifer/' + base)
    for file in env.glob('ModelComponents/' + base + '/*.xml'):
        env.Install(dst, file)
