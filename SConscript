import os
Import('env')

env = env.Copy()
env.AppendUnique(CPPPATH=[env.get_build_path('include/glucifer')])

# Collect our inputs from the directory structure.
bases = ['Base', 'Windowing', 'RenderingEngines', 'OutputFormats',
         'InputFormats', 'DrawingObjects', 'WindowInteractions']
src_objs = []
suite_hdrs = []
suite_objs = []
for base in bases:
    env.build_files(env.glob(base + '/src/*.def'), 'include/glucifer/' + base)
    env.build_headers(env.glob(base + '/src/*.h'), 'include/glucifer/' + base)
    src_objs += env.build_sources(env.glob(base + '/src/*.c'), 'glucifer/' + base)
    src_objs += env.build_metas(env.glob(base + '/src/*.meta'), 'glucifer/' + base)
    suite_hdrs += env.glob(base + '/tests/*Suite.h')
    suite_objs += env.build_sources(env.glob(base + '/tests/*Suite.c'), 'glucifer/' + base)

env.build_headers(env.glob('libglucifer/src/*.h'), 'include/glucifer')
src_objs += env.build_sources(env.glob('libglucifer/src/*.c'), 'glucifer/libglucifer')
src_objs += env.build_sources(env.glob('libglucifer/src/*.meta'), 'glucifer/libglucifer')

# Build library.
if env['static_libraries']:
    env.Library(env.get_build_path('lib/glucifer'), src_objs)
if env['shared_libraries']:
    env.SharedLibrary(env.get_build_path('lib/glucifer'), src_objs)

# Build plugins.
if env['shared_libraries']:
    plgn_bases = ['lucPlugin']
    for base in plgn_bases:
        env.build_headers(env.glob('plugins/' + base + '/*.h'), 'include/glucifer/' + base.split('/')[-1])
        objs = env.build_sources(env.glob('plugins/' + base + '/*.c'), 'glucifer/' + base)
        name = base.split('/')[-1] + 'module'
        env.SharedLibrary(env.get_build_path('lib/' + name), objs,
                          SHLIBPREFIX='',
                          LIBPREFIXES=[env['LIBPREFIXES']] + [''],
                          LIBS=['glucifer'] + env.get('LIBS', []))

# Build unit test runner.
env['PCURUNNERINIT'] = ''
env['PCURUNNERSETUP'] = """StGermain_Init( &argc, &argv );
   StgDomain_Init( &argc, &argv );
   StgFEM_Init( &argc, &argv );"""
env['PCURUNNERTEARDOWN'] = """StgFEM_Finalise();
   StgDomain_Finalise();
   StGermain_Finalise();"""
runner_src = env.PCUSuiteRunner(env.get_build_path('glucifer/testglucifer.c'), suite_hdrs)
runner_obj = env.SharedObject(runner_src)
env.Program(env.get_build_path('bin/testglucifer'),
            runner_obj + suite_objs,
            LIBS=['glucifer', 'pcu'] + env.get('LIBS', []))

# Copy over XML files.
xml_bases = ['']
for base in xml_bases:
    dst = env.get_build_path('lib/StGermain/glucifer/' + base)
    for file in env.glob('ModelComponents/' + base + '/*.xml'):
        env.Install(dst, file)
