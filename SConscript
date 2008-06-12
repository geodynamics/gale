import os
Import('env')

env = env.Copy()
env.AppendUnique(CPPPATH=[env.get_build_path('include/PICellerator')])

# Collect our inputs from the directory structure.
bases = ['Voronoi', 'PopulationControl', 'Weights', 'MaterialPoints', 'Utils']
src_objs = []
suite_hdrs = []
suite_objs = []
for base in bases:
    env.build_files(env.glob(base + '/src/*.def'), 'include/PICellerator/' + base)
    env.build_headers(env.glob(base + '/src/*.h'), 'include/PICellerator/' + base)
    src_objs += env.build_sources(env.glob(base + '/src/*.c'), 'PICellerator/' + base)
    src_objs += env.build_metas(env.glob(base + '/src/*.meta'), 'PICellerator/' + base)
    suite_hdrs += env.glob(base + '/tests/*Suite.h')
    suite_objs += env.build_sources(env.glob(base + '/tests/*Suite.c'), 'PICellerator/' + base)

env.build_headers(env.glob('libPICellerator/src/*.h'), 'include/PICellerator')
src_objs += env.build_sources(env.glob('libPICellerator/src/*.c'), 'PICellerator/libPICellerator')
src_objs += env.build_sources(env.glob('libPICellerator/src/*.meta'), 'PICellerator/libPICellerator')

# Build library.
if env['static_libraries']:
    env.Library(env.get_build_path('lib/PICellerator'), src_objs)
if env['shared_libraries']:
    env.SharedLibrary(env.get_build_path('lib/PICellerator'), src_objs)

# Build toolbox.
if env['shared_libraries']:
    objs = env.build_sources(env.glob('libPICellerator/Toolbox/*.c'),
                             'PICellerator/libPICellerator/Toolbox')
    objs += env.build_metas(env.glob('libPICellerator/Toolbox/*.meta'),
                            'PICellerator/libPICellerator/Toolbox')
    env.SharedLibrary(env.get_target_name('lib/PICellerator_Toolboxmodule'), objs,
                      SHLIBPREFIX='',
                      LIBPREFIXES=env.make_list(env['LIBPREFIXES']) + [''],
                      LIBS=['PICellerator'] + env.get('LIBS', []))

# Build plugins.
if env['shared_libraries']:
    plgn_bases = ['CalculateParticleDisplacement', 'MaterialCentroid']
    for base in plgn_bases:
        env.build_headers(env.glob('plugins/' + base + '/*.h'), 'include/PICellerator/' + base.split('/')[-1])
        objs = env.build_sources(env.glob('plugins/' + base + '/*.c'), 'PICellerator/' + base)
        name = 'PICellerator_' + base.split('/')[-1] + 'module'
        env.SharedLibrary(env.get_build_path('lib/' + name), objs,
                          SHLIBPREFIX='',
                          LIBPREFIXES=env.make_list(env['LIBPREFIXES']) + [''],
                          LIBS=['PICellerator'] + env.get('LIBS', []))

# Build unit test runner.
env['PCURUNNERINIT'] = ''
env['PCURUNNERSETUP'] = """StGermain_Init( &argc, &argv );
   StgDomain_Init( &argc, &argv );
   StgFEM_Init( &argc, &argv );
   PICellerator_Init( &argc, &argv );"""
env['PCURUNNERTEARDOWN'] = """PICellerator_Finalise();
   StgFEM_Finalise();
   StgDomain_Finalise();
   StGermain_Finalise();"""
runner_src = env.PCUSuiteRunner(env.get_build_path('PICellerator/testPICellerator.c'), suite_hdrs)
runner_obj = env.SharedObject(runner_src)
env.Program(env.get_build_path('bin/testPICellerator'),
            runner_obj + suite_objs,
            LIBS=['PICellerator', 'pcu'] + env.get('LIBS', []))

# Copy over XML files.
xml_bases = ['']
for base in xml_bases:
    dst = env.get_build_path('lib/StGermain/PICellerator/' + base)
    for file in env.glob('Apps/src/' + base + '/*.xml'):
        env.Install(dst, file)
