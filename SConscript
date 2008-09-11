import os
Import('env')

#
# Prepare the construction environment by copying the one we
# were given.
env = env.Copy()
env.project_name = 'PICellerator'
env.AppendUnique(CPPPATH=[env.get_build_path('include/' + env.project_name)])
env.src_objs = []
env.suite_hdrs = []
env.suite_objs = []

#
# Build standard stg directories.
env.build_directory('Voronoi')
env.build_directory('PopulationControl')
env.build_directory('Weights')
env.build_directory('MaterialPoints')
env.build_directory('Utils')

#
# Need to handle libPICellerator differently.
env.build_headers(env.glob('libPICellerator/src/*.h'), 'include/PICellerator')
env.src_objs += env.build_sources(env.glob('libPICellerator/src/*.c'), 'PICellerator/libPICellerator')
env.src_objs += env.build_metas(env.glob('libPICellerator/src/*.meta'), 'PICellerator/libPICellerator')

#
# Build shared library.
if env['shared_libraries']:
    env.SharedLibrary(env.get_build_path('lib/PICellerator'), env.src_objs)

#
# Build toolbox.
env.build_toolbox('libPICellerator/Toolbox')

#
# Build plugins. Note that this must happen after the shared library
# has been built.
env.build_plugin('plugins/CalculateParticleDisplacement')
env.build_plugin('plugins/Output/MaterialCentroid')

#
# Build static library.
if env['static_libraries']:
    env.Library(env.get_build_path('lib/PICellerator'), env.src_objs)

#
# Build unit test runner.
if not env.get('dir_target', ''):
    env['PCURUNNERINIT'] = ''
    env['PCURUNNERSETUP'] = """StGermain_Init( &argc, &argv );
   StgDomain_Init( &argc, &argv );
   StgFEM_Init( &argc, &argv );
   PICellerator_Init( &argc, &argv );"""
    env['PCURUNNERTEARDOWN'] = """PICellerator_Finalise();
   StgFEM_Finalise();
   StgDomain_Finalise();
   StGermain_Finalise();"""
    runner_src = env.PCUSuiteRunner(env.get_build_path('PICellerator/testPICellerator.c'), env.suite_hdrs)
    runner_obj = env.SharedObject(runner_src)
    env.Program(env.get_build_path('bin/testPICellerator'),
            runner_obj + env.suite_objs,
            LIBS=['PICellerator', 'pcu'] + env.get('LIBS', []))

#
# Copy over XML files.
xml_bases = ['']
for base in xml_bases:
    dst = env.get_build_path('lib/StGermain/PICellerator/' + base)
    for file in env.glob('Apps/src/' + base + '/*.xml'):
        if env.check_dir_target(file):
            env.Install(dst, file)

#
# Return any module code we need to build into a static binary.
module = (env.get('STGMODULEPROTO', ''), env.get('STGMODULECODE', ''))
Return('module')
