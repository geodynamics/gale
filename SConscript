import os
Import('env')

#
# Prepare the construction environment by copying the one we
# were given.
env = env.Copy()
env.project_name = 'StgFEM'
env.AppendUnique(CPPPATH=[env.get_build_path('include/' + env.project_name)])
env.src_objs = []
env.suite_hdrs = []
env.suite_objs = []

#
# Build standard stg directories.
env.build_directory('Discretisation')
# env.build_directory('SLE/LinearAlgebra')
env.build_directory('SLE/SystemSetup')
env.build_directory('SLE/ProvidedSystems/AdvectionDiffusion')
env.build_directory('SLE/ProvidedSystems/Energy')
env.build_directory('SLE/ProvidedSystems/StokesFlow')
env.build_directory('SLE/ProvidedSystems')
env.build_directory('SLE')
env.build_directory('Assembly')

#
# Need to handle libStgFEM differently.
env.build_headers(env.glob('libStgFEM/src/*.h'), 'include/StgFEM')
env.src_objs += env.build_sources(env.glob('libStgFEM/src/*.c'), 'StgFEM/libStgFEM')
env.src_objs += env.build_metas(env.glob('libStgFEM/src/*.meta'), 'StgFEM/libStgFEM')

#
# Build shared library.
if env['shared_libraries']:
    env.SharedLibrary(env.get_build_path('lib/StgFEM'), env.src_objs)

#
# Build toolbox.
env.build_toolbox('libStgFEM/Toolbox')

#
# Build plugins. Note that this must happen after the shared library
# has been built.
env.build_plugin('plugins/CompareFeVariableAgainstReferenceSolution')
env.build_plugin('plugins/Document')
env.build_plugin('plugins/FeVariableImportExporters/FeVariable_ImportExport_ABAQUS')
env.build_plugin('plugins/FeVariableImportExporters/FeVariable_ImportExport_SpecRidge2D')
env.build_plugin('plugins/FileAnalyticSolution')
env.build_plugin('plugins/Output/CPUTime')
env.build_plugin('plugins/Output/FrequentOutput')
env.build_plugin('plugins/Output/FeVariableList')
env.build_plugin('plugins/Output/PeakMemory')
env.build_plugin('plugins/Output/PrintFeVariableDiscreteValues')
env.build_plugin('plugins/Output/PrintFeVariableDiscreteValues_2dBox')
env.build_plugin('plugins/StandardConditionFunctions')
env.build_plugin('Apps/StokesMomentumUzawa/tests/LinearVelocityAnalytic')
env.build_plugin('Apps/StokesMomentumUzawa/tests/LidDrivenIsoviscousAnalytic')
env.build_plugin('Apps/StokesMomentumUzawa/tests/SimpleShearAnalytic')
env.build_plugin('Apps/StokesMomentumUzawa/tests/LidDrivenStokesAnalytic')

#
# Build static library.
if env['static_libraries']:
    env.Library(env.get_build_path('lib/StgFEM'), env.src_objs)

#
# Build unit test runner.
if not env.get('dir_target', ''):
    env['PCURUNNERINIT'] = ''
    env['PCURUNNERSETUP'] = """StGermain_Init( &argc, &argv );
   StgDomain_Init( &argc, &argv );
   StgFEM_Init( &argc, &argv );"""
    env['PCURUNNERTEARDOWN'] = """StgFEM_Finalise();
   StgDomain_Finalise();
   StGermain_Finalise();"""
    runner_src = env.PCUSuiteRunner(env.get_build_path('StgFEM/testStgFEM.c'), env.suite_hdrs)
    runner_obj = env.SharedObject(runner_src)
    env.Program(env.get_build_path('bin/testStgFEM'),
                runner_obj + env.suite_objs,
                LIBS=['StgFEM', 'pcu'] + env.get('LIBS', []))

#
# Copy over XML files.
xml_bases = ['']
for base in xml_bases:
    dst = env.get_build_path('lib/StGermain/StgFEM/' + base)
    for file in env.glob('Apps/StgFEM_Components/' + base + '/*.xml'):
        if env.check_dir_target(file):
            env.Install(dst, file)

#
# Return any module code we need to build into a static binary.
module = (env.get('STGMODULEPROTO', ''), env.get('STGMODULECODE', ''))
Return('module')
