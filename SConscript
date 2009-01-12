import os
Import('env')

#
# Prepare the construction environment by copying the one we
# were given.
env = env.Copy()
env.project_name = 'Underworld'
env.AppendUnique(CPPPATH=[env.get_build_path('include/' + env.project_name)])
env.src_objs = []
env.suite_hdrs = []
env.suite_objs = []

#
# Build standard stg directories.
env.build_directory('Rheology')
env.build_directory('Utils')

#
# Need to handle libUnderworld differently.
env.build_headers(env.glob('libUnderworld/src/*.h'), 'include/Underworld')
env.src_objs += env.build_sources(env.glob('libUnderworld/src/*.c'), 'Underworld/libUnderworld')
env.src_objs += env.build_metas(env.glob('libUnderworld/src/*.meta'), 'Underworld/libUnderworld')

#
# Build shared library.
if env['shared_libraries']:
    env.SharedLibrary(env.get_build_path('lib/Underworld'), env.src_objs)

#
# Build toolbox.
env.build_toolbox('libUnderworld/Toolbox')

#
# Build plugins. Note that this must happen after the shared library
# has been built.
env.build_plugin('plugins/EulerDeform')
env.build_plugin('plugins/ExtractPetscObjects')
env.build_plugin('plugins/IncompressibleExtensionBC')
env.build_plugin('plugins/Output/BuoyancyIntegrals')
env.build_plugin('plugins/DensityChangeAtDepth')
env.build_plugin('plugins/MaterialThermalDiffusivity')
env.build_plugin('plugins/MeshAdvectionCorrection')
env.build_plugin('plugins/VariableConditions/ShapeFemIC')
env.build_plugin('plugins/Output/Vrms')
env.build_plugin('plugins/Output/Plateness')
env.build_plugin('plugins/Output/Mobility')
env.build_plugin('plugins/Output/Nusselt')
env.build_plugin('plugins/Output/VTKOutput')
env.build_plugin('plugins/Output/MaxTemperature')
env.build_plugin('plugins/Output/MaxVelocity')
env.build_plugin('plugins/Output/BoundaryLayers')
env.build_plugin('plugins/Output/AverageTemperature')
env.build_plugin('plugins/ScalingChecks/Ra_Scaling')
env.build_plugin('SysTest/AnalyticPlugins/VelicIC')
env.build_plugin('SysTest/AnalyticPlugins/Velic_solA')
env.build_plugin('SysTest/AnalyticPlugins/Velic_solB')
env.build_plugin('SysTest/AnalyticPlugins/Velic_solCx')
env.build_plugin('SysTest/AnalyticPlugins/Velic_solHA')
env.build_plugin('SysTest/AnalyticPlugins/Velic_solKz')
env.build_plugin('SysTest/AnalyticPlugins/Velic_solS')

#
# Build static library.
if env['static_libraries']:
    env.Library(env.get_build_path('lib/Underworld'), env.src_objs)

#
# Build unit test runner.
if not env.get('dir_target', ''):
    env['PCURUNNERINIT'] = ''
    env['PCURUNNERSETUP'] = """StGermain_Init( &argc, &argv );
   StgDomain_Init( &argc, &argv );
   StgFEM_Init( &argc, &argv );
   PICellerator_Init( &argc, &argv );
   Underworld_Init( &argc, &argv );"""
    env['PCURUNNERTEARDOWN'] = """Underworld_Finalise();
   PICellerator_Finalise();
   StgFEM_Finalise();
   StgDomain_Finalise();
   StGermain_Finalise();"""
    runner_src = env.PCUSuiteRunner(env.get_build_path('Underworld/testUnderworld.c'), env.suite_hdrs)
    runner_obj = env.SharedObject(runner_src)
    env.Program(env.get_build_path('bin/testUnderworld'),
                runner_obj + env.suite_objs,
                LIBS=['Underworld', 'pcu'] + env.get('LIBS', []))

#
# Copy over XML files.
xml_bases = ['', 'BaseApps', 'VariableConditions', 'Viewports']
for base in xml_bases:
    dst = env.get_build_path('lib/StGermain/Underworld/' + base)
    for file in env.glob('InputFiles/src/' + base + '/*.xml'):
        if env.check_dir_target(file):
            env.Install(dst, file)

#
# Return any module code we need to build into a static binary.
module = (env.get('STGMODULEPROTO', ''), env.get('STGMODULECODE', ''))
Return('module')
