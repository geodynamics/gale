import os
Import('env')

# Need to make a copy because SCons uses the environment
# at it's final state, so StGermain ends up depending on
# StgDomain, etc.
env = env.Clone()

# Inside each project we will be accessing headers without the
# project name as a prefix, so we need to let SCons know how to
# find those headers.
env.Append(CPPPATH=env['build_dir'] + '/include/Underworld')

# Keep a list of all the objects we build so we can make a library
# afterwards.
objs = []
suites = []
tst_exp = []
tst_input = []

# Process each directory uniformly.
dirs = Split('Rheology Utils libUnderworld')
for d in dirs:

    # Need the module name, which is just the directory.
    mod_name = env['ESCAPE']('"' + ''.join(d.split('/')) + '"')
    cpp_defs = [('CURR_MODULE_NAME', mod_name)] + env.get('CPPDEFINES', [])

    # Setup where to look for files.
    src_dir = d + '/src'
    inc_dir = 'include/Underworld/' + d
    tst_dir = d + '/tests'
    tst_exp_dir = tst_dir + '/expected'
    tst_input_dir = tst_dir + '/input'
    tst_install_dir = 'tests/Underworld/' + d

    # Install the headers and '.def' files.
    hdrs = env.Install(inc_dir, Glob(src_dir + '/*.h'))
    defs = env.Install(inc_dir, Glob(src_dir + '/*.def'))

    # Build our source files.
    srcs = Glob(src_dir + '/*.c')
    srcs = [s for s in srcs if s.path.find('-meta.c') == -1]
    objs += env.SharedObject(srcs, CPPDEFINES=cpp_defs)

    # Build any meta files.
    objs += env.stgSharedMeta(Glob(src_dir + '/*.meta'), CPPDEFINES=cpp_defs)

    # If we found any '.def' files make sure to register them as
    # explicit dependencies.
    if defs:
        env.Depends(hdrs + objs, defs)

    # Build any test suites we might find.
    suites += env.Object(Glob(tst_dir + '/*Suite.c'))

    # Install any test expected and input files
    tst_exp += env.Install(tst_install_dir + '/expected', Glob(tst_exp_dir + '/*'))
    tst_input += env.Install(tst_install_dir + '/input', Glob(tst_input_dir + '/*'))

# Need to install headers from libUnderworld.
env.Install('include/Underworld', Glob('libUnderworld/src/*.h'))

# Build libraries.
if env['shared_libs']:
    env.SharedLibrary('lib/Underworld', objs)

# Need to include the Underworld library for binaries.
libs = ['Underworld'] + env.get('LIBS', [])
# create Underworld executable -> it's just a copy of StGermain executable
env.Command('../bin/Underworld', '../bin/StGermain', 'cp $SOURCES $TARGET')

# Test runner program.
env.PCUTest('tests/testUnderworld', suites,
            PCU_LIBHEADERS="#include <StGermain/StGermain.h>\n#include <StgDomain/StgDomain.h>\n" \
                "#include <StgFEM/StgFEM.h>\n#include <PICellerator/PICellerator.h>\n" \
                "#include <Underworld/Underworld.h>", 
            PCU_SETUP="StGermain_Init(&argc, &argv);StgDomain_Init(&argc, &argv);" \
                "StgFEM_Init(&argc, &argv);PICellerator_Init(&argc, &argv);" \
                "Underworld_Init(&argc, &argv);",
            PCU_TEARDOWN="Underworld_Finalise();PICellerator_Finalise();StgFEM_Finalise();" \
                "StgDomain_Finalise();StGermain_Finalise();",
            LIBS=libs,
            PCU_EXP=tst_exp,
            PCU_INPUT=tst_input,
            PROJECT="Underworld")

# Build plugins.
dirs = ['libUnderworld/Toolbox',
        'plugins/EulerDeform',
        'plugins/IncompressibleExtensionBC',
        'plugins/Output/BuoyancyIntegrals',
        'plugins/DensityChangeAtDepth',
        'plugins/MaterialThermalDiffusivity',
        'plugins/MeshAdvectionCorrection',
        'plugins/VariableConditions/ShapeFemIC',
        'plugins/Output/Vrms',
        'plugins/Output/Mobility',
        'plugins/Output/Nusselt',
        'plugins/Output/VTKOutput',
        'plugins/Output/XDMFGenerator',
        'plugins/Output/MaxTemperature',
        'plugins/Output/MaxVelocity',
        'plugins/Output/BoundaryLayers',
        'plugins/Output/AverageTemperature',
        'plugins/ScalingChecks/Ra_Scaling',
        'SysTest/AnalyticPlugins/VelicIC',
        'SysTest/AnalyticPlugins/Velic_solA',
        'SysTest/AnalyticPlugins/Velic_solB',
        'SysTest/AnalyticPlugins/Velic_solCx',
        'SysTest/AnalyticPlugins/Velic_solHA',
        'SysTest/AnalyticPlugins/Velic_solKz',
        'SysTest/AnalyticPlugins/Velic_solS',
        'SysTest/AnalyticPlugins/LateralViscosityAnalytic',
        'SysTest/AnalyticPlugins/NonNewtonianShearSolution']
for d in dirs:

    name = 'Underworld_' + d.split('/')[-1] + 'module'
    mod_name = env['ESCAPE']('"' + ''.join(d.split('/')) + '"')
    cpp_defs = [('CURR_MODULE_NAME', mod_name)] + env.get('CPPDEFINES', [])

    env.Install('include/Underworld/' + d.split('/')[-1], Glob(d + '/*.h'))

    srcs = Glob(d + '/*.c')
    srcs = [s for s in srcs if s.path.find('-meta.c') == -1]
    objs = env.SharedObject(srcs, CPPDEFINES=cpp_defs)
    objs += env.stgSharedMeta(Glob(d + '/*.meta'), CPPDEFINES=cpp_defs)

    if env['shared_libs']:
        lib_pre = env['LIBPREFIXES']
        if not isinstance(lib_pre, list):
            lib_pre = [lib_pre]
        env.SharedLibrary('lib/' + name, objs,
                          SHLIBPREFIX='',
                          LIBPREFIXES=lib_pre + [''],
                          LIBS=libs)

# Install XML input files.
env.Install('lib/StGermain/Underworld', Glob('InputFiles/Underworld_Components/*.xml'))
env.Install('lib/StGermain/Underworld/BaseApps', Glob('InputFiles/BaseApps/*.xml'))
env.Install('lib/StGermain/Underworld/VariableConditions',
            Glob('InputFiles/VariableConditions/*.xml'))
env.Install('lib/StGermain/Underworld/Viewports', Glob('InputFiles/Viewports/*.xml'))
