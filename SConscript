import os
Import('env')

# Need to make a copy because SCons uses the environment
# at it's final state, so StGermain ends up depending on
# StgDomain, etc.
env = env.Clone()

# Inside each project we will be accessing headers without the
# project name as a prefix, so we need to let SCons know how to
# find those headers.
env.Append(CPPPATH=env['build_dir'] + '/include/StgFEM')

# Keep a list of all the objects we build so we can make a library
# afterwards.
objs = []
suites = []

# Process each directory uniformly.
dirs = Split('Discretisation SLE/SystemSetup SLE/ProvidedSystems/AdvectionDiffusion ' \
                 'SLE/ProvidedSystems/AdvDiff_CharGalerkin SLE/ProvidedSystems/Energy ' \
                 'SLE/ProvidedSystems/StokesFlow SLE/ProvidedSystems SLE ' \
                 'Assembly libStgFEM')
for d in dirs:

    # Need the module name, which is just the directory.
    mod_name = env['ESCAPE']('"' + ''.join(d.split('/')) + '"')
    cpp_defs = [('CURR_MODULE_NAME', mod_name)] + env.get('CPPDEFINES', [])

    # Setup where to look for files.
    src_dir = d + '/src'
    inc_dir = 'include/StgFEM/' + d
    tst_dir = d + '/tests'

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

# Need to install headers from libStgFEM.
env.Install('include/StgFEM', Glob('libStgFEM/src/*.h'))

# Build libraries.
if env['shared_libraries']:
    env.SharedLibrary('lib/StgFEM', objs)

# Need to include the StgFEM library for binaries.
libs = ['StgFEM'] + env.get('LIBS', [])

# Test runner program.
env.PCUTest('tests/testStgFEM', suites,
            PCU_SETUP="StGermain_Init(&argc, &argv);StgDomain_Init(&argc, &argv);" \
                "StgFEM_Init(&argc, &argv);",
            PCU_TEARDOWN="StgFEM_Finalise();StgDomain_Finalise();" \
                "StGermain_Finalise();",
            LIBS=libs)

# Build plugins.
dirs = ['libStgFEM/Toolbox',
     'plugins/CompareFeVariableAgainstReferenceSolution',
     'plugins/Document',
     'plugins/FeVariableImportExporters/FeVariable_ImportExport_ABAQUS',
     'plugins/FeVariableImportExporters/FeVariable_ImportExport_SpecRidge2D',
     'plugins/FileAnalyticSolution',
     'plugins/Output/CPUTime',
     'plugins/Output/FrequentOutput',
     'plugins/Output/FeVariableList',
     'plugins/Output/SwarmVariableList',
     'plugins/Output/PeakMemory',
     'plugins/Output/PrintFeVariableDiscreteValues',
     'plugins/Output/PrintFeVariableDiscreteValues_2dBox',
     'plugins/StandardConditionFunctions',
     'Apps/StokesMomentumUzawa/tests/LinearVelocityAnalytic',
     'Apps/StokesMomentumUzawa/tests/LidDrivenIsoviscousAnalytic',
     'Apps/StokesMomentumUzawa/tests/SimpleShearAnalytic',
     'Apps/StokesMomentumUzawa/tests/LidDrivenStokesAnalytic']
for d in dirs:

    name = 'StgFEM_' + d.split('/')[-1] + 'module'
    mod_name = env['ESCAPE']('"' + ''.join(d.split('/')) + '"')
    cpp_defs = [('CURR_MODULE_NAME', mod_name)] + env.get('CPPDEFINES', [])

    env.Install('include/StgFEM/' + d.split('/')[-1],
                Glob(d + '/*.h'))

    srcs = Glob(d + '/*.c')
    srcs = [s for s in srcs if s.path.find('-meta.c') == -1]
    objs = env.SharedObject(srcs, CPPDEFINES=cpp_defs)
    objs += env.stgSharedMeta(Glob(d + '/*.meta'), CPPDEFINES=cpp_defs)

    if env['shared_libraries']:
        lib_pre = env['LIBPREFIXES']
        if not isinstance(lib_pre, list):
            lib_pre = [lib_pre]
        env.SharedLibrary('lib/' + name, objs,
                          SHLIBPREFIX='',
                          LIBPREFIXES=lib_pre + [''],
                          LIBS=libs)

# Install XML input files.
env.Install('lib/StGermain/StgFEM', Glob('Apps/StgFEM_Components/*.xml'))
