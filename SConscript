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

# Need to include the StgFEM library for binaries/plugins.
libs = ['StgFEM'] + env.get('LIBS', [])

# Keep a list of all the objects we build so we can make a library
# afterwards.
objs = []
suites = []
tst_exp = []
tst_input = []

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
    tst_exp_dir = tst_dir + '/expected'
    tst_input_dir = tst_dir + '/input'
    tst_install_dir = 'tests/StgFEM/' + d

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

# Need to install headers from libStgFEM.
env.Install('include/StgFEM', Glob('libStgFEM/src/*.h'))

#
# Build plugins.
#

dirs = ['libStgFEM/Toolbox',
     'plugins/CompareFeVariableAgainstReferenceSolution',
     'plugins/Document',
     'plugins/FeVariableImportExporters/FeVariable_ImportExport_ABAQUS',
     'plugins/FeVariableImportExporters/FeVariable_ImportExport_SpecRidge2D',
#  'plugins/FileAnalyticSolution',
     'plugins/Output/CPUTime',
     'plugins/Output/FrequentOutput',
     'plugins/Output/FeVariableList',
     'plugins/Output/SwarmVariableList',
     'plugins/Output/PeakMemory',
     'plugins/Output/PrintFeVariableDiscreteValues',
     'plugins/Output/PrintFeVariableDiscreteValues_2dBox',
     'plugins/StandardConditionFunctions',
     'SysTest/AnalyticPlugins/CosineHillRotate',
     'Apps/TempDiffusion/tests/LinearTemperatureField',
     'Apps/StokesMomentumUzawa/tests/LinearVelocityAnalytic',
     'Apps/StokesMomentumUzawa/tests/LidDrivenIsoviscousAnalytic',
     'Apps/StokesMomentumUzawa/tests/SimpleShearAnalytic',
     'Apps/StokesMomentumUzawa/tests/LidDrivenStokesAnalytic']
pl_objs = []
pl_regs = []
for d in dirs:

    name = 'StgFEM_' + d.split('/')[-1] + 'module'
    mod_name = env['ESCAPE']('"' + ''.join(d.split('/')) + '"')
    cpp_defs = [('CURR_MODULE_NAME', mod_name)] + env.get('CPPDEFINES', [])

    env.Install('include/StgFEM/' + d.split('/')[-1], Glob(d + '/*.h'))

    srcs = Glob(d + '/*.c')
    srcs = [s for s in srcs if s.path.find('-meta.c') == -1]
    cur_objs = env.SharedObject(srcs, CPPDEFINES=cpp_defs)
    cur_objs += env.stgSharedMeta(Glob(d + '/*.meta'), CPPDEFINES=cpp_defs)

    # If we have shared libraries, build the dynamic plugin.
    if env['shared_libs']:
        lib_pre = env['LIBPREFIXES']
        if not isinstance(lib_pre, list):
            lib_pre = [lib_pre]
        env.SharedLibrary('lib/' + name, cur_objs,
                          SHLIBPREFIX='',
                          LIBPREFIXES=lib_pre + [''],
                          LIBS=libs)

    # If we are building static libs we need to construct a C file
    # registering each plugin.
    if env['static_libs']:
        pl_regs += [name[:-6] + '_Register(pluginsManager);']

    # Keep track of all the plugin objects.
    pl_objs += cur_objs

#
# Build libraries.
#

# Shared library.
if env['shared_libs']:
    env.SharedLibrary('lib/StgFEM', objs)

# Static library.
if env['static_libs']:
    # Need to include all the plugins aswell.
    l = env.StaticLibrary(env['build_dir'] + '/lib/StgFEM', objs + pl_objs)
    env.Install(env['prefix'] + '/lib', l)

"""
#
# Build static binary.
#

if env['static_libs']:

    # Select an appropriate name.
    if env['shared_libs']:
        name = 'StgFEM_static'
    else:
        name = 'StgFEM'

    # Build the registry code.
    reg_code = '#include <StGermain/StGermain.h>\n\n'
    reg_code += 'void register_plugins(PluginsManager *pluginsManager) {\n'
    reg_code += '  ' + '\n  '.join(pl_regs)
    reg_code += '\n}\n'
    reg_filename = os.path.join(env['build_dir'], 'StgFEM', 'reg.c')
    reg_file = open(reg_filename, 'w')
    reg_file.write(reg_code)
    reg_file.close()
    reg_obj = env.Object(reg_filename)

    # Build the static library.
    env.Program('bin/' + name, ['src/main_static.c', reg_obj])
"""

#
# Test runner program.
#

env.PCUTest('tests/testStgFEM', suites,
            PCU_LIBHEADERS="#include <StGermain/StGermain.h>\n#include <StgDomain/StgDomain.h>\n" \
                "#include <StgFEM/StgFEM.h>",
            PCU_SETUP="StGermain_Init(&argc, &argv);StgDomain_Init(&argc, &argv);" \
                "StgFEM_Init(&argc, &argv);",
            PCU_TEARDOWN="StgFEM_Finalise();StgDomain_Finalise();" \
                "StGermain_Finalise();",
            LIBS=libs,
            PCU_EXP=tst_exp,
            PCU_INPUT=tst_input,
            PROJECT="StgFEM")

# Install XML input files.
env.Install('lib/StGermain/StgFEM', Glob('Apps/StgFEM_Components/*.xml'))

#env.PCUSysTest('SysTest/RegressionTests/testAll.sh')
