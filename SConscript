import os
Import('env')

env = env.Copy()
env.AppendUnique(CPPPATH=[env.get_build_path('include/StgFEM')])

bases = ['Discretisation', 'SLE/LinearAlgebra', 'SLE/SystemSetup',
         'SLE/ProvidedSystems/AdvectionDiffusion', 'SLE/ProvidedSystems/Energy',
         'SLE/ProvidedSystems/StokesFlow', 'SLE/ProvidedSystems',
         'SLE', 'Assembly']
src_objs = []
suite_hdrs = []
suite_objs = []
for base in bases:
    env.build_files(env.glob(base + '/src/*.def'), 'include/StgFEM/' + base)
    env.build_headers(env.glob(base + '/src/*.h'), 'include/StgFEM/' + base)
    src_objs += env.build_sources(env.glob(base + '/src/*.c'), 'StgFEM/' + base)
    src_objs += env.build_metas(env.glob(base + '/src/*.meta'), 'StgFEM/' + base)
    suite_hdrs += env.glob(base + '/tests/*Suite.h')
    suite_objs += env.build_sources(env.glob(base + '/tests/*Suite.c'), 'StgFEM/' + base)

env.build_headers(env.glob('libStgFEM/src/*.h'), 'include/StgFEM')
src_objs += env.build_sources(env.glob('libStgFEM/src/*.c'), 'StgFEM/libStgFEM')
src_objs += env.build_sources(env.glob('libStgFEM/src/*.meta'), 'StgFEM/libStgFEM')

# Build library.
if env['static_libraries']:
    env.Library(env.get_build_path('lib/StgFEM'), src_objs)
if env['shared_libraries']:
    env.SharedLibrary(env.get_build_path('lib/StgFEM'), src_objs)

# Build toolbox.
if env['shared_libraries']:
    objs = env.build_sources(env.glob('libStgFEM/Toolbox/*.c'),
                             'StgFEM/libStgFEM/Toolbox')
    objs += env.build_metas(env.glob('libStgFEM/Toolbox/*.meta'),
                            'StgFEM/libStgFEM/Toolbox')
    env.SharedLibrary(env.get_target_name('lib/StgFEM_Toolboxmodule'), objs,
                      SHLIBPREFIX='',
                      LIBPREFIXES=env.make_list(env['LIBPREFIXES']) + [''],
                      LIBS=['StgFEM'] + env.get('LIBS', []))

# Build plugins.
if env['shared_libraries']:
    plgn_bases = ['plugins/CompareFeVariableAgainstReferenceSolution',
                  'plugins/Document',
                  'plugins/FeVariableImportExporters/FeVariable_ImportExport_ABAQUS',
                  'plugins/FeVariableImportExporters/FeVariable_ImportExport_SpecRidge2D',
                  'plugins/FileAnalyticSolution',
                  'plugins/Output/CPUTime',
                  'plugins/Output/FrequentOutput',
                  'plugins/Output/PeakMemory',
                  'plugins/Output/PrintFeVariableDiscreteValues',
                  'plugins/Output/PrintFeVariableDiscreteValues_2dBox',
                  'plugins/StandardConditionFunctions',
                  ('Apps/StokesMomentumUzawa/tests/LinearVelocityAnalytic',
                   'LinearVelocityAnalytic'),
                  ('Apps/StokesMomentumUzawa/tests/LidDrivenIsoviscousAnalytic',
                   'LidDrivenIsoviscousAnalytic'),
                  ('Apps/StokesMomentumUzawa/tests/SimpleShearAnalytic',
                   'SimpleShearAnalytic')]
    for base in plgn_bases:
        if isinstance(base, tuple):
            name = 'StgFEM_' + base[1] + 'module'
            base = base[0]
        else:
            name = 'StgFEM_' + base.split('/')[-1] + 'module'
        env.build_headers(env.glob(base + '/*.h'), 'include/StgFEM/' + base.split('/')[-1])
        objs = env.build_sources(env.glob(base + '/*.c'), 'StgFEM/' + base)
        env.SharedLibrary(env.get_build_path('lib/' + name), objs,
                          SHLIBPREFIX='',
                          LIBPREFIXES=env.make_list(env['LIBPREFIXES']) + [''],
                          LIBS=['StgFEM'] + env.get('LIBS', []))

# Build unit test runner.
env['PCURUNNERINIT'] = ''
env['PCURUNNERSETUP'] = """StGermain_Init( &argc, &argv );
   StgDomain_Init( &argc, &argv );
   StgFEM_Init( &argc, &argv );"""
env['PCURUNNERTEARDOWN'] = """StgFEM_Finalise();
   StgDomain_Finalise();
   StGermain_Finalise();"""
runner_src = env.PCUSuiteRunner(env.get_build_path('StgFEM/testStgFEM.c'), suite_hdrs)
runner_obj = env.SharedObject(runner_src)
env.Program(env.get_build_path('bin/testStgFEM'),
            runner_obj + suite_objs,
            LIBS=['StgFEM', 'pcu'] + env.get('LIBS', []))

# Copy over XML files.
xml_bases = ['']
for base in xml_bases:
    dst = env.get_build_path('lib/StGermain/StgFEM/' + base)
    for file in env.glob('Apps/StgFEM_Components/' + base + '/*.xml'):
        env.Install(dst, file)
