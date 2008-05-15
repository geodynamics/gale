import os
Import('env')

env = env.Copy()
env.AppendUnique(CPPPATH=[env.get_build_path('include/Underworld')])

bases = ['Rheology', 'Utils']
src_objs = []
suite_hdrs = []
suite_objs = []
for base in bases:
    env.build_files(env.glob(base + '/src/*.def'), 'include/Underworld/' + base)
    env.build_headers(env.glob(base + '/src/*.h'), 'include/Underworld/' + base)
    src_objs += env.build_sources(env.glob(base + '/src/*.c'), 'Underworld/' + base)
    src_objs += env.build_metas(env.glob(base + '/src/*.meta'), 'Underworld/' + base)
    suite_hdrs += env.glob(base + '/tests/*Suite.h')
    suite_objs += env.build_sources(env.glob(base + '/tests/*Suite.c'), 'Underworld/' + base)

env.build_headers(env.glob('libUnderworld/src/*.h'), 'include/Underworld')
src_objs += env.build_sources(env.glob('libUnderworld/src/*.c'), 'Underworld/libUnderworld')
src_objs += env.build_sources(env.glob('libUnderworld/src/*.meta'), 'Underworld/libUnderworld')

# Build library.
if env['static_libraries']:
    env.Library(env.get_build_path('lib/Underworld'), src_objs)
if env['shared_libraries']:
    env.SharedLibrary(env.get_build_path('lib/Underworld'), src_objs)

# Build toolbox.
if env['shared_libraries']:
    objs = env.build_sources(env.glob('libUnderworld/Toolbox/*.c'),
                             'Underworld/libUnderworld/Toolbox')
    objs += env.build_metas(env.glob('libUnderworld/Toolbox/*.meta'),
                            'Underworld/libUnderworld/Toolbox')
    env.SharedLibrary(env.get_target_name('lib/Underworld_Toolboxmodule'), objs,
                      SHLIBPREFIX='',
                      LIBPREFIXES=[env['LIBPREFIXES']] + [''],
                      LIBS=['Underworld'] + env.get('LIBS', []))

# Build plugins.
if env['shared_libraries']:
    plgn_bases = ['EulerDeform', 'ExtractPetscObjects',
                  'IncompressibleExtensionBC', 'MaterialThermalDiffusivity',
                  'VariableConditions/ShapeTemperatureIC',
                  'Output/Vrms', 'Output/Nusselt']
    for base in plgn_bases:
        env.build_headers(env.glob('plugins/' + base + '/*.h'), 'include/Underworld/' + base.split('/')[-1])
        objs = env.build_sources(env.glob('plugins/' + base + '/*.c'), 'Underworld/' + base)
        name = 'Underworld_' + base.split('/')[-1] + 'module'
        env.SharedLibrary(env.get_build_path('lib/' + name), objs,
                          SHLIBPREFIX='',
                          LIBPREFIXES=[env['LIBPREFIXES']] + [''],
                          LIBS=['Underworld'] + env.get('LIBS', []))

# Build unit test runner.
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
runner_src = env.PCUSuiteRunner(env.get_build_path('Underworld/testUnderworld.c'), suite_hdrs)
runner_obj = env.SharedObject(runner_src)
env.Program(env.get_build_path('bin/testUnderworld'),
            runner_obj + suite_objs,
            LIBS=['Underworld', 'pcu'] + env.get('LIBS', []))

# Copy over XML files.
xml_bases = ['', 'BaseApps', 'VariableConditions', 'Viewports']
for base in xml_bases:
    dst = env.get_build_path('lib/StGermain/Underworld/' + base)
    for file in env.glob('InputFiles/src/' + base + '/*.xml'):
        env.Install(dst, file)
