Import('env')

env = env.Copy()
env.AppendUnique(CPPPATH=[os.path.join(env['buildPath'], # Add Underworld include path.
                                       'include',
                                       'Underworld')])
env.project_name = 'Underworld' # Need a project name.
env.clear_all() # ... so that our structures are ready.

env.build_directory('Rheology')
env.build_directory('Utils')

env.build_headers(env.glob('libUnderworld/src/*.h'), 'Underworld')
env.build_objects(env.glob('libUnderworld/src/*.c'), 'libUnderworld')
env.build_objects(env.glob('libUnderworld/Toolbox/*.c'), 'Toolbox')
env.build_metadata(env.glob('libUnderworld/src/*.meta'), 'libUnderworld')
env.build_metadata(env.glob('libUnderworld/Toolbox/*.meta'), 'Toolbol')

env.build_library(env.get_hnodes(env.SharedObject), 'Underworld')

env.build_library(env.get_hnodes(env.SharedObject, 'Toolbox'),
                  'Underworld_Toolboxmodule', ['Underworld'],
                  True)

env.build_tests(env.glob('libUnderworld/tests/test*.c'),
                'Underworld', libs='Underworld')

env.build_plugin('plugins/EulerDeform', 'EulerDeform')
env.build_plugin('plugins/ExtractPetscObjects', 'ExtractPetscObjects')
env.build_plugin('plugins/IncompressibleExtensionBC',
                 'IncompressibleExtensionBC')
env.build_plugin('plugins/MaterialThermalDiffusivity',
                 'MaterialThermalDiffusivity')
env.build_plugin('plugins/VariableConditions/ShapeTemperatureIC',
                 'VariableConditions/ShapeTemperatureIC')
env.build_plugin('plugins/Output/Vrms', 'Vrms')

env.build_xmls(env.glob('InputFiles/src/*.xml'), 'StGermain/Underworld')
env.build_xmls(env.glob('InputFiles/src/BaseApps/*.xml'),
               'StGermain/Underworld/BaseApps')
env.build_xmls(env.glob('InputFiles/src/VariableConditions/*.xml'),
               'StGermain/Underworld/VariableConditions')
env.build_xmls(env.glob('InputFiles/src/Viewports/*.xml'),
               'StGermain/Underworld/Viewports')
