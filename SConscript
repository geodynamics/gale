Import('env')

env = env.Copy()
env.project_name = 'StgFEM' # Need a project name.
env.clear_all() # ... so that our structures are ready.

env.build_directory('Discretisation')
env.build_directory('SLE/LinearAlgebra')
env.build_directory('SLE/SystemSetup')
env.build_directory('SLE/ProvidedSystems/AdvectionDiffusion')
env.build_directory('SLE/ProvidedSystems/Energy')
env.build_directory('SLE/ProvidedSystems/StokesFlow')
env.build_directory('SLE/ProvidedSystems')
env.build_directory('SLE')
env.build_directory('Assembly')

env.build_headers(env.glob('libStgFEM/src/*.h'), 'StgFEM')
env.build_objects(env.glob('libStgFEM/src/*.c'), 'libStgFEM')
env.build_objects(env.glob('libStgFEM/Toolbox/*.c'), 'Toolbox')
env.build_metadata(env.glob('libStgFEM/src/*.meta'), 'libStgFEM')
env.build_metadata(env.glob('libStgFEM/Toolbox/*.meta'), 'Toolbol')

env.build_library(env.get_hnodes(env.SharedObject), 'StgFEM')

env.build_library(env.get_hnodes(env.SharedObject, 'Toolbox'),
                  'StgFEM_Toolboxmodule', ['StgFEM'],
                  True)

env.build_tests(env.glob('libStgFEM/tests/test*.c'),
                'StgFEM', libs='StgFEM')

env.build_plugin('plugins/CompareFeVariableAgainstReferenceSolution',
                 'CompareFeVariableAgainstReferenceSolution')
env.build_plugin('plugins/Document', 'Document')
env.build_plugin('plugins/FeVariableImportExporters/FeVariable_ImportExport_ABAQUS',
                 'FeVariableImportExporters/FeVariable_ImportExport_ABAQUS')
env.build_plugin('plugins/FeVariableImportExporters/FeVariable_ImportExport_SpecRidge2D',
                 'FeVariableImportExporters/FeVariable_ImportExport_SpecRidge2D')
env.build_plugin('plugins/FileAnalyticSolution', 'FileAnalyticSolution')
env.build_plugin('plugins/Output/CPUTime', 'CPUTime')
env.build_plugin('plugins/Output/FrequentOutput', 'FrequentOutput')
env.build_plugin('plugins/Output/PeakMemory', 'Output/PeakMemory')
env.build_plugin('plugins/Output/PrintFeVariableDiscreteValues',
                 'Output/PrintFeVariableDiscreteValues')
env.build_plugin('plugins/Output/PrintFeVariableDiscreteValues_2dBox',
                 'Output/PrintFeVariableDiscreteValues_2dBox')
env.build_plugin('plugins/StandardConditionFunctions',
                 'StandardConditionFunctions')

env.build_xmls(env.glob('Apps/src/*.xml'), 'StGermain/StgFEM')
