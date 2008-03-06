import os
Import('env')

env = env.Copy()
env.AppendUnique(CPPPATH=[os.path.join(env['buildPath'], # Add PICellerator include path.
                                       'include',
                                       'PICellerator')])
env.project_name = 'PICellerator' # Need a project name.
env.clear_all() # ... so that our structures are ready.

env.build_directory('Voronoi')
env.build_directory('PopulationControl')
env.build_directory('Weights')
env.build_directory('MaterialPoints')
env.build_directory('Utils')

env.build_headers(env.glob('libPICellerator/src/*.h'), 'PICellerator')
env.build_objects(env.glob('libPICellerator/src/*.c'), 'libPICellerator')
env.build_objects(env.glob('libPICellerator/Toolbox/*.c'), 'Toolbox')
env.build_metadata(env.glob('libPICellerator/src/*.meta'), 'libPICellerator')
env.build_metadata(env.glob('libPICellerator/Toolbox/*.meta'), 'Toolbol')

env.build_library(env.get_hnodes(env.SharedObject), 'PICellerator')

env.build_library(env.get_hnodes(env.SharedObject, 'Toolbox'),
                  'PICellerator_Toolboxmodule', ['PICellerator'],
                  True)

env.build_tests(env.glob('libPICellerator/tests/test*.c'),
                'PICellerator', libs='PICellerator')

env.build_plugin('plugins/CalculateParticleDisplacement',
                 'CalculateParticleDisplacement')
env.build_plugin('plugins/Output/MaterialCentroid',
                 'Output/MaterialCentroid')

env.build_xmls(env.glob('Apps/src/*.xml'), 'StGermain/PICellerator')
