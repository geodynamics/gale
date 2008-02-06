Import('env')

# Need to copy the environment for use here.
env = env.Copy()

# Add some extra stuff.
env.proj = 'PICellerator'
env.Append(CPPPATH=['#build/include/' + env.proj])

#
# Target specification section.
#

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
                  'PICellerator_Toolboxmodule',
                  True)

env.build_tests(env.glob('libPICellerator/tests/test*.c'),
                'PICellerator', libs='PICellerator')

env.build_plugin('plugins/CalculateParticleDisplacement',
                 'CalculateParticleDisplacement')
env.build_plugin('plugins/Output/MaterialCentroid',
                 'Output/MaterialCentroid')

env.build_xmls(env.glob('Apps/src/*.xml'), 'StGermain/PICellerator')
