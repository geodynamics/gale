import os, glob

# Check versions of some things.
EnsurePythonVersion(2, 5)
EnsureSConsVersion(0, 97)

#
# Setup our option database.
#

opts = Options('config.cache')
opts.AddOptions(
    BoolOption('debug', 'Enable debugging version', 1),
    BoolOption('useMpiRecord', 'Don''t know what this does...', 0),
    PathOption('prefix', 'Installation path',
               '/usr/local', PathOption.PathIsDirCreate),
    PathOption('mpichDir', 'MPI installation path',
               None, PathOption.PathIsDir),
    PathOption('mpichIncDir', 'MPICH header installation path',
               None, PathOption.PathIsDir),
    PathOption('mpichLibDir', 'MPICH library installation path',
               None, PathOption.PathIsDir),
    PathOption('petscDir', 'PETSc installation path',
               None, PathOption.PathIsDir),
    PathOption('libxml2Dir', 'libXML2 installation path',
               None, PathOption.PathIsDir),
    PathOption('libxml2IncDir', 'libXML2 header installation path',
               None, PathOption.PathIsDir),
    PackageOption('csoap', 'Enable use of the package CSoap', 'no'))

#
# Create our substitution environment.
#

env = Environment(CC='cc', ENV=os.environ, options=opts)
env['CPPPATH'] = env['CPPPATH'] if 'CPPPATH' in env._dict else []
env['LIBPATH'] = env['LIBPATH'] if 'LIBPATH' in env._dict else []
env['CPPDEFINES'] = env['CPPDEFINES'] if 'CPPDEFINES' in env._dict else []
env['LIBS'] = env['LIBS'] if 'LIBS' in env._dict else []
env['RPATH'] = env['RPATH'] if 'RPATH' in env._dict else []

# Add any variables that get used throughout the whole build.
env.proj = 'StgFEM'
if env['debug']:
    env.Append(CCFLAGS='-g')
env.Append(CPPPATH=['#build/include'])
env.Append(CPPPATH=['#build/include/' + env.proj])
env.Append(LIBPATH=['#build/lib'])

# Add any helper functions we may need.
SConscript('StgSCons', exports='env')

# Setup any additional builds.
env.Alias('install', env['prefix'])

#
# Configuration section.
#

SConscript('SConfigure', exports='env opts')

#
# Target specification section.
#

env.BuildDir('build', '.', duplicate=0)

env.build_directory('Discretisation')
env.build_directory('SLE/LinearAlgebra')
env.build_directory('SLE/SystemSetup')
env.build_directory('SLE/ProvidedSystems/AdvectionDiffusion')
env.build_directory('SLE/ProvidedSystems/Energy')
env.build_directory('SLE/ProvidedSystems/StokesFlow')
env.build_directory('SLE/ProvidedSystems')
env.build_directory('SLE')
env.build_directory('Assembly')

env.build_headers(glob.glob('libStgFEM/src/*.h'), 'StgFEM')
env.build_objects(glob.glob('libStgFEM/src/*.c'), 'libStgFEM')
env.build_objects(glob.glob('libStgFEM/Toolbox/*.c'), 'Toolbox')
env.build_metadata(glob.glob('libStgFEM/src/*.meta'), 'libStgFEM')
env.build_metadata(glob.glob('libStgFEM/Toolbox/*.meta'), 'Toolbol')

env.build_library(env.obj_nodes['.'], 'StgFEM')

env.build_library(env.obj_nodes['Toolbox']['.'],
                  'StgFEM_Toolboxmodule',
                  True)

env.build_tests(glob.glob('libStgFEM/tests/test*.c'),
                'StgFEM', libs='StgFEM')

env.build_plugin('plugins/CompareFeVariableAgainstReferenceSolution',
                 'CompareFeVariableAgainstReferenceSolution')
env.build_plugin('plugins/Document', 'Document')
env.build_plugin('plugins/FeVariableImportExporters/FeVariable_ImportExport_ABAQUS',
                 'FeVariableImportExporters/FeVariable_ImportExport_ABAQUS')
env.build_plugin('plugins/FeVariableImportExporters/FeVariable_ImportExport_SpecRidge2D',
                 'FeVariableImportExporters/FeVariable_ImportExport_SpecRidge2D')
env.build_plugin('plugins/FileAnalyticSolution', 'FileAnalyticSolution')
env.build_plugin('plugins/Output/CPUTime', 'Output/CPUTime')
env.build_plugin('plugins/Output/FrequentOutput', 'Output/FrequentOutput')
env.build_plugin('plugins/Output/PeakMemory', 'Output/PeakMemory')
env.build_plugin('plugins/Output/PrintFeVariableDiscreteValues',
                 'Output/PrintFeVariableDiscreteValues')
env.build_plugin('plugins/Output/PrintFeVariableDiscreteValues_2dBox',
                 'Output/PrintFeVariableDiscreteValues_2dBox')
env.build_plugin('plugins/StandardConditionFunctions',
                 'StandardConditionFunctions')
