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
    PathOption('mpiDir', 'MPI installation path',
               None, PathOption.PathIsDir),
    PathOption('petscDir', 'PETSc installation path',
               None, PathOption.PathIsDir),
    PathOption('libxml2Dir', 'libXML2 installation path',
               None, PathOption.PathIsDir),
    PathOption('libxml2IncDir', 'libXML2 installation path',
               None, PathOption.PathIsDir),
    PackageOption('csoap', 'Enable use of the package CSoap', 'no'))

#
# Create our substitution environment.
#

env = Environment(CC='cc', ENV=os.environ, options=opts)
if env['debug']:
    env.Append(CCFLAGS='-g')
env.Append(CPPPATH=['#build/include'])
env.Append(CPPPATH=['#build/include/StGermain'])
env.Append(LIBPATH=['#build/lib'])
env.Alias('install', env['prefix'])

# Add any variables that get used throughout the whole build.
env.proj = 'StGermain'

# Add any helper functions we may need.
SConscript('StgSCons', exports='env')

#
# Configuration section.
#

SConscript('SConfigure', exports='env opts')

#
# Target specification section.
#

# Specify build and source directories.
env.BuildDir('build', '.', duplicate=0)

# Recurse into subdirectories to build source.
objs = SConscript('build/Regresstor/SConscript', exports='env')
objs += SConscript('build/Base/SConscript', exports='env')
objs += SConscript('build/Utils/SConscript', exports='env')

# Build the StGermain library.
SConscript('build/libStGermain/SConscript', exports='env objs')
SConscript('build/src/SConscript', exports='env')
