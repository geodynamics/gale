import os, glob as pyglob
from SCons.Script.SConscript import SConsEnvironment

# Check versions of some things.
EnsurePythonVersion(2, 5)
EnsureSConsVersion(0, 97)

#
# Setup some basic path utilities.
#

# Globbing for hierarchical SCons scripts.
def glob(env, ptrn):
    if not os.path.isabs(ptrn):
        old = os.getcwd()
        os.chdir(Dir('.').srcnode().abspath)
        res = pyglob.glob(ptrn)
        os.chdir(old)
    else:
        res = pyglob.glob(ptrn)
    return res

# Normalise an SCons path to the project root.
def norm_path(env, path):
    cur_dir = env.Dir('.').srcnode().path
    if path[0] != '#':
        path = os.path.join(cur_dir, path)
    else:
        path = path[1:]
    return path

# Copy a file immediately.
def copy_file(env, dst, src):
    dst = env.norm_path(dst)
    old = os.getcwd()
    os.chdir(env.GetLaunchDir())
    if not File(dst).current():
        src = env.norm_path(src)
        dst_dir = os.path.dirname(dst)
        if not os.path.exists(dst_dir):
            Execute(Mkdir(dst_dir))
            Execute(Copy(dst, src))
        os.chdir(old)

# Add to the SCons main class.
SConsEnvironment.glob = glob
SConsEnvironment.norm_path = norm_path
SConsEnvironment.copy_file = copy_file

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
    PathOption('libxml2Dir', 'libXML2 installation path',
               None, PathOption.PathIsDir),
    PathOption('libxml2IncDir', 'libXML2 header installation path',
               None, PathOption.PathIsDir),
    PackageOption('csoap', 'Enable use of the package CSoap', 'no'))

#
# Create our substitution environment.
#

env = Environment(CC='cc', ENV=os.environ, options=opts)

# We need to add these guys so we don't get dictionary key errors.
env['CPPPATH'] = env['CPPPATH'] if 'CPPPATH' in env._dict else []
env['LIBPATH'] = env['LIBPATH'] if 'LIBPATH' in env._dict else []
env['CPPDEFINES'] = env['CPPDEFINES'] if 'CPPDEFINES' in env._dict else []
env['LIBS'] = env['LIBS'] if 'LIBS' in env._dict else []
env['RPATH'] = env['RPATH'] if 'RPATH' in env._dict else []

# Add any variables that get used throughout the whole build.
if env['debug']:
    env.Append(CCFLAGS='-g')
env.Append(CPPPATH=['#build/include'])
env.Append(LIBPATH=['#build/lib'])

# Setup any additional build targets.
env.Alias('install', env['prefix'])

#
# Configuration section.
#

SConscript('SConfigure', exports='env opts')

#
# Call target SConscripts.
#

env.BuildDir('build', '.', duplicate=0)
SConscript('build/SConscript', exports='env')
