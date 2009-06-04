import sys, os, subprocess
import config

EnsureSConsVersion(0, 98)

#
# CUSTOMISE THE ENVIRONMENT HERE.
#

env = Environment(ENV=os.environ,
                  tools=['default', 'pcutest', 'stg', 'dist'],
                  toolpath=['StGermain/pcu/script', 'StGermain/script',
                            'script'])

# Ludicrous-speed!
env.Decider("MD5-timestamp")

# Load configuration.
values = {}
execfile("output.cfg", globals(), values)
env._dict.update(values)
env.MergeFlags(env.get("RAWFLAGS", ""))

# Need to manipulate the build directory to keep SCons happy. Because of SCons' target
# rules we need to make the build directory a default target.
env["build_dir"] = os.path.join(env.GetLaunchDir(), env["build_dir"])
env["prefix"] = os.path.join(env.GetLaunchDir(), env["prefix"])
env["INST_BUILD_DIR"] = env["build_dir"]
env["INST_PREFIX"] = env["prefix"]
env.Default(env["build_dir"])

# Scan command line targets to see if the user has specified any directories
# they want singled out to be built alone.<
env['dir_target'] = ARGUMENTS.get('target', '')
if len(env['dir_target']) and env['dir_target'][-1] == os.path.sep:
    env['dir_target'] = env['dir_target'][:-1]

# Add the build directory's include path.
env.AppendUnique(CPPPATH=env['build_dir'] + '/include')

# Need to define the extension for shared libraries as well
# as the library directory.
ext = env['ESCAPE']('"' + env['SHLIBSUFFIX'][1:] + '"')
lib_dir = env['ESCAPE']('"' + os.path.abspath(env['build_dir']) + '/lib' + '"')
env.AppendUnique(CPPDEFINES=[('MODULE_EXT', ext), ('LIB_DIR', lib_dir)])

# Include the library path.
env.AppendUnique(LIBPATH=env['build_dir'] + '/lib')
env.AppendUnique(RPATH=env.Dir(env['build_dir'] + '/lib').abspath)

# If we have no shared libraries, include a pre-processor definition to
# prevent modules from trying to load dynamically.
if not env['shared_libraries']:
    env.AppendUnique(CPPDEFINES=['NOSHARED', 'SINGLE_EXE'])

# Need to extract some kind of hg version number.
subp = subprocess.Popen("hg identify",
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        shell=True)
out, err = subp.communicate()
res = subp.wait()
if res:
    print "Failed to extract hg revision number."
    sys.exit()
env.AppendUnique(CPPDEFINES=[("VERSION", env["ESCAPE"]('"' + out.split()[0].strip() + '"'))])

# Need to insert some 'HAVE_*' definitions based on what packages we
# found during configuration.
if "OpenGL" in env["found_packages"]:
    env.AppendUnique(CPPDEFINES=["HAVE_GL"])
if "libpng" in env["found_packages"]:
    env.AppendUnique(CPPDEFINES=["HAVE_PNG"])
if "libjpeg" in env["found_packages"]:
    env.AppendUnique(CPPDEFINES=["HAVE_JPEG"])
if "libfame" in env["found_packages"]:
    env.AppendUnique(CPPDEFINES=["HAVE_LIBFAME"])
if "SDL" in env["found_packages"]:
    env.AppendUnique(CPPDEFINES=["HAVE_SDL"])
if "OSMesa" in env["found_packages"]:
    env.AppendUnique(CPPDEFINES=["HAVE_OSMESA"])
if "PETScExt" in env["found_packages"]:
    env.AppendUnique(CPPDEFINES=["HAVE_PETSCEXT"])
if "PETSc" in env["found_packages"]:
    env.AppendUnique(CPPDEFINES=["HAVE_PETSC"])
if "X11" in env["found_packages"]:
    env.AppendUnique(CPPDEFINES=["HAVE_X11"])
elif "Carbon" in env["found_packages"]:
    env.AppendUnique(CPPDEFINES=["HAVE_CARBON"])
if "HDF5" in env["found_packages"]:
    env.AppendUnique(CPPDEFINES=["READ_HDF5", "WRITE_HDF5"])
if "Python" in env["found_packages"]:
    env.AppendUnique(CPPDEFINES=["Python"])

#
# Make sure 'install' has a proper target.
#

env.Alias("install", env["prefix"])

#
# INSERT TARGETS HERE.
#

Export('env')

SConscript('StGermain/SConscript',
           variant_dir=env['build_dir'] + '/StGermain',
           duplicate=0)
env.Prepend(LIBS=['StGermain'])

SConscript('StgDomain/SConscript',
           variant_dir=env['build_dir'] + '/StgDomain',
           duplicate=0)
env.Prepend(LIBS=['StgDomain'])

SConscript('StgFEM/SConscript',
           variant_dir=env['build_dir'] + '/StgFEM',
           duplicate=0)
env.Prepend(LIBS=['StgFEM'])

SConscript('PICellerator/SConscript',
           variant_dir=env['build_dir'] + '/PICellerator',
           duplicate=0)
env.Prepend(LIBS=['PICellerator'])

SConscript('Underworld/SConscript',
           variant_dir=env['build_dir'] + '/Underworld',
           duplicate=0)
env.Prepend(LIBS=['Underworld'])

SConscript('gLucifer/SConscript',
           variant_dir=env['build_dir'] + '/gLucifer',
           duplicate=0)
env.Prepend(LIBS=['glucifer'])

# Adding in documentation.
env.Alias("doc", None, env.Action(File("config/tools/createDocs.py").abspath))
env.AlwaysBuild("doc")

# Dump package config.
filename = env['build_dir'] + '/lib/pkgconfig/stgermain.pc'
# env.write_pkgconfig(filename, 'StGermain', 'The StGermain Framework')

env.Dist("underworld-%s"%env.GetOption("dist_version"),
         ["configure.py", "SConstruct", "config", "script", "StGermain",
          "StgDomain", "StgFEM", "PICellerator", "Underworld",
          "gLucifer"])
