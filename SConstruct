import sys, os, subprocess

EnsureSConsVersion(0, 98)

#
# CUSTOMISE THE ENVIRONMENT HERE.
#

env = Environment(ENV=os.environ,
                  tools=['default', 'pcutest', 'stg', 'dist'],
                  toolpath=['StGermain/pcu/script', 'StGermain/script',
                            'script'])

# Needed for Darwin.
env['_abspath'] = lambda x: File(x).abspath

# Ludicrous-speed!
env.Decider("MD5-timestamp")

# Load configuration.
values = {}
execfile("config.cfg", globals(), values)
env._dict.update(values)

# Need to manipulate the build directory to keep SCons happy. Because of SCons' target
# rules we need to make the build directory a default target.
env["build_dir"] = os.path.join(env.GetLaunchDir(), env["build_dir"])
env["prefix"] = os.path.join(env.GetLaunchDir(), env["prefix"])
env["INST_BUILD_DIR"] = env["build_dir"]
env["INST_PREFIX"] = env["prefix"]
env.Default(env["build_dir"])

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
if not env['shared_libs']:
    env.AppendUnique(CPPDEFINES=['NOSHARED'])

# Need to extract some kind of hg version number.
subp = subprocess.Popen("hg identify",
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        shell=True)
out, err = subp.communicate()
res = subp.wait()
if res:
    print "Failed to extract hg revision number."
    hg_ver = 'unknown'
else:
    hg_ver = out.split()[0].strip()
env.AppendUnique(CPPDEFINES=[("VERSION", env["ESCAPE"]('"' + hg_ver + '"'))])

# Need to insert some 'HAVE_*' definitions based on what packages we
# found during configuration.
if 'HAVE_HDF5' in env['CPPDEFINES']:
    env.AppendUnique(CPPDEFINES=["READ_HDF5", "WRITE_HDF5"])

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
env.Prepend(LIBS=['pcu'])
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

if env['with_glucifer']:
    SConscript('gLucifer/SConscript',
               variant_dir=env['build_dir'] + '/gLucifer',
               duplicate=0)
    env.Prepend(LIBS=['glucifer'])

#
# Build static version of StGermain.
#

if env['static_libs']:
    env.Program('bin/StGermain',
                ['StGermain/src/main.c',
                 File(env['build_dir'] + '/StGermain/stg_static_modules.c').abspath])

# Adding in documentation.
env.Alias("doc", None, env.Action(File("StGermain/script/createDocs.py").abspath))
env.AlwaysBuild("doc")

# Dump package config.
filename = env['build_dir'] + '/lib/pkgconfig/stgermain.pc'
# env.write_pkgconfig(filename, 'StGermain', 'The StGermain Framework')

env.Dist("underworld-%s"%env.GetOption("dist_version"),
         ["configure.py", "SConstruct", "config", "script", "StGermain",
          "StgDomain", "StgFEM", "PICellerator", "Underworld",
          "Experimental", "gLucifer"])
