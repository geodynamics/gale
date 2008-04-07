import sys, os

sys.path.insert(0, os.path.abspath(os.path.join(os.getcwd(), 'config')))
import SConfig
SConscript('config/SConfig/SConscript')

#
# CUSTOMISE THE ENVIRONMENT HERE.
#

env = Environment(ENV=os.environ)
env['_abspath'] = lambda x: File(x).abspath # Needed by Darwin.

# Determine whether we are configuring, helping or building.
if 'config' in COMMAND_LINE_TARGETS or 'help' in COMMAND_LINE_TARGETS:

    #
    # INSERT CONFIGURATION HERE.
    #

    proj = env.Package(SConfig.Project)
    proj.opts.AddOptions(
        BoolOption('with_glucifer', 'Use gLucifer visualisation', 1)
        )
    proj.dependency(SConfig.packages.cmath)
    proj.dependency(SConfig.packages.libXML2)
    proj.dependency(SConfig.packages.MPI)
    proj.dependency(SConfig.packages.SVNRevision)
    proj.dependency(SConfig.packages.BlasLapack)
    proj.dependency(SConfig.packages.HDF5, False)
    petsc = proj.dependency(SConfig.packages.PETSc)
    petsc.have_define = 'HAVE_PETSC'
    if env['with_glucifer']:
        proj.dependency(SConfig.packages.OpenGL)
        proj.dependency(SConfig.packages.OSMesa, False)
        proj.dependency(SConfig.packages.X11, False)
        proj.dependency(SConfig.packages.SDL, False)
        proj.dependency(SConfig.packages.libavcodec, False)
        proj.dependency(SConfig.packages.libFAME, False)
        proj.dependency(SConfig.packages.libPNG, False)
        proj.dependency(SConfig.packages.libJPEG, False)
        proj.dependency(SConfig.packages.libTIFF, False)
    env.configure_packages()

    # Need to define the extension for shared libraries as well
    # as the library directory.
    ext = env['ESCAPE']('"' + env['SHLIBSUFFIX'][1:] + '"')
    lib_dir = env['ESCAPE']('"' + os.path.abspath(env['build_dir']) + '/lib' + '"')
    env.AppendUnique(CPPDEFINES=[('MODULE_EXT', ext), ('LIB_DIR', lib_dir)])

    # Save results.
    env.save_config()

else:
    # Load configuration.
    env.load_config()
    SConscript('StGermain/pcu/script/scons.py', exports='env')
    SConscript('StGermain/script/scons.py', exports='env')
    env.Append(CPPDEFINES=['HAVE_SCONS']) # Needed for old VMake.

    #
    # INSERT TARGETS HERE.
    #

    SConscript('StGermain/SConscript', exports='env')
    env.Prepend(LIBS='StGermain')
    SConscript('StgDomain/SConscript', exports='env')
    env.Prepend(LIBS='StgDomain')
    SConscript('StgFEM/SConscript', exports='env')
    env.Prepend(LIBS='StgFEM')
    SConscript('PICellerator/SConscript', exports='env')
    env.Prepend(LIBS='PICellerator')
    SConscript('Underworld/SConscript', exports='env')
    if env['with_glucifer']:
        SConscript('gLucifer/SConscript', exports='env')
