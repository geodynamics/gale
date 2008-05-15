import sys, os

sys.path.insert(0, os.path.abspath(os.path.join(os.getcwd(), 'config')))
import SConfig
SConscript('config/SConfig/SConscript')

#
# CUSTOMISE THE ENVIRONMENT HERE.
#

import platform
if platform.platform().find('ia64') != -1:
    # Hack needed for APAC, damn it all!
    env = Environment(ENV=os.environ, tools=['gcc', 'gnulink'])
    import SCons.Tool
    SCons.Tool.createStaticLibBuilder(env)
else:
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
    proj.dependency(SConfig.packages.HDF5, False, have_define='HAVE_HDF5')
    proj.dependency(SConfig.packages.PETSc, have_define='HAVE_PETSC')
    if env['with_glucifer']:
        proj.dependency(SConfig.packages.OpenGL, have_define='HAVE_GL')
        proj.dependency(SConfig.packages.OSMesa, False, have_define='HAVE_OSMESA')
        proj.dependency(SConfig.packages.X11, False, have_define='HAVE_X11')
        proj.dependency(SConfig.packages.SDL, False, have_define='HAVE_SDL')
        proj.dependency(SConfig.packages.libavcodec, False, have_define='HAVE_LIBAVCODEC')
        proj.dependency(SConfig.packages.libFAME, False, have_define='HAVE_LIBFAME')
        proj.dependency(SConfig.packages.libPNG, False, have_define='HAVE_PNG')
        proj.dependency(SConfig.packages.libJPEG, False, have_define='HAVE_JPEG')
        proj.dependency(SConfig.packages.libTIFF, False, have_define='HAVE_TIFF')
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
    env.Prepend(LIBS='Underworld')
    if env['with_glucifer']:
        SConscript('gLucifer/SConscript', exports='env')
