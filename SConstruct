import sys, os

sys.path.insert(0, os.path.abspath(os.path.join(os.getcwd(), 'config')))
import SConfig

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
SConscript('config/SConfig/SConscript', exports='env')

# Determine whether we are configuring, helping or building.
if 'config' in COMMAND_LINE_TARGETS or 'help' in COMMAND_LINE_TARGETS:

    #
    # INSERT CONFIGURATION HERE.
    #

    proj = env.Package(SConfig.Project)
    proj.opts.AddOptions(
        BoolOption('with_glucifer', 'Use gLucifer visualisation', 1),
        BoolOption('read_hdf5', 'Read from HDF5 files when restarting from checkpoint', 0),
        BoolOption('write_hdf5', 'Write checkpoint files as HDF5', 0)
        )
    proj.dependency(SConfig.packages.cmath)
    proj.dependency(SConfig.packages.libXML2)
    proj.dependency(SConfig.packages.MPI)
    proj.dependency(SConfig.packages.HGRevision)
    proj.dependency(SConfig.packages.BlasLapack)
    
    if env['read_hdf5'] and env['write_hdf5']:
        proj.dependency(SConfig.packages.HDF5, False, have_define=['READ_HDF5','WRITE_HDF5'])
    elif env['write_hdf5']:
        proj.dependency(SConfig.packages.HDF5, False, have_define='WRITE_HDF5')
    elif env['read_hdf5']:
        proj.dependency(SConfig.packages.HDF5, False, have_define='READ_HDF5')
    
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
        if platform.system() == 'Darwin':
          proj.dependency(SConfig.packages.CoreServices)
          proj.dependency(SConfig.packages.Carbon, False, have_define='HAVE_CARBON')
          proj.dependency(SConfig.packages.AGL, False, have_define='HAVE_AGL')
    env.configure_packages()

    # Need to define the extension for shared libraries as well
    # as the library directory.
    ext = env['ESCAPE']('"' + env['SHLIBSUFFIX'][1:] + '"')
    lib_dir = env['ESCAPE']('"' + os.path.abspath(env['build_dir']) + '/lib' + '"')
    env.AppendUnique(CPPDEFINES=[('MODULE_EXT', ext), ('LIB_DIR', lib_dir)])

    # If we have no shared libraries, include a pre-processor definition to
    # prevent modules from trying to load dynamically.
    if not env['shared_libraries']:
        env.AppendUnique(CPPDEFINES=['NOSHARED', 'SINGLE_EXE'])

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

    module = SConscript('StGermain/SConscript', exports='env')
    module_proto = module[0]
    module_code = module[1]
    env.Prepend(LIBS='StGermain')
    module = SConscript('StgDomain/SConscript', exports='env')
    module_proto += module[0]
    module_code += module[1]
    env.Prepend(LIBS='StgDomain')
    module = SConscript('StgFEM/SConscript', exports='env')
    module_proto += module[0]
    module_code += module[1]
    env.Prepend(LIBS='StgFEM')
    module = SConscript('PICellerator/SConscript', exports='env')
    module_proto += module[0]
    module_code += module[1]
    env.Prepend(LIBS='PICellerator')
    module = SConscript('Underworld/SConscript', exports='env')
    module_proto += module[0]
    module_code += module[1]
    env.Prepend(LIBS='Underworld')
    if env['with_glucifer']:
        module += SConscript('gLucifer/SConscript', exports='env')
        env.Prepend(LIBS='glucifer')
        module_proto += module[0]
        module_code += module[1]

    # If we're building only with static libraries, we need to build the
    # StGermain executable here, after we've finished building all the
    # other static libraries.
    if not env['shared_libraries'] and env['static_libraries']:
        module_src = env.build_module_register(env.get_build_path('StGermain/ModuleRegister.c'), module_proto, module_code)
        module_obj = env.SharedObject(module_src)
        module_lib = env.Library(env.get_build_path('lib/ModuleRegister'), module_obj)
        libs = env.get('LIBS', [])
        ind = libs.index('StGermain')
        libs.insert(ind + 1, module_lib)

        src = 'StGermain/src/main.c'
        obj_dst = env.get_build_path('StGermain/' + src[:-2])
        obj = env.SharedObject(obj_dst, src)
        prog_dst = env.get_build_path('bin/StGermain')
        env.Program(prog_dst, obj, LIBS=libs)

        src = 'StGermain/Base/FlattenXML/src/main.c'
        obj_dst = env.get_build_path('StGermain/' + src[:-2])
        obj = env.SharedObject(obj_dst, src)
        prog_dst = env.get_build_path('bin/FlattenXML')
        env.Program(prog_dst, obj, LIBS=libs)

    # Dump package config.
    filename = env.get_build_path('lib/pkgconfig/stgermain.pc')
    env.write_pkgconfig(filename, 'StGermain', 'The StGermain Framework')
