import os, platform
import SConfig
import SCons.Script

class stgUnderworld(SConfig.Node):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Node.__init__(self, scons_env, scons_opts, required)
        self.dependency(SConfig.packages.cmath)
        self.dependency(SConfig.packages.libXML2)
        self.dependency(SConfig.packages.MPI)
        self.svnrevision = self.dependency(SConfig.packages.SVNRevision)
        self.dependency(SConfig.packages.BlasLapack)
        self.dependency(SConfig.packages.PETSc)
        self.dependency(SConfig.packages.HDF5, False)
        if self.env['with_glucifer']:
            self.dependency(SConfig.packages.OpenGL)
            self.dependency(SConfig.packages.OSMesa, False)
            self.dependency(SConfig.packages.X11, False)
            self.dependency(SConfig.packages.SDL, False)
            self.dependency(SConfig.packages.libavcodec, False)
            self.dependency(SConfig.packages.libFAME, False)
            self.dependency(SConfig.packages.libPNG, False)
            self.dependency(SConfig.packages.libJPEG, False)
            self.dependency(SConfig.packages.libTIFF, False)

    def setup_options(self):
        SConfig.Node.setup_options(self)
        self.opts.AddOptions(
            SCons.Script.BoolOption('debug',
                                    'Enable debugging', 1),
            SCons.Script.BoolOption('static_libraries',
                                    'Build static libraries', 0),
            SCons.Script.BoolOption('shared_libraries',
                                    'Build shared libraries', 1),
            SCons.Script.BoolOption('with_glucifer', 'Enable the gLucifer module', 1),
            ('buildPath', 'Temporary build path', 'build')
            )
        self.option_map = {'debug': None,
                           'static_libraries': None,
                           'shared_libraries': None,
                           'with_glucifer': None}

    def enable(self, scons_env):
        SConfig.Node.enable(self, scons_env)

        # Setup the build path.
        if not os.path.isabs(scons_env['buildPath']):
            scons_env['buildPath'] = '#' + scons_env['buildPath']

        # Setup LIB_DIR.
        lib_dir = os.path.join(scons_env['buildPath'], 'lib')
        abs_lib_dir = scons_env.Dir(lib_dir).abspath
        scons_env.AppendUnique(CPPDEFINES=[('LIB_DIR',
                                           scons_env['ESCAPE']('"' + abs_lib_dir + '"'))])
        scons_env.PrependUnique(LIBPATH=[lib_dir])

        # Setup the RPATH.
        scons_env.PrependUnique(RPATH=[scons_env.Dir(abs_lib_dir).abspath])

        # Setup the module extension.
        ext = scons_env['ESCAPE']('"' + scons_env['SHLIBSUFFIX'][1:] + '"')
        scons_env.Append(CPPDEFINES=[('MODULE_EXT', ext)])

        # Setup the include paths.
        inc_path = os.path.join(scons_env['buildPath'], 'include')
        scons_env.AppendUnique(CPPPATH=[inc_path])

        # Setup debugging.
        if scons_env['debug']:
            scons_env.MergeFlags('-g')

        # Setup the revision number.
        ver = scons_env['ESCAPE']('"' + str(self.svnrevision.revision) + '"')
        scons_env.Append(CPPDEFINES=[('VERSION', ver)])
