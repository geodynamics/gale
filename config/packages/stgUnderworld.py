import os, platform
import SConfig
import SCons.Script

class stgUnderworld(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options, True)
        self.checks = [self.setup_environment]
        self.dependency(SConfig.packages.cmath)
        self.dependency(SConfig.packages.libXML2)
        self.dependency(SConfig.packages.MPI)
        self.dependency(SConfig.packages.SVNRevision)
        self.dependency(SConfig.packages.BlasLapack)
        self.dependency(SConfig.packages.PETSc)
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
        SConfig.Package.setup_options(self)
        if self.opts:
            self.opts.AddOptions(
                SCons.Script.BoolOption('debug',
                                        'Enable debugging', 1),
                SCons.Script.BoolOption('staticLibraries',
                                        'Build static libraries only', 0),
                SCons.Script.BoolOption('with_glucifer', 'Enable the gLucifer module', 1),
                ('buildPath', 'Temporary build path', 'build'))

    def setup_environment(self):
        # Setup the build path.
        if not os.path.isabs(self.env['buildPath']):
            self.env['buildPath'] = '#' + self.env['buildPath']

        # Setup LIB_DIR.
        lib_dir = os.path.join(self.env['buildPath'], 'lib')
        abs_lib_dir = self.env.Dir(lib_dir).abspath
        self.env.AppendUnique(CPPDEFINES=[('LIB_DIR',
                                           self.env['ESCAPE']('"' + abs_lib_dir + '"'))])
        self.env.PrependUnique(LIBPATH=[lib_dir])

        # Setup the RPATH.
        self.env.PrependUnique(RPATH=[self.env.Dir(abs_lib_dir).abspath])

        # Setup the module extension.
        ext = self.env['ESCAPE']('"' + self.env['SHLIBSUFFIX'][1:] + '"')
        self.env.Append(CPPDEFINES=[('MODULE_EXT', ext)])

        # Setup the include paths.
        inc_path = os.path.join(self.env['buildPath'], 'include')
        self.env.AppendUnique(CPPPATH=[inc_path])

        # Setup debugging.
        if self.env['debug']:
            self.env.MergeFlags('-g')

        return [1, '', '']
