import os, platform
import SConfig
import SCons.Script

class stgUnderworld(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.checks = [self.setup_environment]
        self.require(SConfig.packages.cmath)
        self.require(SConfig.packages.libXML2)
        self.require(SConfig.packages.MPI)
        self.require(SConfig.packages.SVNRevision)
        self.require(SConfig.packages.BlasLapack)
        self.require(SConfig.packages.PETSc)
        if self.env['with_glucifer']:
            self.require(SConfig.packages.OpenGL)
            mesa = self.require(SConfig.packages.OSMesa)
            mesa.required = False

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
        self.env.AppendUnique(CPPPATH=[os.path.join(inc_path, 'StGermain')])
        self.env.AppendUnique(CPPPATH=[os.path.join(inc_path, 'StgDomain')])
        self.env.AppendUnique(CPPPATH=[os.path.join(inc_path, 'StgFEM')])
        self.env.AppendUnique(CPPPATH=[os.path.join(inc_path, 'PICellerator')])
        self.env.AppendUnique(CPPPATH=[os.path.join(inc_path, 'Underworld')])

        # Setup debugging.
        if self.env['debug']:
            self.env.MergeFlags('-g')

        # Setup 64 bit builds.
        #if platform.architecture()[0].find('64') != -1:
        #    self.env.MergeFlags('-m64')

        return [1, '', '']
