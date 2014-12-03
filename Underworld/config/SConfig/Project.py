import os
import SCons.Script
import SConfig

class Project(SConfig.Node):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Node.__init__(self, scons_env, scons_opts, required)
        self.checks += [self.check_libs, self.print_results]

    def setup_options(self):
        self.opts.AddOptions(
            SCons.Script.BoolOption('with_debug',
                                    'Generate debugging symbols', 1),
            SCons.Script.BoolOption('static_libraries',
                                    'Build static libraries', 1),
            SCons.Script.BoolOption('shared_libraries',
                                    'Build shared libraries', 1),
            ('build_dir', 'Temporary build directory', 'build')
            )

    def check_libs(self):
        if not self.env['static_libraries'] and not self.env['shared_libraries']:
            self.ctx.Display("      Both static and shared libraries disabled!\n")
            return False
        return True

    def print_results(self):
        self.ctx.Display("  Static libraries: %s\n" % str(bool(self.env['static_libraries'])))
        self.ctx.Display("  Shared libraries: %s\n" % str(bool(self.env['shared_libraries'])))
        self.ctx.Display("  Using build directory: %s\n" % self.env['build_dir'])
        self.ctx.Display("  Debugging symbols: %s\n" % str(bool(self.env['with_debug'])))
        return True

    def setup(self):
        SConfig.Node.setup(self)
        modified = []
        if self.env['shared_libraries']:
            for pkg in self.env.package_list:
                if isinstance(pkg, SConfig.Package):
                    pkg.require_shared = True
                    modified += [pkg]
        return modified

    def enable(self, scons_env, old_state=None):
        SConfig.Node.enable(self, scons_env, old_state)

        # Setup debugging flags.
        if self.env['with_debug']:
            scons_env.MergeFlags('-g')

        # Setup the include paths.
        inc_dir = self.env.get_build_path('include')
        self.backup_variable(scons_env, 'CPPPATH', old_state)
        scons_env.PrependUnique(CPPPATH=[inc_dir])

        # Setup LIB_DIR.
        lib_dir = self.env.get_build_path('lib')
        self.backup_variable(scons_env, 'LIBPATH', old_state)
        scons_env.PrependUnique(LIBPATH=[lib_dir])

        # Setup the RPATH.
        self.backup_variable(scons_env, 'RPATH', old_state)
        scons_env.PrependUnique(RPATH=[scons_env.Dir(lib_dir).abspath])
