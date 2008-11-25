import os
import SConfig

class OpenGL(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False, **kw):
        SConfig.Package.__init__(self, scons_env, scons_opts, required, **kw)
        self.base_dirs += ['/usr/X11R6']
        self.header_sub_dir = ['GL']
        self.headers = [['gl.h', 'glu.h']]
        self.libraries = [['GL', 'GLU']]
        self.frameworks = ['OpenGL']
        self.candidate_checks += [self.check_glx, self.check_mac]

        self.require_glx = kw.get('require_glx', False)
        self.x11 = self.dependency(SConfig.packages.X11, required=self.require_glx)

    def check_glx(self, cfg):
        """Check if this configuration comes with a glX implementation."""

        # For now we just have a silly header check.
        paths = cfg.inst.find_header('glx.h')
        cfg.has_glx = (len(paths) > 0)
        if cfg.has_glx:
            cfg.hdrs += ['glx.h']
            cfg.outputs += [('Has glx: %s', 'has_glx')]
            cfg.cpp_defines += ['HAVE_GLX']

        if self.require_glx and not paths:
            return False
        return True

    def check_mac(self, cfg):
        if self.platform.system == "Darwin":
            self.env.AppendUnique(SHLINKFLAGS=["-dylib_file /System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib:/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib"])
