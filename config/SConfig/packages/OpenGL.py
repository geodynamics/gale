import os
import SConfig

class OpenGL(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False, **kw):
        SConfig.Package.__init__(self, scons_env, scons_opts, required, **kw)
        self.header_sub_dir = 'GL'
        self.headers = [['gl.h', 'glu.h']]
        self.libraries = [['GL', 'GLU']]
        self.frameworks = [['OpenGL']]
        self.have_define = 'HAVE_GL'
