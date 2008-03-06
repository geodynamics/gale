import os
import SConfig

class OpenGL(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.header_sub_dir = 'GL'
        self.headers = [['gl.h', 'glu.h'],
                        ['OpenGL/gl.h', 'OpenGL/glu.h']] # For framework.
        self.libraries = [['GL', 'GLU']]
        self.frameworks = [['OpenGL']]
        self.have_define = 'HAVE_OPENGL'
