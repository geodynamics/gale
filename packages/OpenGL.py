import os, platform
from config import Package

class OpenGL(Package):

    def gen_locations(self):
        yield ('/usr', [], [])
        yield ('/usr', [], ['/usr/X11R6'])
        yield ('/usr/X11R6', [], [])
        yield ('/usr/local', [], [])

    def gen_base_extensions(self):
        for e in Package.gen_base_extensions(self):
            yield e
            yield ([os.path.join(h, 'GL') for h in e[0]], e[1])

    def gen_envs(self, loc):
        for env in Package.gen_envs(self, loc):
            self.headers = ['gl.h', 'glu.h']

	if platform.system() == "Darwin":
            env.AppendUnique(CPPPATH=['/System/Library/Frameworks/OpenGL.framework/Headers'])
            env.AppendUnique(FRAMEWORKS=['OpenGL'])
        else:
            env.PrependUnique(LIBS=['GL', 'GLU'])

        yield env
