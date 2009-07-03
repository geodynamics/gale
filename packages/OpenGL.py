import platform
from config import Package

class OpenGL(Package):

    def gen_locations(self):
        yield ('/usr', ['/usr/include/GL'], ['/usr/lib'])
        yield ('/usr', ['/usr/include/GL'], ['/usr/X11R6/lib'])
        yield ('/usr', ['/usr/include/GL'], ['/usr/X11R6/lib64'])
        yield ('/usr/X11R6', ['/usr/X11R6/include'], ['/usr/X11R6/lib'])
        yield ('/usr/X11R6', ['/usr/X11R6/include/GL'], ['/usr/X11R6/lib'])
        yield ('/usr/X11R6', ['/usr/X11R6/include'], ['/usr/X11R6/lib64'])
        yield ('/usr/X11R6', ['/usr/X11R6/include/GL'], ['/usr/X11R6/lib64'])
        yield ('/usr/local', ['/usr/local/include'], ['/usr/local/lib'])
        yield ('/usr/local', ['/usr/local/include/GL'], ['/usr/local/lib'])
        yield ('/usr/local', ['/usr/local/include'], ['/usr/local/lib64'])
        yield ('/usr/local', ['/usr/local/include/GL'], ['/usr/local/lib64'])

    def gen_envs(self, loc):
        for env in Package.gen_envs(self, loc):
            env['pkg_headers'] = ['gl.h', 'glu.h']

	if platform.system() == "Darwin":
            env.AppendUnique(CPPPATH=['/System/Library/Frameworks/OpenGL.framework/Headers'])
            env.AppendUnique(FRAMEWORKS=['OpenGL'])
        else:
            env.PrependUnique(LIBS=['GL', 'GLU'])

        yield env
