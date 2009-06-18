import platform
from Package import Package

class OpenGL(Package):

    def gen_locations(self):
        yield ('/usr', ['/usr/include/GL'], ['/usr/lib'])
        yield ('/usr/X11R6', ['/usr/X11R6/include'], ['/usr/X11R6/lib'])
        yield ('/usr/X11R6', ['/usr/X11R6/include/GL'], ['/usr/X11R6/lib'])
        yield ('/usr/local', ['/usr/local'], ['/usr/local'])
        yield ('/usr/local', ['/usr/local/GL'], ['/usr/local'])

    def gen_envs(self, loc):
        for env in Package.gen_envs(self, loc):
            env['pkg_headers'] = ['gl.h', 'glu.h']
            env.PrependUnique(LIBS=['GL', 'GLU'])
            yield env

	if platform.system() == "Darwin":
            env.AppendUnique(CPPPATH=['/System/Library/Frameworks/OpenGL.framework/Headers'])
            env.AppendUnique(FRAMEWORKS=['OpenGL'])
            yield env
