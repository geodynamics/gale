import platform
from config import Package

class SDL(Package):

    def setup_dependencies(self):
	    if platform.system() == "Darwin":
            self.cocoa = self.add_dependency(Cocoa, required=True)

    def gen_locations(self):
        yield ('/usr', ['/usr/include/SDL'], ['/usr/lib'])
        yield ('/usr/local', ['/usr/local'], ['/usr/local'])
        yield ('/usr/local', ['/usr/local/SDL'], ['/usr/local'])

    def gen_envs(self, loc):
        for env in Package.gen_envs(self, loc):
            env['pkg_headers'] = ['SDL.h']
            lib_env = env.Clone()
            lib_env.PrependUnique(LIBS=['SDL'])
            yield lib_env
            lib_env = env.Clone()
            lib_env.PrependUnique(LIBS=['SDL', 'SDLmain'])
            yield lib_env

	if platform.system() == "Darwin":
            env.AppendUnique(CPPPATH=['/System/Library/Frameworks/SDL.framework/Headers'])
            env.AppendUnique(FRAMEWORKS=['SDL'])
