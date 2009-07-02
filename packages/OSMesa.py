from config import Package

class OSMesa(Package):

    def gen_locations(self):
        yield ('/usr', ['/usr/include'], ['/usr/lib'])
        yield ('/usr/X11R6', ['/usr/X11R6/include'], ['/usr/X11R6/lib'])
        yield ('/usr/local', ['/usr/local/include'], ['/usr/local'])

    def gen_envs(self, loc):
        for env in Package.gen_envs(self, loc):
            env['pkg_headers'] = ['GL/glu.h', 'GL/osmesa.h']
            if self.find_libraries(loc[2], 'OSMesa'):
                env.PrependUnique(LIBS=['OSMesa', 'GLU'])
                yield env
