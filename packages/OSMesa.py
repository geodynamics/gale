from config import Package

class OSMesa(Package):

    def gen_locations(self):
        yield ('/usr', ['/usr/include/GL'], ['/usr/lib'])
        yield ('/usr', ['/usr/include/GL'], ['/usr/lib64'])
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
            env['pkg_headers'] = ['glu.h', 'osmesa.h']
            if self.find_libraries(loc[2], 'OSMesa'):
                env.PrependUnique(LIBS=['OSMesa', 'GLU'])
                yield env
