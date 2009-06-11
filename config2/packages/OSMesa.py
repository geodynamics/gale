from Package import Package

class OSMesa(Package):

    def gen_locations(self):
        yield ('/usr', ['/usr/include/GL'], ['/usr/lib'])
        yield ('/usr/X11R6', ['/usr/X11R6/include'], ['/usr/X11R6/lib'])
        yield ('/usr/X11R6', ['/usr/X11R6/include/GL'], ['/usr/X11R6/lib'])
        yield ('/usr/local', ['/usr/local'], ['/usr/local'])
        yield ('/usr/local', ['/usr/local/GL'], ['/usr/local'])

    def gen_envs(self, loc):
        for env in Package.gen_envs(self, loc):
            env['pkg_headers'] = ['osmesa.h', 'glu.h']
            env.PrependUnique(LIBS=['OSMesa', 'GLU'])
            yield env
