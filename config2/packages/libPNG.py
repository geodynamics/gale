from Package import Package

class libPNG(Package):

    def gen_locations(self):
        yield ('/usr/local', ['/usr/local/include'], ['/usr/local/lib'])

    def gen_envs(self, loc):
        for env in Package.gen_envs(self, loc):
            env['pkg_headers'] = ['png.h']
            if self.find_libraries(loc[2], 'png'):
                env.PrependUnique(LIBS=['png'])
                yield env
