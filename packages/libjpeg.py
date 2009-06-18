from config import Package

class libJPEG(Package):

    def gen_locations(self):
        yield ('/usr/local', ['/usr/local/include'], ['/usr/local/lib'])

    def gen_envs(self, loc):
        for env in Package.gen_envs(self, loc):
            env['pkg_headers'] = ['jpeglib.h']
            if self.find_libraries(loc[2], 'jpeg'):
                env.PrependUnique(LIBS=['jpeg'])
                yield env
