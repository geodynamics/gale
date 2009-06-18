from Package import Package

class libFAME(Package):

    def gen_locations(self):
        yield ('/usr/local', ['/usr/local/include'], ['/usr/local/lib'])

    def gen_envs(self, loc):
        for env in Package.gen_envs(self, loc):
            env['pkg_headers'] = ['fame.h']
            if self.find_libraries(loc[2], 'fame'):
                env.PrependUnique(LIBS=['fame'])
                yield env
