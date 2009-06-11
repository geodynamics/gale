from Package import Package

class HDF5(Package):

    def gen_locations(self):
        yield ('/usr/local', ['/usr/local'], ['/usr/local'])

    def gen_envs(self, loc):
        for env in Package.gen_envs(self, loc):
            env['pkg_headers'] = ['hdf5.h']
            env.PrependUnique(LIBS=['hdf5'])
            yield env
