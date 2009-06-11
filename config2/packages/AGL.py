from Package import Package

class Cocoa(Package):

    def gen_envs(self, loc):
        env = self.env.Clone()
        env['pkg_headers'] = ['AGL.h']
        env.AppendUnique(FRAMEWORKS=['AGL'])
        yield env
