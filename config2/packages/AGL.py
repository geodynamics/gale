from Package import Package

class AGL(Package):

    def gen_envs(self, loc):
        env = self.env.Clone()
        env['pkg_headers'] = ['AGL/agl.h']
        env.AppendUnique(FRAMEWORKS=['AGL'])
        yield env
