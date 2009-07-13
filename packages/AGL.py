from config import Package

class AGL(Package):

    def gen_envs(self, loc):
        env = self.env.Clone()
        self.headers = ['AGL/agl.h']
        env.AppendUnique(FRAMEWORKS=['AGL'])
        yield env
