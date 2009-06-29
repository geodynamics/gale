from config import Package

class Cocoa(Package):

    def gen_envs(self, loc):
        env = self.env.Clone()
        env.AppendUnique(FRAMEWORKS=['Cocoa'])
        yield env
