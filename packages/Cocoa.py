from config import Package

class Cocoa(Package):

    def gen_envs(self, loc):
        env = self.env.Clone()
        env['pkg_headers'] = ['Cocoa/Cocoa.h']
        env.AppendUnique(FRAMEWORKS=['Cocoa'])
        yield env
