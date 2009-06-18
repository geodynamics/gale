from Package import Package
from AGL import AGL

class Carbon(Package):

    def setup_dependencies(self):
        self.agl = self.add_dependency(AGL, required=False)

    def gen_envs(self, loc):
        env = self.env.Clone()
        env['pkg_headers'] = ['Carbon/Carbon.h']
        env.AppendUnique(FRAMEWORKS=['Carbon'])
        yield env
