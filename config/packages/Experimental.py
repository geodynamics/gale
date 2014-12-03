import os
from config import Package
from Underworld import Underworld

class Experimental(Package):

    def setup_dependencies(self):
        self.underworld = self.add_dependency(Underworld, required=True)

    def gen_locations(self):
        yield ('/usr', [], [])
        yield ('/usr/local', [], [])

    def gen_envs(self, loc):
        for env in Package.gen_envs(self, loc):
            self.headers = [os.path.join('StGermain', 'StGermain.h'),
                            os.path.join('StgDomain', 'StgDomain.h'),
                            os.path.join('StgFEM', 'StgFEM.h'),
                            os.path.join('PICellerator', 'PICellerator.h'),
                            os.path.join('Underworld', 'Underworld.h'),
                            os.path.join('Experimental', 'Experimental.h')]
            if self.find_libraries(loc[2], 'Underworld'):
                env.PrependUnique(LIBS=['Underworld'])
                yield env
