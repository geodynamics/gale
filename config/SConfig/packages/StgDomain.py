import os
import SConfig

class StgDomain(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Package.__init__(self, scons_env, scons_opts, required)
        self.dependency(SConfig.packages.StGermain)
        self.dependency(SConfig.packages.BlasLapack)
        self.dependency(SConfig.packages.HDF5, False)
        self.base_patterns = ['StgDomain*']
        self.headers = [[os.path.join('StgDomain', 'StgDomain.h')]]
        self.libraries = [['StgDomain']]
        self.symbols = [(['StgDomain_Init', 'StgDomain_Finalise'], '')]
        self.symbol_setup = '''MPI_Init(&argc, &argv);
StGermain_Init(&argc, &argv);
'''
        self.symbol_teardown = '''StGermain_Finalise();
MPI_Finalize();
'''
        self.symbol_calls = ['%s(&argc, &argv);', '%s();']
