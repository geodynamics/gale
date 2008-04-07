import os
import SConfig

class PICellerator(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Package.__init__(self, scons_env, scons_opts, required)
        self.dependency(SConfig.packages.StgFEM)
        self.base_patterns = ['PICellerator*']
        self.headers = [[os.path.join('PICellerator', 'PICellerator.h')]]
        self.libraries = [['PICellerator']]
        self.symbols = [(['PICellerator_Init', 'PICellerator_Finalise'], '')]
        self.symbol_setup = '''MPI_Init(&argc, &argv);
StGermain_Init(&argc, &argv);
StgDomain_Init(&argc, &argv);
StgFEM_Init(&argc, &argv);'''
        self.symbol_teardown = '''StgFEM_Finalise();
StgDomain_Finalise();
StGermain_Finalise();
MPI_Finalize();'''
        self.symbol_calls = ['%s(&argc, &argv);', '%s();']
