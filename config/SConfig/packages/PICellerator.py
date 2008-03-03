import os
import SConfig

class PICellerator(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.setup_search_defaults()
        self.setup_options()
        self.dependencies = [SConfig.packages.StGermain]
        self.base_patterns = ['PICellerator*']
        self.headers = [os.path.join('PICellerator', 'PICellerator.h')]
        self.dependency_headers = [os.path.join('StGermain', 'StGermain.h'),
                                   os.path.join('StgDomain', 'StgDomain.h'),
                                   os.path.join('StgFEM', 'StgFEM.h')]
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
        self.require_shared = True
        self.use_rpath = True

    def get_run_error_message(self, console):
        if len(console):
            return 'Incompatible libraries, check \'config.log\'.'
