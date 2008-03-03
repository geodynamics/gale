import os
import SConfig

class StgDomain(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.setup_search_defaults()
        self.setup_options()
        self.dependencies = [SConfig.packages.StGermain]
        self.base_patterns = ['StgDomain*']
        self.headers = [os.path.join('StgDomain', 'StgDomain.h')]
        self.dependency_headers = [os.path.join('StGermain', 'StGermain.h')]
        self.libraries = [['StgDomain']]
        self.symbols = [(['StgDomain_Init', 'StgDomain_Finalise'], '')]
        self.symbol_setup = '''MPI_Init(&argc, &argv);
StGermain_Init(&argc, &argv);
'''
        self.symbol_teardown = '''StGermain_Finalise();
MPI_Finalize();
'''
        self.symbol_calls = ['%s(&argc, &argv);', '%s();']
        self.require_shared = True
        self.use_rpath = True

    def get_run_error_message(self, console):
        if len(console):
            return 'Incompatible libraries, check \'config.log\'.'
