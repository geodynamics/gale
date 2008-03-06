import os
import SConfig

class StGermain(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.require(SConfig.packages.cmath)
        self.require(SConfig.packages.libXML2)
        self.require(SConfig.packages.MPI)
        self.base_patterns = ['StGermain*']
        self.headers = [[os.path.join('StGermain', 'StGermain.h')]]
        self.libraries = [['StGermain']]
        self.symbols = [(['StGermain_Init', 'StGermain_Finalise'], '')]
        self.symbol_setup = 'MPI_Init(&argc, &argv);'
        self.symbol_teardown = 'MPI_Finalize();'
        self.symbol_calls = ['%s(&argc, &argv);', '%s();']
        self.require_shared = True
        self.use_rpath = True

    def get_run_error_message(self, console):
        if len(console):
            return 'Incompatible libraries, check \'config.log\'.'
