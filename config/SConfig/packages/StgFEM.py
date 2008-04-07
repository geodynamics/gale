import os
import SConfig

class StgFEM(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Package.__init__(self, scons_env, scons_opts, required)
        self.dependency(SConfig.packages.StgDomain)
        petsc = self.dependency(SConfig.packages.PETSc)
        petsc.have_define = 'HAVE_PETSC'
        self.base_patterns = ['StgFEM*']
        self.headers = [[os.path.join('StgFEM', 'StgFEM.h')]]
        self.libraries = [['StgFEM']]
        self.symbols = [(['StgFEM_Init', 'StgFEM_Finalise'], '')]
        self.symbol_setup = '''MPI_Init(&argc, &argv);
StGermain_Init(&argc, &argv);
StgDomain_Init(&argc, &argv);'''
        self.symbol_teardown = '''StgDomain_Finalise();
StGermain_Finalise();
MPI_Finalize();'''
        self.symbol_calls = ['%s(&argc, &argv);', '%s();']
