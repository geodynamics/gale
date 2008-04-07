import os
import SConfig

class StGermain(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Package.__init__(self, scons_env, scons_opts, required)
        self.dependency(SConfig.packages.cmath)
        self.dependency(SConfig.packages.libXML2)
        self.dependency(SConfig.packages.MPI)
        self.dependency(SConfig.packages.SVNRevision)
        self.base_patterns = ['StGermain*']
        self.headers = [[os.path.join('StGermain', 'StGermain.h')]]
        self.libraries = [['StGermain']]
        self.symbols = [(['StGermain_Init', 'StGermain_Finalise'], '')]
        self.symbol_setup = 'MPI_Init(&argc, &argv);'
        self.symbol_teardown = 'MPI_Finalize();'
        self.symbol_calls = ['%s(&argc, &argv);', '%s();']

    def enable(self, scons_env, old_state=None):
        SConfig.Package.enable(self, scons_env, old_state)
        if self.base_dir:
            script = os.path.join(self.base_dir, 'script', 'pcu', 'scons.py')
            if os.path.exists(script):
                self.backup_variable(scons_env, 'CONFIGSCRIPTS', old_state)
                scons_env.AppendUnique(CONFIGVARS=['CONFIGSCRIPTS'])
                scons_env.AppendUnique(CONFIGSCRIPTS=[script])
                env = scons_env
                scons_env.SConscript(script, 'env')

            script = os.path.join(self.base_dir, 'script', 'StGermain', 'scons.py')
            if os.path.exists(script):
                self.backup_variable(scons_env, 'CONFIGSCRIPTS', old_state)
                scons_env.AppendUnique(CONFIGVARS=['CONFIGSCRIPTS'])
                scons_env.AppendUnique(CONFIGSCRIPTS=[script])
                env = scons_env
                scons_env.SConscript(script, 'env')
