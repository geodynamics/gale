import os
import SConfig

class pcu(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Package.__init__(self, scons_env, scons_opts, required)
        self.dependency(SConfig.packages.MPI)
        self.headers = [[os.path.join('pcu', 'pcu.h')]]
        self.libraries = [['pcu']]
        self.checks += [self.check_scons_script]

        self.scons_script = ''

    def check_scons_script(self):
        if self.base_dir:
            script = os.path.join(self.base_dir, 'script', 'pcu', 'scons.py')
            if os.path.exists(script):
                self.scons_script = script
                self.ctx.Display('      Found SCons builder script.\n')
        return True

    def enable(self, scons_env, old_state=None):
        SConfig.Package.enable(self, scons_env, old_state)
        if self.scons_script:
            self.backup_variable(scons_env, 'CONFIGSCRIPTS', old_state)
            scons_env.AppendUnique(CONFIGVARS=['CONFIGSCRIPTS'])
            scons_env.AppendUnique(CONFIGSCRIPTS=[self.scons_script])
            env = scons_env
            scons_env.SConscript(self.scons_script, 'env')
