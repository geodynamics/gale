import platform
import SCons.Script
import SConfig

class Platform(SConfig.Node):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Node.__init__(self, scons_env, scons_opts, required)

        # Will be set after successful configuration.
        self.system = ''
        self.bits = 0

        # We need to do these now.
        self.check_system()
        self.check_bits()

    def setup_options(self):
        self.opts.AddOptions(
            SCons.Script.BoolOption('with_32bit', 'Generate 32bit code', 0),
            SCons.Script.BoolOption('with_64bit', 'Generate 64bit code', 0)
            )

    def check_system(self):
        self.system = platform.system()
        if not self.system or self.system in ['Linux', 'Unix']:
            self.system = '*ix'

    def check_bits(self):
        if (platform.platform().find('x86_64') != -1 or \
                platform.platform().find('ppc64') != -1 or \
                platform.architecture()[0].find('64') != -1 or \
                self.env['with_64bit']) and \
                not self.env['with_32bit']:
            self.bits = 64
        else:
            self.bits = 32
