import os, platform
import SCons.Script
import SConfig

class CompilerFlags(SConfig.Node):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Node.__init__(self, scons_env, scons_opts, required)
        self.checks = [self.check_bit_flags,
                       self.check_architecture]

    def setup_options(self):
        SConfig.Node.setup_options(self)
        self.opts.AddOptions(
            SCons.Script.BoolOption('with_32bit', 'Generate 32bit code', 0),
            SCons.Script.BoolOption('with_64bit', 'Generate 64bit code', 0)
            )

    def check_architecture(self):
        if (platform.platform().find('x86_64') != -1 or \
            platform.platform().find('ppc64') != -1 or \
            platform.architecture()[0].find('64') != -1 or \
            self.env['with_64bit']) and \
            not self.env['with_32bit']:
            self.bits = 64
            if self.flag_64bit:
                self.env.MergeFlags(self.flag_64bit)
                if self.env.subst('$CC') == self.env.subst('$LINK'):
                    self.env.AppendUnique(LINKFLAGS=[self.flag_64bit])
        else:
            self.bits = 32
            if self.flag_32bit:
                self.env.MergeFlags(self.flag_32bit)
                if self.env.subst('$CC') == self.env.subst('$LINK'):
                    self.env.AppendUnique(LINKFLAGS=[self.flag_32bit])
        return True

    def check_bit_flags(self):
        if self.try_flag('-m32')[0]:
            self.flag_32bit = '-m32'
        elif self.try_flag('-q32')[0]:
            self.flag_32bit = '-q32'
        else:
            self.flag_32bit = ''
        if self.try_flag('-m64')[0]:
            self.flag_64bit = '-m64'
        elif self.try_flag('-q64')[0]:
            self.flag_64bit = '-q64'
        else:
            self.flag_64bit = ''
        return True

    def try_flag(self, flag):
        state = self.env.ParseFlags(flag)
        old = self.push_state(state)
        result = self.run_scons_cmd(self.ctx.TryCompile, '', '.c')
        self.pop_state(old)
        if result[0] and (result[1].find('not recognized') != -1 or
                          result[1].find('not recognised') != -1 or
                          result[1].find('unknown option') != -1):
            result[0] = 0
        return [result[0], result[1], '']
