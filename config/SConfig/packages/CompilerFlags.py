import os, platform
import SCons.Script
import SConfig

class CompilerFlags(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.checks = [self.check_bit_flags,
                       self.check_architecture]

    def setup_options(self):
        SConfig.Package.setup_options(self)
        if not self.opts:
            return
        self.command_options += ['with_32bit']
        self.opts.AddOptions(
            SCons.Script.BoolOption('with_32bit', 'Generate 32bit code', 0),
            SCons.Script.BoolOption('with_64bit', 'Generate 64bit code', 0))

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
        return [1, '', '']

    def check_bit_flags(self):
        if self.try_flag('-m32')[0]:
            self.flag_32bit = '-m32'
        else:
            self.flag_32bit = ''
        if self.try_flag('-m64')[0]:
            self.flag_64bit = '-m64'
        else:
            self.flag_64bit = ''
        return [1, '', '']

    def try_flag(self, flag):
        state = self.env.ParseFlags(flag)
        old = self.push_state(state)
        result = self.run_scons_cmd(self.ctx.TryCompile, '', '.c')
        self.pop_state(old)
        return result
