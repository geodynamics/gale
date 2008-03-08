import os, platform
import SCons.Script
import SConfig

class CompilerFlags(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.checks = [self.check_architecture]

    def setup_options(self):
        SConfig.Package.setup_options(self)
        if not self.opts:
            return
        self.command_options += ['with_32bit']
        self.opts.AddOptions(
            SCons.Script.BoolOption('with_32bit', 'Generate 32bit code', 0))

    def check_architecture(self):
        if (platform.platform().find('x86_64') != -1 or \
                platform.platorm().find('ppc64') != -1 or \
                platform.architecture()[0].find('64') != -1) and \
                not self.env['with_32bit']:
            self.bits = 64
            #self.env.MergeFlags('-m64') Fails on apac.
        else:
            self.bits = 32
            self.env.MergeFlags('-m32')
            if self.env.subst('$CC') == self.env.subst('$LINK'):
                self.env.AppendUnique(LINKFLAGS='-m32')
        return [1, '', '']
