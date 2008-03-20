import os
import SConfig

class X11(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Package.__init__(self, scons_env, scons_opts, required)
        self.header_sub_dir = 'X11'
        self.headers = [['Xlib.h']]
        self.libraries = [['X11', 'Xmu']]
        self.symbols = [(['XOpenDisplay'], '')]
        self.symbol_setup = 'void* display;'
        self.symbol_calls = ['display = %s(NULL);']
