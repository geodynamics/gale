import os
import SConfig

class X11(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.header_sub_dir = 'X11'
        self.headers = [['Xlib.h']]
        self.libraries = [['X11', 'Xmu']]
        self.have_define = 'HAVE_X11'
        self.symbols = [(['XOpenDisplay'], '')]
        self.symbol_setup = 'void* display;'
        self.symbol_calls = ['display = %s(NULL);']
