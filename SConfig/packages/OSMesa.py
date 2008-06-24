import os
import SConfig

class OSMesa(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False, **kw):
        SConfig.Package.__init__(self, scons_env, scons_opts, required, **kw)
        self.dependency(SConfig.packages.OpenGL)
        self.header_sub_dir = ['GL']
        self.headers = [['osmesa.h']]
        self.libraries = [['OSMesa']]
        self.symbols = [(['OSMesaCreateContext', 'OSMesaDestroyContext'], '')]
        self.symbol_setup = 'void* ctx;'
        self.symbol_calls = ['ctx = %s(OSMESA_RGBA, NULL);', '%s(ctx);']
