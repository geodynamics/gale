import os
import SConfig

class libavcodec(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False, **kw):
        SConfig.Package.__init__(self, scons_env, scons_opts, required, **kw)
        self.header_sub_dir = 'ffmpeg'
        self.headers = [['avcodec.h']]
        self.libraries = [['avcodec']]
        self.have_define = 'HAVE_AVCODEC'
