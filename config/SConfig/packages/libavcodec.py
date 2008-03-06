import os
import SConfig

class libavcodec(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.header_sub_dir = 'ffmpeg'
        self.headers = [['avcodec.h']]
        self.libraries = [['avcodec']]
        self.have_define = 'HAVE_AVCODEC'
