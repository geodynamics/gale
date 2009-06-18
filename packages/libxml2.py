import os
from config import Package

class libXML2(Package):

    def gen_locations(self):
        yield ('', ['/usr/include/libxml2'], [])
        yield ('/usr/local', ['/usr/local/include'], ['/usr/local/lib'])
        yield ('/usr/local', ['/usr/local/include/libxml2'], ['/usr/local/lib'])

    def gen_envs(self, loc):
        for env in Package.gen_envs(self, loc):
            env['pkg_headers'] = [os.path.join('libxml', 'parser.h')]
            if self.find_libraries(loc[2], 'libxml'):
                env.PrependUnique(LIBS=['libxml2'])
                yield env
