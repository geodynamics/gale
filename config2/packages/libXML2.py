import os
from config2 import Package

class libXML2(Package):

    def gen_locations(self):
        yield ('', ['/usr/include/libxml2'], [])
        yield ('', ['/usr/local/include'], ['/usr/local/lib'])
        yield ('', ['/usr/local/include/libxml2'], ['/usr/local/lib'])

    def check(self, conf):
        return conf.CheckLibWithHeader('libxml2', os.path.join('libxml', 'parser.h'), 'c')
