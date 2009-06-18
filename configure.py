#!/usr/bin/env python
import sys, subprocess, shutil, os
shutil.copy('test-SConstruct', 'SConstruct')
shutil.copy('config-SConstruct', 'config/config2/SConstruct')
subp = subprocess.Popen(
    'scons/scons.py --config=force -C config/config2 ' + ' '.join(sys.argv[1:]), shell=True
)
subp.wait()
