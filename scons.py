#!/usr/bin/env python
import sys, subprocess
subp = subprocess.Popen(
    'config/scons/scons.py ' + ' '.join(sys.argv[1:]), shell=True
)
subp.wait()
