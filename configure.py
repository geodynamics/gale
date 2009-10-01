#!/usr/bin/env python
import sys, subprocess, shutil, os
subp = subprocess.Popen(
    'config/scons/scons.py --config=force -f SConfigure ' + ' '.join(sys.argv[1:]), shell=True
)
subp.wait()
