#!/usr/bin/env python
import os, sys, subprocess
subp = subprocess.Popen(
    'config/scons/scons.py ' + ' '.join(sys.argv[1:]), shell=True
)
subp.wait()

# A manual way of displaying the result of the integration and 
# convergence tests from ./scons.py check*.
if len( sys.argv ) > 1:
    if sys.argv[1] >= "check":
        filename = "summary.dat"
        if os.path.exists( filename ):
            FILE = open( filename, "r" )
            print "--------------------------------------------------------\n" + \
                  "[SYS] SYSTEM TESTS SUMMARY \n" + \
                  "--------------------------------------------------------"
            print FILE.read()
            FILE.close()
            os.remove( filename )
