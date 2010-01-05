#!/usr/bin/env python 
"""
This scipt will update Underworld repositories to a specific branch.
The one argument this script uses is the branch name, eg.
	./updateRepos.py v1.3.x

"""

from mercurial import hg, ui, util
import urllib2
import ConfigParser
import os, errno
import sys

deps = [ 
	 'gLucifer', 'PICellerator', 'StgDomain' , 'StGermain' , 'StgFEM' , 'Underworld' , 'Experimental' \
		] #, 'Experimental/PDERework/config', 'Experimental/Magma/config' ]

cwd = os.getcwd()

# check if there are 
if len(sys.argv) != 2:
   print "ERROR - must supply one argument only (the branch name), currently\n" 
   print sys.argv
   print "\nExample usage: ./updateRepos.py v1.3.x\n\n"
   sys.exit()

for dep in deps:
   if not os.path.exists(cwd + "/" + dep ):
      continue

   os.chdir(cwd + "/" + dep)
   # get the branch name
   branch = os.popen('hg branch').readlines()
   branch = branch[0].replace("\n","")
   # check branch
   os.system("hg up -C " + sys.argv[1])
   print "updating " + dep + " to branch " + sys.argv[1]

