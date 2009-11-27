#!/usr/bin/env python

from mercurial import hg, ui, util
import urllib2
import ConfigParser
import os, errno
import sys

# hacky way to find out what the name of the branch is
branch = os.popen('hg branch').readlines()
branch = branch[0].replace("\n","")
cwd = os.getcwd()

if len(sys.argv) > 1:
	deps = [ \
		['https://' + (sys.argv[1] + '@') + 'csd.vpac.org/hg/hgforest', '.hg/forest' ], \
		['https://' + (sys.argv[1] + '@') + 'csd.vpac.org/hg/SConfigure', 'config' ], \
		['https://' + (sys.argv[1] + '@') + 'www.mcc.monash.edu.au/hg/gLucifer', 'gLucifer' ], \
		['https://' + (sys.argv[1] + '@') + 'csd.vpac.org/hg/PICellerator', 'PICellerator' ], \
		['https://' + (sys.argv[1] + '@') + 'csd.vpac.org/hg/StgDomain', 'StgDomain' ], \
		['https://' + (sys.argv[1] + '@') + 'csd.vpac.org/hg/StGermain', 'StGermain' ], \
		['https://' + (sys.argv[1] + '@') + 'csd.vpac.org/hg/StgFEM', 'StgFEM' ], \
		['https://' + (sys.argv[1] + '@') + 'www.mcc.monash.edu.au/hg/Underworld', 'Underworld' ], \
		['https://' + (sys.argv[1] + '@') + 'www.mcc.monash.edu.au/hg/Experimental', 'Experimental' ], \
		['https://' + (sys.argv[1] + '@') + 'www.mcc.monash.edu.au/hg/Experimental', 'Experimental/PDERework/config' ], \
		['https://' + (sys.argv[1] + '@') + 'www.mcc.monash.edu.au/hg/Experimental', 'Experimental/Magma/config' ], \
		['https://' + (sys.argv[1] + '@') + 'www.mcc.monash.edu.au/hg/Experimental', 'Experimental/Geothermal/config' ] ]
else:
	deps = [ \
		['https://csd.vpac.org/hg/hgforest', '.hg/forest' ], \
		['https://csd.vpac.org/hg/SConfigure', 'config' ], \
		['https://www.mcc.monash.edu.au/hg/gLucifer', 'gLucifer' ], \
		['https://csd.vpac.org/hg/PICellerator', 'PICellerator' ], \
		['https://csd.vpac.org/hg/StgDomain', 'StgDomain' ], \
		['https://csd.vpac.org/hg/StGermain', 'StGermain' ], \
		['https://csd.vpac.org/hg/StgFEM', 'StgFEM' ], \
		['https://www.mcc.monash.edu.au/hg/Underworld', 'Underworld' ], \
		['https://www.mcc.monash.edu.au/hg/Experimental', 'Experimental' ] ]
	

# Make sure the '.hg' directory exists
try:
	os.mkdir('.hg')
except OSError, e:
	if( e.errno != errno.EEXIST ):
		raise OSError( e )  # if the error is the directory exists already, keep going


# Download the dependancies...
u = ui.ui()
for dep in deps:
	os.chdir(cwd)
	try:
		print dep[0], '-->', dep[1], '...'
		hg.clone( u, dep[0], dep[1] );
	except util.Abort, e:
		c = ConfigParser.ConfigParser()
		try:
			c.readfp( open( dep[1] + '/.hg/hgrc', 'r' ) )
			if( c.get( 'paths', 'default' ) != dep[0] ):
				print 'Creation failed - ', e, ' but points to another repository'
			else:
				print dep[1], 'already present'
		except:
			print 'Creation failed - ', e, ' and does not seem to be a valid repository'
	except urllib2.URLError, e:
		print 'Download failed - ', e

	# and heres some more hackyness to get hg to switch branches
        # I couldn't find any functionality to do this in the python bindings.
        # JS - 12/11/2008
        os.chdir(dep[1])
        try:
                os.system("hg up -C "+branch)
        except:
                print "fail"

# Tell this root repository that its has forest
os.chdir(cwd)
c=ConfigParser.ConfigParser()
c.read('.hg/hgrc')
try:
	c.add_section('extensions')
except ConfigParser.DuplicateSectionError:
	pass    # If the error is that it exists already, keep going
c.set( 'extensions', 'hgext.forest', os.getcwd() + '/.hg/forest/forest.py' )
c.write( open( '.hg/hgrc', 'w' ) )

