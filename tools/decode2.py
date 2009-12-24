#!/usr/bin/env python

# Copyright (C) 2009 Bill Broadley
# Modified by Walter Landry

#  This software is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation; either version 2 of the
#  License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this library; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
#  02110-1301 USA



import sys
import xml.etree.ElementTree as ElementTree
root = ElementTree.parse(open(sys.argv[1]))
iter = root.getiterator()

if (len(sys.argv)>2):
        sys.stdout=open(sys.argv[2],"w")

x=[]
y=[]
z=[]
fields={}
components={}
for element in iter:
        if(element.tag=="DataArray"):
                field=""
                comps=0
		for name, value in element.items():
			if (name == 'NumberOfComponents'):
                                comps=int(value)
                        elif (name=='Name' and value!="offsets" and value!="types"):
                                fields[value]=[]
                                field=value
                if field:
                        if comps:
                                components[field]=comps
                        else:
                                components[field]=1
                if element.text and (field or comps):
                        text = element.text
                        text_list = text.split()
                        aDict=dict()
                        aDict.clear()	
                        for i in range (0,len(text_list)):
                                if (not field):
                                        if (i%3 == 0):
                                                x.append(text_list[i])
                                        if (i%3 == 1):
                                                y.append(text_list[i])
                                        if (i%3 == 2):
                                                z.append(text_list[i])
                                elif (field):
                                        fields[field].append(text_list[i])

print "x y z",
for j in components:
        print j,components[j],
print
for i in range(0,len(x)):
	print x[i],y[i],z[i],
        for j in fields:
                for n in range(0,components[j]):
                        print fields[j][i*components[j]+n],
        print
