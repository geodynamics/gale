# vtk2csv
# 
# Reads in Gale pvts and/or pvtu data and combines it into a single
# csv file that can be read by vts2matlab Requires pvpython (included
# with paraview) to be installed.  Invoke it with something like
#
#   pvpython vtk2csv.py output/*.pvt[su]
#
# Tested with paraview 3.6.1
#
# It creates a csv file corresponding to each pvts and pvtu file.
#
# Written by Walter Landry with contributions from Bill Broadley and
# Mark Fleharty
#
# Copyright (C) 2009 California Institute of Technology, University of
# New Mexico, and the Regents of the University of California
#
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

def decode(vtk_file):
    import xml.etree.ElementTree as ElementTree
    root = ElementTree.parse(vtk_file)
    iter = root.getiterator()

    outfile=open(os.path.splitext(vtk_file)[0]+".csv","w")

    fields={}
    components={}
    extents=""
    for element in iter:
        # Get the extents for use by external scripts
        if(element.tag=="Piece"):
            for name, value in element.items():
                if(name=="Extent" or name=="NumberOfPoints"):
                    extents=value
        # Get the data for each individual variable
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
                if element.text:
                    text = element.text
                    text_list = text.split()
                    aDict=dict()
                    aDict.clear()	
                    for i in range (0,len(text_list)):
                        fields[field].append(text_list[i])

    # Print everything out
    if extents:
        outfile.write("# %s\n" % extents)
    outfile.write("# x, y, z, ")
    for j in components:
        if j!="Points":
            if components==1:
                outfile.write("%s, " % j)
            else:
                for k in range(0,components[j]):
                    outfile.write("%s%d, " % (j,k))
    outfile.write("\n")
    for i in range(0,len(fields["Points"])/3):
        for n in range(0,3):
            outfile.write("%s, " % fields["Points"][i*3+n])
        for j in fields:
            if j!="Points":
                for n in range(0,components[j]):
                    outfile.write("%s, " % fields[j][i*components[j]+n])
        outfile.write("\n")
    # Remove the old vtk file
    os.remove(vtk_file)


# First open the parallel vtk xml file and write a serial xml file

from paraview.servermanager import *
Connect()

for i in range(1,len(sys.argv)):
   filename=sys.argv[i]
   print filename

   # Structured Grid
   if os.path.splitext(filename)[1]=='.pvts':
       vts_name=os.path.splitext(filename)[0] + ".vts"
       reader = sources.XMLPStructuredGridReader(FileName=filename)
       writer=writers.XMLStructuredGridWriter(Input=reader,
                                              DataMode=0,
                                              FileName=vts_name)
       writer.UpdatePipeline()
       decode(vts_name)
       
   # Unstructured Grid (particles)
   elif os.path.splitext(filename)[1]=='.pvtu':
       vtu_name=os.path.splitext(filename)[0] + ".vtu"
       reader = sources.XMLPUnstructuredGridReader(FileName=filename)
       writer=writers.XMLUnstructuredGridWriter(Input=reader,
                                                DataMode=0,
                                                FileName=vtu_name)
       writer.UpdatePipeline()
       decode(vtu_name)
   else:
       print "Skipping non-parallel VTK file:",filename


