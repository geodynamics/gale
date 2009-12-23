
import sys
import xml.etree.ElementTree as ElementTree
print (sys.argv[1])
root = ElementTree.parse(open(sys.argv[1]))
iter = root.getiterator()

x=[]
y=[]
z=[]
visl=[]
denl=[]
indl=[]
for element in iter:
	xyz=0
	vis=0
	den=0
	ind=0
	print ("Element:", element.tag)
	if element.keys():
		print ("\tAttributes:")
		for name, value in element.items():
			if ((name == 'NumberOfComponents') and ( value == '3' )):
				xyz=1
			if ((name == 'Name') and ( value == 'Viscosity')):
				vis=1
			if ((name == 'Name') and ( value == 'Density')):
				den=1
			if ((name == 'Name') and ( value == 'Material_Index')):
				ind=1
			print ("\t\tName: '%s', Value: '%s'"%(name, value))
	print ("\tChildren:")
	if element.text:
		text = element.text
#        text = len(text) > 40 and text[:40] + "..." or text
		list = text.split()
		print ("Length of list = ",len(list))
		aDict=dict()
		aDict.clear()	
#		minlen=min(len(list),100)
#		print "list=",list
		if (len(list)>0):
			for i in range (0,len(list)):
				if (xyz):
					index=0
					if (index>0):
					 	index=int(index/3)
#					print i,xyz,"x,y,z=",list[i*3],list[i*3+1],list[i*3+2]," index=",index
					if (i%3 == 0):
						x.append(list[i])
					if (i%3 == 1):
						y.append(list[i])
					if (i%3 == 2):
						z.append(list[i])
				elif (vis):
					visl.append(list[i])
#					print i,"vis=",list[i]
				elif (den):
					denl.append(list[i])
#					print i,"den=",list[i]
				elif (ind):
					indl.append(list[i])
#					print i,"ind=",list[i]
#			else:
#				print i,vis,list[i]
#            if (aDict.get(list[i])):
#                aDict[list[i]]=aDict[list[i]]+1
#            else:
#                aDict[list[i]]=1	
#        if (len(aDict)<100):
#        	print (aDict)
#	for key in aDict:
#		print "this key =",key," qty=",aDict(key)

#        print "\t\tText:", repr(text)
	if element.getchildren():
		for child in element:
			print ("\t\tElement", child.tag)
			if child.tail:
				text = child.tail
#                text = len(text) > 40 and text[:40] + "..." or text
				print ("\t\tText: ", repr(text))  

print "len x =", len(x)
#print "x =", x
for i in range(0,len(x)):
	print x[i],y[i],z[i],visl[i],denl[i],indl[i]
