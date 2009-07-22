#!/usr/bin/python
import os, platform
import getopt
import glob
import os.path
import shutil
import string
import sys
import token
import tokenize
import string
from createDocument import Meta

## #######################
# \file This file is designed to act as a doxygen filter to
# change Underworld headers to a C++ format to be parsed
# correctly  into the doxygen html pages, preserving
# inheritance and class
#########################

## Base input class
class Input():
    def __init__(self, filename):
		self.filename = str(filename)
		self.file=""
		self.output=""
		self.functionName = ""
		self.parentName = ""
		self.variablesString = ""
		self.lines = []
		self.addedText=""
    ## Add any meta file data to comments at top of function  
    def addMetaFile(self):
        # check for meta file
        
        filename = os.path.realpath(self.filename)
        nameList = string.rsplit(str(filename), ".")        
        metaFilename = nameList[0] + ".meta"
        
        #open meta file
        if os.path.isfile(metaFilename):
            #parse out meta information using Steve's scripts
            fileData =open(metaFilename)
            metaText = fileData.read()

            metaFile = Meta(metaText)
            metaFile.createMetaDictionaryDtd()
            #print metaFile.dictionary

            # put into C friendly format at top of function
            # allows for dictionaries embedded in lists in the dictionary
            metaText = "/** "
            description = (metaFile.dictionary).pop('Description')
            latexNum = description.count("$")
            # This forces it to 'verbatim" which won't display correctly in latex version
            # TODO make this section 'latex friendly'
            for i in range(latexNum/2):
                description =(description.replace("$", " \\verbatim ", 1)).replace("$", " \\endverbatim ", 1)
            metaText += "\\brief "+description+".  <ul> "
            for key, val in (metaFile.dictionary).iteritems():
                if (isinstance(val, list) == True):
                    if (len(val) <= 0):
                        metaText +=  ""
                    else: 
                        metaText +=  "<li> <b> "+ key +" </b> : \\n <table> "
                        keyList = []
                        for item in val:
                            if (isinstance(item, dict) == True):    
                                keyListTemp = item.keys()
                                if len(keyListTemp) > len(keyList) :
                                    keyList = keyListTemp
                        if len(keyList) > 0 :
                            metaText+= " <tr> \n "
                            for key2 in keyList:
                                metaText += "<th> "+ key2 + "</th>"
                            metaText += "  </tr> \n "
                        for item in val:
                            if (isinstance(item, dict) == True) :
                                metaText += "<tr>"
                                for key2 in keyList :
                                    val2 = item.get(key2)
                                    newVal2 =(val2.replace("$", " \\verbatim ", 1)).replace("$", " \\endverbatim ")
                                    metaText+= " <td> "+ newVal2 +" </td> "                            
                            else:
                                metaText+= "<td> " +str(item) + " </td> "
                            metaText += " </tr>"
                        metaText += " </table> \n </li>"
                                         
                elif (isinstance(val, dict) == True):
                    metaText += "<ul> \n "
                    for key2, val2 in val.iteritems():
                        metaText+= " <li> " +key2 + ": "+ val2 +" </li> "
                    metaText += "</ul> \n "     
                elif (key == "Example") or (key == "Equation") :
                    metaText+= "<li><b> "+key + "</b>: "  +" \\verbatim "+ val +" \\endverbatim " +" </li> \n "
                else:
                    metaText+= "<li><b> "+key + "</b>: " + str(val) +" </li> \n "
            metaText += "</ul> */\n "
            # finding the initial location of the struct for this file
            if self.parentName != "":
			inheritanceText = " : "+self.parentName
            else:
                        inheritanceText = ""
	    searchText = "struct "+ self.functionName + inheritanceText +"\n {\n "
            stringNum = (self.output).find(searchText)

            # Now to create new string
            outText = ""
            stringEnd = len(self.output)
            # add in text before stringNum
            if (stringNum >=0):
                for i in range(stringNum-1):
                    outText +=self.output[i]
                # add in meta text
                outText += metaText.encode('ascii','xmlcharrefreplace')
                # add in rest of file
                for i in range(stringNum, stringEnd-1):
                    outText += self.output[i]
                #print outText
                # save to self.output
                self.output = outText
            else:
                outText =  metaText.encode('ascii','xmlcharrefreplace') + self.output

    ## Convert any "TODO's" into "\todos"
    def convertTodos(self, text):
        #search text for Todos
        self.output = string.replace(text, "TODO", "\\todo")

    ## This function converts the header information into a C++ format.
    def convertHeaders(self,  text) :
        #print text
        stringStart = -1
        stringEnd = -1
        lineCount = 0
        for line in text:
			# find the function Name of the macro:
			# Find line that starts with: "#define __" and ends with " \"
			stringStart = (str(line)).find("#define __" )
			if stringStart >= 0:
				stringEnd =(str(line)).find("\\")
				if stringEnd >=0:
					#print line[stringStart: stringStart+10]
					self.lines = [lineCount]
					self.functionName = ((line[stringStart+10: stringEnd]).lstrip()).rstrip()
					stringNextStart = -1
					stringNextEnd = -1
					stringTemp = ""
					# find the parentName if one exists:
					# Find line that starts with "__" and ends with "\\"
					for i in range(1,15):
					    if (lineCount + i) < (len(text)-1):
					        stringNextStart = (str(text[lineCount +i])).find("__")
					    else:
					        stringNextStart = (str(text[len(text)-1])).find("__")
                                                
					    if stringNextStart >=0:
                                                #print "FOUND __ for parent"
                                                if (lineCount + i) < (len(text)-1): 
						    stringNextEnd =(str(text[lineCount+i])).find("\\")
                                                else:
                                                    stringNextEnd = (str(text[len(text)-1])).find("\\")
						if stringNextEnd >= 0:
                                                    if (lineCount + i) < (len(text)-1):
							self.parentName = ((text[lineCount +i][stringNextStart + 2:stringNextEnd]).lstrip()).rstrip()
                                                    else:
                                                        self.parentName = ((text[len(text)-1][stringNextStart + 2:stringNextEnd]).lstrip()).rstrip()
							#print self.parentName
					    stringNextStart = -1
					    stringNextEnd = -1
			stringStart = -1
			stringEnd = -1
			
			#Find line that creates the inherited struct:
			# Find line that looks like "	struct ",self.functionName," { __",self.FunctionName," };"
			myString = ""
			myString = "struct "
			myNextString = "__"+self.functionName
			stringStart = (str(line)).find(myString )
			if stringStart >=0:
				stringEnd = (str(line)).find(myNextString)
				if stringEnd >=0:
					#print line
					#print lineCount
					(self.lines).append(lineCount)
			
			
			lineCount = lineCount + 1		
	#print self.functionName
	#print self.parentName		
	#print self.lines
		
		# Now, if there is text to replace:
        if self.functionName != "":
			# replace all text between lines [a, b] with new struct structure:
			inheritanceText = ""
			# if no parent function, add all inbetween text
			if self.parentName != "":
				inheritanceText = " : "+self.parentName
			self.addedText = "struct "+ self.functionName + inheritanceText +"\n {\n "
			lineNum = -1
			#print self.lines
			if len(self.lines) == 2:
				for lineNum in range(self.lines[0]+1, self.lines[1]-1):
					if self.parentName == "":
						self.addedText += text[lineNum].split("\\")[0] + "\n "
					else:
						# add inbetween text removing old parent function line
						myText = 0
						myText = text[lineNum].count(self.parentName)
						#print myText
						if myText ==  0:
							self.addedText += text[lineNum].split("\\")[0] + "\n "
				
				#print self.addedText
		
				# replace old text section with new 'addedText'
				for lines in range(0, self.lines[0]-1):
					self.output += text[lines]
				self.output += self.addedText

                                #Add in public: statement just in case.
                                self.output += "public: \n" 
				#print lines
				#print len(text)

                # Find #endif statement and insert bracket before it.
				for lines in range(self.lines[1]+1, len(text)):
					if ((text[lines]).find("#endif") >= 0):
                                           #print "FOUND IT"
					   self.output +=  "};\n" + text[lines]
					else:
					    self.output += text[lines]
				#print "*********"

                
		# if no alterations are needed, then write output to screen as is	
        elif self.functionName == "":
			for lines in range(0, len(text)):
				self.output += text[lines] 
        self.output+="  "
		
    ## main function to tie all class functions together in a run
    def main(self):  	
		# Read in files
		#print self.filename
		self.file = open(self.filename, "r")
		
		text = []
		text = (self.file).readlines()
		#print text[1]
		self.convertHeaders(text)
                self.convertTodos(self.output)       
                self.addMetaFile()

		if self.output != "":
			print "/* *****************************************************************/"
			print "/* THIS FILE HAS BEEN ALTERED TO LOOK LIKE C++ FOR DOXYGEN PARSING */"
			print "/* IT IS NOT THE REAL SOURCE CODE. REFER TO CHECKOUTS OF CODE FOR  */"
			print "/* ACCURATE C CODE DETAILS.                                        */"
			print "/* *****************************************************************/"	
			print self.output
		else:
			print self.file
		

    
## code to run as a standalone from terminal
if __name__ == "__main__":
	# script to run function:
	# sys.argv[1] is the filename to parse.
	values = sys.argv
        if len(values)>1:
	    inputValues = Input(sys.argv[1])

	    inputValues.main()
        else:
            print "NO INPUT FILE!!!!!!!!!!!"            

