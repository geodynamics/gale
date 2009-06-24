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

#########################
# This file is designed to act as a doxygen filter to
# change Underworld headers to a C++ format to be parsed
# correctly  into the doxygen html pages, preserving
# inheritance and class
#########################
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
#New function
    def convertMacros(self,  text) :
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
						stringNextEnd =(str(text[lineCount+i])).find("\\")
						if stringNextEnd >= 0:
							self.parentName = ((text[lineCount +i][stringNextStart + 2:stringNextEnd]).lstrip()).rstrip()
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
		

    def main(self):  	
		# Read in files
		#print self.filename
		self.file = open(self.filename, "r")
		
		text = []
		text = (self.file).readlines()
		#print text[1]
		self.convertMacros(text)
       

		if self.output != "":
			print "/* *****************************************************************/"
			print "/* THIS FILE HAS BEEN ALTERED TO LOOK LIKE C++ FOR DOXYGEN PARSING */"
			print "/* IT IS NOT THE REAL SOURCE CODE. REFER TO CHECKOUTS OF CODE FOR  */"
			print "/* ACCURATE C CODE DETAILS.                                        */"
			print "/* *****************************************************************/"	
			print self.output
		else:
			print self.file
		

    
#####
if __name__ == "__main__":
	# script to run function:
	# sys.argv[1] is the filename to parse.
	values = sys.argv
        if len(values)>1:
	    inputValues = Input(sys.argv[1])

	    inputValues.main()
        else:
            print "NO INPUT FILE!!!!!!!!!!!"            

