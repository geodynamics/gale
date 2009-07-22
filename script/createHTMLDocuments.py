#!/usr/bin/python
import os, os.path, sys,  datetime,  string,  math,  subprocess
#Modules created for this project
import convert,  util
import createDocument
from createDocument import createDocuments
from createDocument import Project
from createDocument import Meta
## The createHTMLDocuments class. Inherits from createDocuments. Contains all the information
# needed to create a series of HTML pages documenting the Component Codex. Should be flexible
# enough to adapt to being used for different displays of the meta files.
class createHTMLDocuments(createDocuments):
    def __init__(self, documentName, documentBlurb,  projectList,  documentPath, styleSheetList,  externalWebpages,  contactEmail, scripts,  pictureDirectory,  indexFileName,  searchBarType,  searchPath):
        createDocuments.__init__(self, documentName,  projectList, documentPath)
        #Set universal page values
        self.styleSheetList = styleSheetList
        self.externalWebpages = externalWebpages
        self.contactEmail = contactEmail
        self.scripts = scripts
        self.blurb = documentBlurb
        self.div= []
        self.pictureSubDirectoryName = pictureDirectory
        self.indexFileName = indexFileName
        self.searchBarType = searchBarType
        self.searchPath = searchPath
    ## Create the HTML Pages.    
    def createHTMLPages(self) :
        #check output path
        createDocument.checkOutputPath(self.path)
        #create project pages
        for count in range(len(self.projectList)):
            self.createHTMLPage(count)
            
        #create Index page    
        self.createIndexPage()

    ## Creates the HTML Page for each project.    
    def createHTMLPage(self, count):
        # create div ids list
        self.div = DivIds()
        self.div.createDivIds(self.projectList[count].name)
        projectName = self.projectList[count].name
        print "Creating HTML page for project: " + str(projectName)
        # create string for content of page
        htmlPage =u'<html>\n'
        
        htmlPage += unicode(self.createHeader(projectName))
        htmlPage +=u'<body>\n'
        htmlPage += unicode(self.createTopMenu())
        
        htmlPage += u'<div id="'+unicode(self.div.container)+u'">\n'
        htmlPage += unicode(self.createBreadcrumbTrail(projectName))
        htmlPage += unicode(self.createNavitab(projectName))
        htmlPage += unicode(self.createSidebar(count))
        
        htmlPage += u'<div id="'+unicode(self.div.main)+u'">\n'         
        htmlPage += unicode(self.createBlurb(projectName))       
        htmlPage += unicode(self.createMain(count))
        htmlPage +=u'</div>\n'
        htmlPage += unicode(self.createFooter())
        htmlPage += u'\n</div>'
        htmlPage +=unicode(self.addScripts())
        htmlPage +=u'\n</body>\n</html>'
        htmlPage = unicode(htmlPage)
        # write page to file
        filename = self.path + '/'+ self.projectList[count].name + '.html'
        page = open(filename,  'w')
        page.write(htmlPage)
        page.close()
    ## The Main page.It's contents are mostly specified by an external file.    
    def createIndexPage(self):
        self.div = DivIds()
        self.div.createDivIds('')
        print "Creating index Page"
        projectName='index'
        htmlPage ='<html>\n'
        
        htmlPage += self.createHeader("")
        htmlPage +='<body>\n'
        htmlPage += self.createTopMenu()
        
        htmlPage += '<div id="'+self.div.container+'">\n'
        htmlPage += self.createIndexBreadcrumbTrail()
        htmlPage += self.createNavitab(projectName)
        htmlPage += '<div id="'+self.div.sidebar+'">\n'
        # Set search path for searching through the webpages.
        myPath = self.searchPath
        if myPath == "":
            myPath = self.path
        #print self.searchPath, myPath
        if self.searchBarType == 'Google':
            htmlPage += self.addSearchBarGoogle(myPath)
        elif self.searchBarType == 'None':
            htmlPage += ''
        else:
            print "Search Bar Type: " + self.searchBarType + " not available."
        htmlPage+= '</div>\n'
        htmlPage += '<div id="'+self.div.main+'">\n'         
        htmlPage += self.createIndexBlurb("Main page for "+ self.name)      
       # Extract out content from the main file        
        htmlPage += self.createIndexContent()
        htmlPage +='</div>\n'
        htmlPage += self.createFooter()
        htmlPage += '\n</div>'
        htmlPage +=self.addScripts()
        htmlPage +='\n</body>\n</html>'
        
        # write page to file
        filename = self.path + '/index.html'
        page = open(filename,  'w')
        page.write(htmlPage)
        page.close()
    ## Extract content from file for main page
    def createIndexContent(self):
        text = ''
        if not(self.indexFileName == ''):
            f = open(self.indexFileName,  'r')
            text += f.read()
            f.close()
            
        return text
    ## Add a google search bar for specified path.
    def addSearchBarGoogle(self,  path):
        text = ''
        text += '<!-- SiteSearch Google: Adds a Google search bar for the Codex.\n' 
        text += 'From: http://www.google.com.au/searchcode.html  -->\n'

        text += '<FORM method=GET action="http://www.google.com.au/search">\n'
        text += '<input type=hidden name=ie value=UTF-8>\n'
        text += '<input type=hidden name=oe value=UTF-8>\n'
        text += '<TABLE bgcolor="#FFFFFF"><tr><td>\n'
        text += '<A HREF="http://www.google.com.au/">'
        text += '<IMG SRC="http://www.google.com.au/logos/Logo_40wht.gif"' 
        text += 'border="0" ALT="Google"></A>\n'
        text += '</td></tr>\n'
        text += '<tr><td>\n'
        text += '<INPUT TYPE=text name=q size=31 maxlength=255 value="">\n'
        text += '<INPUT type=submit name=btnG VALUE="Google Search">\n'
        text += '<font size=-1>\n'
        text += '<input type=hidden name=domains'+\
        ' value="'+path+'"><br><input type=radio name=sitesearch '+\
        'value=""> WWW <input type=radio name=sitesearch'+\
        ' value="'+path+'" checked> '+path+' <br>\n'
        text +='</font>\n'
        text +='</td></tr></TABLE>\n'
        text +='</FORM>\n'
        text +='<!-- SiteSearch Google -->\n'


        return text
    ## Header data for webpages.    
    def createHeader(self,  projectName):
        text = '<head>\n<META http-equiv="Content-Type" content="text/html; charset=UTF-8">' 
        text +='<title>' +projectName +' '+ self.name+'</title>\n'
        for stylesheetPath in self.styleSheetList:
            (path, stylesheet) = os.path.split(stylesheetPath) 
            if string.find(stylesheet,  'print.css') > -1:
                text +='<link rel="stylesheet" type="text/css" href="'+ stylesheet + '" media="print" >\n'
            else:
               text +='<link rel="stylesheet" type="text/css" href="'+ stylesheet + '" media="screen" >\n'
            
        text += '</head>\n'
        return text
    ## Top menu to appear along very top of webpages.
    def createTopMenu(self):
        text = '<div id="'+self.div.toptab+'">\n'
        for webpage in self.externalWebpages:
            text +='<a class="'+self.div.toptab+'" href="'+webpage[0]+'">'+webpage[1]+'</a><span class="hide"> | </span>\n'

        text += '</div>\n'
        return text
    ## Specific div from the codex stylesheet
    def createBreadcrumbTrail(self, projectName):
        text ='<div id="'+self.div.breadcrumbtrail+'">\n<h1>\n' 
        text +=        '<a href="'+projectName+'.html'+'">'+ projectName 
        text +=        ' '+self.name+' </a>\n </h1>\n </div>\n'
        return text
        
    ## Specific div from the codex stylesheet for the main page
    def createIndexBreadcrumbTrail(self):
        text ='<div id="'+self.div.breadcrumbtrail+'">\n<h1>\n' 
        text +=        '<a href="index.html'+'">'
        text +=        ''+self.name+' </a>\n </h1>\n </div>\n'
        return text

    ## Create the tabs to navigate between the project pages.
    def createNavitab(self,  projectName):
        text = '<div id="navitab">\n <br>\n '
        if projectName == 'index':
            text +=        '<a class="'+self.div.navitab[1]+'" href="index.html">Main</a><span class="hide"> | </span>'        
        else:
            text +=        '<a class="'+self.div.navitab[0]+'" href="index.html">Main</a><span class="hide"> | </span>'
        for project in self.projectList:
            if project.name ==projectName:
                text  += '<a class="'+self.div.navitab[1]+'" href="'+project.name+'.html">'+project.name+'</a><span class="hide"> | </span>'
            else:
                text  += '<a class="'+self.div.navitab[0]+'" href="'+project.name+'.html">'+project.name+'</a><span class="hide"> | </span>'
        text += '\n<hr>\n </div>\n'
        return text
        
    ## Create the footer with contact info, and page creation time.
    def createFooter(self):
        text = '<div id="'+self.div.footer+'">\n' 
        text +=        '<h3>Last updated: '+str(datetime.datetime.today())+'</h3>\n' 
        text +=        'Email:  <a href="mailto:'+self.contactEmail+ '">'+self.contactEmail+'</a>.\n'  
        text +=        '</div>\n'
        return text

    ## Create the sidebar. This searches though all the project metas and attempts to auto-sort
    # them into a reasonably sized 2-tier menu system. Not perfect.
    def createSidebar(self,  count):
        # Get out name of each meta file in project
        project = self.projectList[count]

        text = '<div id="'+self.div.sidebar+'">\n' 
        text +=       '<a name="'+project.name+'">\n'
        text +=        '<h2>'+project.name+'</h2>\n </a>\n' 
        text +=       '<h3> '+self.name+' Menu </h3>\n'
        text += '<ul class="menu">\n'
        #Menu items
        names = []
        for meta in project.metas:
            if project.metaFlag == 'xsd':
                names.append(meta.dictionary['info']['title'])
            elif project.metaFlag =='dtd':
                names.append(meta.dictionary['Name'])
        names.sort()
        #Maximum size for 1st tier on  2-tier  sorting algorithm. (quite arbitrary):
        sidebarLengthMax = 10
        if len(names) <= sidebarLengthMax:
            #single tier adequate
            for name in names:
                #Add menu item
                text +=	'<li><a class="'+self.div.sidebarTiers[1]+'" href="./'+ project.name+'.html#'+name+'">'+name+'<br/></a></li>\n'
        else:
                #2-tier sort. This is designed to create a 2-Tier menu which is self generated based
                # on the number of meta files, and the meta file names. It will roughly be 5 tier 1 items
                # with up to 8 meta in each. The string creation is to create usable labels for the Tier 1 menus
                # that have unique ranges based on the ordered set of meta names.
                sidebarMaxTier1Menu = 5
                sidebarMaxTier2Menu = 8
                for ind in range(4):
                    perMenu = math.ceil(len(names)/sidebarMaxTier1Menu)
                    if perMenu > sidebarMaxTier2Menu:
                        sidebarMaxTier1Menu = sidebarMaxTier1Menu + 1
                #check menus not too big.
                for ind in range(4):
                    if (sidebarMaxTier1Menu * (sidebarMaxTier2Menu-1)) > len(names):
                        sidebarMaxTier2Menu = sidebarMaxTier2Menu - 1
                nextString =""
                newStartString =""
                startString =""
                endString =""
                for i in range(sidebarMaxTier1Menu):
                    if ((i+1)*(sidebarMaxTier2Menu)-1) >=  len(names):
                        finish = len(names)-1
                    else:
                        finish =((i+1)*(sidebarMaxTier2Menu)-1)
                    start =    i*sidebarMaxTier2Menu
                    # checks for repitition in A to B and in the next menus C to D
                    # eg "meshM to meshPo" followed by "meshPr to O"
                    startString = string.capitalize(names[start][0] + names[start][1])
                    endString = string.capitalize(names[finish][0] + names[finish][1]) 
                    if (finish + 1) < len(names):
                        altEndString =string.capitalize(names[finish][0] + names[finish][1]) 
                        nextString = string.capitalize(names[finish+1][0] + names[finish+1][1]) 
                    else:
                        altEndString =string.capitalize(names[finish][0] + names[finish][1])                         
                        nextString =""
                    for k in range(2, 6):
                        if startString == endString:
                            startString = string.capitalize(startString+ names[start][k])
                            endString = string.capitalize(endString + names[finish][k])
                    if (len(startString) < len(newStartString)):
                        startString = newStartString
                    for k in range(2, 6):
                        if altEndString == nextString:
                            altEndString = string.capitalize(altEndString+ names[start][k])
                            nextString = string.capitalize(nextString + names[finish][k])
                    if (len(endString) < len(altEndString)):
                        endString = altEndString
                    newStartString = nextString
                    text +='      <li> <a class="'+self.div.sidebarTiers[0]+'">'+startString+' to '+endString+'</a>\n'
                    text += '<ul>\n'
                    for j in range(sidebarMaxTier2Menu):
                        index = i*sidebarMaxTier2Menu + j
                        if index < len(names):
                            text +=	'<li><a class="'+self.div.sidebarTiers[1]+'" href="./'+ project.name+'.html#'+names[index]+'">'+names[index]+'<br/></a></li>\n'
                    text += '</ul></li>\n'
                    
        text +='<li><a href="./'+project.name+'.html#top">Back to top</a>\n</li>'
        text += '</ul>\n'                        
        text += '</div>\n'
        return text
    ## Put small blurb in box at top of webpage. Blurb is given from command line.    
    def createBlurb(self,  projectName):
        text =        '<div id="'+self.div.desc+'">\n'
        text +=        '<h3> '+projectName+' '+self.name+' </h3>\n' 
        text +=        '<p>'+self.blurb+'<br>\n </p>\n </div>\n<br>\n'
        return text
        
    ## Put small blurb on index page, different from project pages. This is fixed
    def createIndexBlurb(self,  blurb):
        text =        '<div id="'+self.div.desc+'">\n'
        text +=        '<h3>'+self.name+'</h3>\n' 
        text +=        '<p>'+blurb+'<br>\n </p>\n </div>\n<br>\n'
        return text 
        
    ## Create the meta-file content part of the pages.    
    def createMain(self,  count):
        text = ''
        #create all meta entries for current Project page
        project =self.projectList[count]
        #Sort meta names into an ordered list
        print "Creating Meta data entries for project: " + str(project.name) 

        names = []
        for meta in project.metas:
            if project.metaFlag == 'xsd':
                names.append(meta.dictionary['info']['title'])
            elif project.metaFlag =='dtd':
                names.append(meta.dictionary['Name'])
        names.sort()
        for name in names:
	    
            #Select meta type:
            if project.metaFlag == 'dtd':
                text += self.createMetaEntryDtd(name, count)
            elif project.metaFlag =='xsd':
                text += self.createMetaEntryXsd(name, count)
        return text
    ## For the simpler componentInfo entries, can easily fit into this function.    
    def addSimpleComponentInfo(self, name, value):
                text = ''
                text += '<div id="'+unicode(self.div.componentInfo[0])+'">\n'
                text += '<div id="'+unicode(self.div.componentInfo[1])+'">\n'
                text += '<b>'+unicode(name)+'</b>: </div>\n'
                text += '<div id="'+self.div.componentInfo[2]+'">'+value+'<br>\n'
                text += '</div>\n</div>\n'        
                return unicode(text)
    ## Example part of meta component            
    def addExampleInfo(self,  name,  value):
        text = ''
        text += '<div id="'+self.div.componentInfo[0]+'">\n'
        text += '<div id="'+self.div.componentInfo[1]+'">\n'
        text += '<b>'+str(name)+'</b>:</div>\n'
        text += '<div id="'+self.div.componentInfo[2]+'">\n<div id="codebox">\n'
        text += '<xmp>'+str(value)+'</xmp>\n'
        text += '<br>\n</div>\n</div>\n</div>\n'

        return text
    ## List entry html code    
    def addListComponentInfo(self,  name, tableId,  value,  componentName):
        text = ''
        text += '<div id="'+self.div.componentInfo[0]+'">\n'
        text += '<div id="'+self.div.componentInfo[1]+'">\n'
        text += '<b>'+name+'</b>\n'
        text += '</div>\n<div id="'+string.lower(name)+'">\n'
        text += '<table id="'+str(tableId)+'">\n'
        
        titleList = []
        #sort the inner dictionary keys of the first entry
        for item in value:
            tempList = item.keys()
            for title in tempList:
                # if title isn't in titleList, add it
                if titleList.count(title) == 0:
                    titleList.append(title)
        titleList.sort()
        # If name, type exist, put them first. otherwise, alphabetical sort.
        if titleList.count('Type') > 0 :
            titleList.remove('Type')
            titleList.insert(0,  "Type")
        if titleList.count('type') > 0 :
            titleList.remove('type')
            titleList.insert(0,  "type")
            
        if titleList.count('Name') > 0:
            titleList.remove('Name')
            titleList.insert(0,  "Name")
            
        if titleList.count('name') > 0:
            titleList.remove('name')
            titleList.insert(0,  "name")
            
        # Now add titles
        text += '<tr>\n'
        for title in titleList: 
            text += '<th>'+str(title)+'</th>'    
        text += '\n</tr>\n'
        
        # Now add items associated with titles
        picsIndex = 0
        
        for item in value:
            text += '<tr>\n'
            for title in titleList:
                if item.has_key(title):
                    # check for latex code to convert Does multiple latex entries in one item.
                    if string.count(item[title], '$') >1:
                        #split string at $ signs
                        titleString = " " + item[title]
                        myList = string.split(titleString,  '$')
                        # all 'even' occurances(index 1,3,5 etc as nos start at 0) in this string will be latex code, even if string starts with a $.
                        for listValue in range(len(myList)/2):
                            latexString = '$' + myList[listValue+1]+'$'
                            #call latex function
                            picName = self.addEquationPicturePng(componentName +name,  latexString,  picsIndex)
                            newString = '<img src="'+self.pictureSubDirectoryName+'/'+picName+'" ><br>\n'
                            myList[listValue+1] = newString
                            picsIndex = picsIndex + 1
                         #turn list back into string   
                        text += '<td>'+ (string.lstrip(string.join(myList))).encode('ascii','xmlcharrefreplace')+'</td>'    
                    else:
                        # add code to table
                        text += '<td>'+(item[title]).encode('ascii','xmlcharrefreplace')+'</td>'
                else:
                    text += '<td> </td>'
            text +='</tr>\n'
        text += '\n</table>\n</div>\n</div>'
        return text
        
    ## Call latex code to convert latex to pictures in html   
    def addEquationPicturePng(self,  name,  value,  index):
        # add to pics sub-directory 
        picturePath = os.path.normpath(self.path) + '/'+ self.pictureSubDirectoryName
        
        if os.path.exists(self.path) and not os.path.exists(picturePath) :
            # create pics directory
            print "Creating directory: "+ picturePath
            os.mkdir(picturePath)
        #create picture name
        pictureName = name + str(index) + '.png'
        #Todo: create equation pic from latex source.
        # create latex file.
        latexCode = "\\documentclass{article}\n"
        latexCode += "\\usepackage{amsmath}\n"
        latexCode += "\\usepackage{amsthm}\n"
        latexCode += "\\usepackage{amssymb}\n"
        latexCode += "\\usepackage{bm}\n"
        latexCode += "\\pagestyle{empty}\n" 
        latexCode += "\\begin{document}\n"
        latexCode += value + "\n"
        latexCode += "\\end{document}"
        latexFileBase = picturePath + "/" + "temp"
        latexFileName = latexFileBase +".tex"
        latexFile = open(latexFileName, 'w')
        latexFile.write(latexCode)
        latexFile.close()
        
        # change working directory to the same as source's
        cwd = os.getcwd()
        os.chdir(picturePath)
        # run latex on file
        #os.system('latex  '+ latexFileName)
        try:
            retcode = subprocess.call('latex  '+ latexFileName + ">>temp.txt", shell=True)
            if retcode < 0:
                print >>sys.stderr, "Child was terminated by signal: latex file could not be run for creating " + pictureName, -retcode
            
        except OSError, e:
            print >>sys.stderr, "Execution failed:", e 
        # Run dvipng on the generated DVI file. Use tight bounding box. 
        # Magnification is set to 1200
        cmd = "dvipng -T tight -z 9 -bg Transparent "\
        + "-o " + self.path+ "/" + self.pictureSubDirectoryName + "/" + pictureName + " " + latexFileBase + '.dvi' + " >>temp.txt"
        #os.system(cmd) 
        try:
            retcode = subprocess.call(cmd, shell=True)
            if retcode < 0:
                print >>sys.stderr, "Child was terminated by signal", -retcode
            #else:
                #print >>sys.stderr, "Picture created for ", pictureName
        except OSError, e:
            print >>sys.stderr, "Execution failed:", e 

# remove temporary files
        
        os.remove(latexFileBase+'.tex')
        os.remove(latexFileBase+'.log')
        os.remove(latexFileBase+'.aux')
        os.remove(latexFileBase+'.dvi')
        os.remove('temp.txt')
        os.chdir(cwd)
        
        return pictureName
     
    ## Add equation list value. 
    def addEquationInfo(self, componentName,  value):

        # add proper coding reference
        text = ''
        text += '<div id="'+self.div.componentInfo[0]+'">\n'
        text += '<div id="'+self.div.componentInfo[1]+'">\n'
        text += '<b>Equation</b>: </div>\n'
        text += '<div id="'+self.div.componentInfo[2]+'">\n'
        # create picture
        pictureName = self.addEquationPicturePng(componentName,  value, 1)
        text += '<img src="'+self.pictureSubDirectoryName+'/'+pictureName+'" align="top"><br>\n'
        text += '</div>\n</div>\n'

        return text
        
    ## Function Puts the parent info in, and creates a html link to that entry.    
    def addParentInfoDtd(self, titleName,  nameTitle, parentName):
        myString  = ''
        #search through all projects to find meta with parentName.
        for project in self.projectList:
            for meta in project.metas:
                name = meta.dictionary[nameTitle]
                
                if name == parentName:
                    #If parent found, add as hyperlink
                    myString += '<a href="./'+ project.name+'.html#'+name+'">'+name+'<br/></a>'
        text  = self.addSimpleComponentInfo(titleName,myString)
        return text
        
    ## Find all children of given meta component.    
    def addChildrenInfoDtd(self,  titleName,  parentTitle, nameTitle, parentName):
        myString  = '<ul>\n'
        #search through all metas to find meta that list parentName as their parent.
        for project in self.projectList:
            for meta in project.metas:
                if meta.dictionary.has_key('Parent'):
                    name = meta.dictionary['Parent']
                else:
                    #print str(meta.dictionary[nameTitle]) + " has no "+parentTitle+" information."
                    name ="" 
                if name == parentName:
                    #If child found, add as hyperlink
                    myString += '<li><a href="./'+ project.name+'.html#'+meta.dictionary[nameTitle]+'">'+meta.dictionary[nameTitle]+'<br/></a></li>\n'
        myString +='</ul>'            
        text  = self.addSimpleComponentInfo(titleName,myString)
        return text
        
    ## Create the individual meta entry, classic Dtd style.    
    def createMetaEntryDtd(self,  name,  count):
        text =""
        #create individual meta entry
        project = self.projectList[count]
        #Select Meta

        for meta in project.metas:
            if meta.dictionary['Name'] == name:
                # Add meta data to file
                #print project.name,  name,  meta.dictionary['Location']
                text +='<h3 align="center">\n'
                
                text +='<a name="'+unicode(meta.dictionary['Name'])+'">'+unicode(meta.dictionary['Name'])+'</a>\n'
                text += '</h3>\n<br>\n'
                if meta.dictionary.has_key('Organisation'):
                    text += self.addSimpleComponentInfo('Organisation', unicode(meta.dictionary['Organisation']))
                
                text += self.addSimpleComponentInfo('Project', unicode(meta.dictionary['Project']))
                text += self.addSimpleComponentInfo('Location', unicode(meta.dictionary['Location']))
                if meta.dictionary.has_key('Project Web'):
                    text += self.addSimpleComponentInfo('Project Web', unicode(meta.dictionary['Project Web']))
                text += self.addSimpleComponentInfo('Copyright', unicode(meta.dictionary['Copyright']))
                text += self.addSimpleComponentInfo('License', unicode(meta.dictionary['License']))
                #This one may need work
                if meta.dictionary.has_key('Parent'):
                    text += self.addParentInfoDtd('Parent', 'Name', unicode(meta.dictionary['Parent']))
                
                text += self.addChildrenInfoDtd('Children', 'Parent', 'Name', (meta.dictionary['Name']).encode('ascii','xmlcharrefreplace'))
                
                text += self.addSimpleComponentInfo('Description', (meta.dictionary['Description']).encode('ascii','xmlcharrefreplace'))
                # Next comes reference, if any
                if meta.dictionary.has_key('Reference'):
                    refValue = unicode(meta.dictionary['Reference'])
                    text += self.addSimpleComponentInfo('Reference', refValue.encode('ascii','xmlcharrefreplace'))
                #next comes Example, if any
                if meta.dictionary.has_key('Example'):
                    text +=self.addExampleInfo('Example',  (meta.dictionary['Example']).encode('ascii','xmlcharrefreplace'))
                #Then equation if any
                if meta.dictionary.has_key('Equation'):
                    text +=self.addEquationInfo(unicode(meta.dictionary['Name']),  (meta.dictionary['Equation']).encode('ascii','xmlcharrefreplace'))
                if meta.dictionary.has_key('Params'):
                # then params if any
                    text +=self.addListComponentInfo('Params',  'paramtable', meta.dictionary['Params'],  meta.dictionary['Name'])
                # then dependencies, if any
                if meta.dictionary.has_key('Dependencies'):
                    text +=self.addListComponentInfo("Dependencies", 'deptable',  meta.dictionary['Dependencies'],  meta.dictionary['Name'])
        return unicode(text)        
        
    ## Find the parent info for meta component, using xsd.    
    def addParentInfoXsd(self, titleName,  infoTitle, nameTitle, parentName):
        myString  = ''
        #search through all projects to find meta with parentName.
        for project in self.projectList:
            for meta in project.metas:
                name = meta.dictionary[infoTitle][nameTitle]
                
                if name == parentName:
                    #If parent found, add as hyperlink
                    myString += '<a href="./'+ project.name+'.html#'+name+'">'+name+'<br/></a>'
        text  = self.addSimpleComponentInfo(titleName,myString)
        return text
        
    ## Find all children of meta component, usiing xsd    
    def addChildrenInfoXsd(self,  titleName,  codeTitle, parentTitle, infoTitle, nameTitle, parentName):
        myString  = '<ul>\n'
        childExists = False
        #search through all metas to find meta that list parentName as their parent.
        
        for project in self.projectList:
            for meta in project.metas:
                if meta.dictionary[codeTitle].has_key(parentTitle):
                    name = meta.dictionary[codeTitle][parentTitle]
                else:
                    name ="" 
                if name == parentName:
                    childExists = True
                    #If child found, add as hyperlink
                    myString += '<li><a href="./'+ project.name+'.html#'+meta.dictionary[infoTitle][nameTitle]+'">'+meta.dictionary[infoTitle][nameTitle]+'<br/></a></li>\n'
        myString +='</ul>'
        if childExists ==True:
            text  = self.addSimpleComponentInfo(titleName,myString)
        else:
            text = ""
        return text

    ## Create individual meta entry, using new Xsd format.
    def createMetaEntryXsd(self,  name,  count):
        #create individual meta entry
        text = ''
        project = self.projectList[count]
        #Select Meta
        for meta in project.metas:
            if meta.dictionary['info']['title'] == name:
                # Add meta data to file
                #print project.name,  name,  meta.dictionary.keys()
                text +='<h3 align="center">\n'
                
                text +='<a name="'+meta.dictionary['info']['title']+'">'+meta.dictionary['info']['title']+'</a>\n'
                text += '</h3>\n<br>\n'
                text += '<h4 align="center">Information</h4>\n'
                if meta.dictionary['info'].has_key('creator'):
                    text += self.addSimpleComponentInfo('creator', meta.dictionary['info']['creator'])                
                text += self.addSimpleComponentInfo('subject', meta.dictionary['info']['subject'])
                text += self.addSimpleComponentInfo('source', meta.dictionary['info']['source'])
                if meta.dictionary['info'].has_key('publisher'):
                    text += self.addSimpleComponentInfo('publisher', meta.dictionary['info']['publisher'])
                text += self.addSimpleComponentInfo('rights', meta.dictionary['info']['rights'])    
                text += self.addSimpleComponentInfo('description', meta.dictionary['info']['description'])
                if meta.dictionary['implements'].has_key('equation') or meta.dictionary['implements'].has_key('reference'): 
                    text += '<h4 align="center">Implements</h4>\n'
                #Then equation and reference if any
                if meta.dictionary['implements'].has_key('equation'):
                    text +=self.addEquationInfo(meta.dictionary['info']['title'],  meta.dictionary['implements']['equation'])
                if meta.dictionary['implements'].has_key('reference'):
                    text += self.addSimpleComponentInfo('reference', (meta.dictionary['implements']['reference']).encode('ascii','xmlcharrefreplace')) 
                text += '<h4 align="center">Code</h4>\n'
                #This one may need work
                if meta.dictionary['code'].has_key('inherits'):
                    text += self.addParentInfoXsd('inherits', 'info','title', meta.dictionary['code']['inherits'])
                
                text += self.addChildrenInfoXsd('Children', 'code','inherits', 'info','title', meta.dictionary['info']['title'])
                
                #next comes Example, if any
                if meta.dictionary['code'].has_key('example-code'):
                    text +=self.addExampleInfo('Example Code',  meta.dictionary['code']['example-code'])
                

                if meta.dictionary.has_key('parameters'):
                # then params if any
                    text +=self.addListComponentInfo('parameters',  'paramtable', meta.dictionary['parameters'],  meta.dictionary['info']['title'])
                # then dependencies, if any
                if meta.dictionary.has_key('associations'):
                    text +=self.addListComponentInfo("associations", 'deptable',  meta.dictionary['associations'],  meta.dictionary['info']['title'])
        return text        

    ## Function adds any scripts needed to be embedded in webpages. This code does not check
    # if it is genuine java code, just embeds it in \<script\> tags.
    def addScripts(self):
        #add scripts
        text = ''
        for script in self.scripts:
            if os.path.isfile(script):
                text += '<!-- This is the script from file: '+script+' -->\n'
                text += '<Script type="text/javascript">\n'
                #open script file and read content.
                f = open(script, 'r')
                scriptString = f.read()
                f.close()
                #put content in script tag
                text += scriptString
                text += "\n</Script>\n"
            else:
                print "" + script + " is not a valid file. Not adding."
        return text
    
    ## Copy pictures in pictures directory to new pictures directory. Actually just copies ALL
    # content from that directory to pictures sub-directory.
    def copyPictures(self,  picPath):
        if not(picPath == ""):
            print "Copying directory of pictures :" + os.path.realpath(picPath) + "/* to directory: " + os.path.join(self.path,  self.pictureSubDirectoryName) 
            # Copy everything there to the moment
            os.system('cp '+ os.path.realpath(picPath) + "/* "+ os.path.join(self.path,  self.pictureSubDirectoryName) )
        else:
            print "No Picture directory to copy."
    ## Copy stylesheet files using given paths to the base documentPath directory.       
    def copyStylesheets(self):
        # copy stylesheet from current location to documentPath.
        print "Copying Stylesheets from current locations to " + self.path
        for stylesheet in self.styleSheetList:
            #strip out path and name

            (path,  file) = os.path.split(os.path.realpath(stylesheet))
            newPath = os.path.join(self.path,  file)

            os.system('cp '+ os.path.realpath(stylesheet)+" " + newPath)
            
## DivIds class. List of names associated with types for creating \<div\> tags in html
# This should be converted to a dictionary later for cgreater flexibility in adding
# webpage div types.
class DivIds():
    def __init__(self):
        self.body ='body'
        self.main  = ''
        self.head = 'head'
        self.container =''
        self.toptab = ''
        self.breadcrumbtrail = ''
        self.navitab = []
        self.desc = ''
        self.sidebar = ''
        self.sidebarMenu = ''
        self.sidebarTiers =[]
        self.footer = ''
        self.componentInfo =[]
        self.params =[]
        self.dependencies =[]
        self.codebox =''
    ## Function that assigns names to div types.    
    def createDivIds(self,  projectName):
        # TODO: eventually might make this search the stylefile for values, or input from external file or argv,
        # For present, just set them.
        self.body ='body'
        self.main  = 'main'
        self.head = 'head'
        self.toptab = 'toptab'
        self.breadcrumbtrail = 'breadcrumbtrail'
        self.navitab = ['navitab',  'activenavitab']
        self.container ='container'
        self.sidebar = 'sidebar'
        self.sidebarMenu = 'menu'
        self.sidebarTiers =['tier1sidelink',  'sidelink']
        self.footer = 'footer'
        self.componentInfo =['componentInfo',  'componentInfoTitle',  'componentInfoData']
        self.params =['params',  'paramtable']
        self.dependencies =['dependencies',  'deptable']
        self.codebox ='codebox'       
        
        if str(projectName) == '':
            self.desc = 'desc'
        else:
            self.desc = 'desc-' + str(projectName)
            
## Preprocesses input data so that it can be converted into a useful python list.         
def StripTextToList(text,  tier):
    myList = []
    if tier == 1:
        myList = string.split(string.lstrip(string.rstrip(text, ']'), '['), ',')
        for count in range(len(myList)):
            myList[count] = str(string.lstrip(string.rstrip(myList[count]) ))
            myList[count] = str(string.lstrip(string.rstrip(myList[count], "'"), "'"))
            myList[count] = str(string.lstrip(string.rstrip(myList[count], '"'), '"'))

    elif tier ==2:
        values = string.split(string.lstrip(string.rstrip(text, "]"), "["), ',[')
    
        for  value in values:
            myList.append(string.split(string.rstrip(value,  ']'), ','))
    
        for count in range(len(myList)):
            for value in range(len(myList[count])):
                myList[count][value] = str(string.lstrip(string.rstrip(myList[count][value])))
                myList[count][value] = str(string.lstrip(string.rstrip(myList[count][value], "'"), "'"))
                myList[count][value] = str(string.lstrip(string.rstrip(myList[count][value], '"'), '"'))
    else:
        sys.exit("Level Number: " +str(tier)+" invalid option for list function" )

    return myList

## Main function run.
if __name__=='__main__':
    from optparse import OptionParser
    #Option parser gives a nice interface to pick the directory to read from
    #TODO: might in future add an option to specify name suffix to search for.
    parser = OptionParser()
    parser.add_option("-t", "--type", dest="type",
                  help="file type to create. options: html",  default="html")
    parser.add_option("-p", "--path", dest="path",
                  help="path to project source/s. Default: '../'",  default="../")
    parser.add_option("-o", "--outputPath", dest="documentPath",
                  help="location to write documents. Default: './DOCS'",  default="./DOCS")
    parser.add_option("-P", "--Projects", dest="Projects",
                  help="List of projects to parse. Should be in format: ['name', 'name'] . Default: All current Underworld Project.",  
                  default="['StGermain','StgDomain','StgFEM', 'PICellerator','Underworld', 'gLucifer']")                  
    parser.add_option("-n", "--documentName", dest="documentName",
                  help="Overall name of document/s. Default:'Component Codex'",  default="Component Codex")
    parser.add_option("-b", "--documentBlurb", dest="documentBlurb",
                  help="General blurb to appear in blurb section on each page. Default:'This is a list of the available components.' ",  default="This is a list of the available components.")
                                    
    parser.add_option("-s", "--stylesheetList", dest="stylesheetList",
                  help="list of stylesheets to use, with path. Example: ['./print.css','./more.css'] , Default: ['./print.css']",  default="['./print.css']")
    parser.add_option("-w", "--extweb", dest="externalWebpages",
                  help="list of external webpages to be listed at top of document. format: [[page, title]] . Default: [['http://www.underworldproject.org','Underworld Home Page']]",  
                  default="[['http://www.underworldproject.org','Underworld Home Page']]")  
    parser.add_option("-e", "--email", dest="email",
                  help="single email contact address, Example(default): underworld-users@vpac.org",  default="underworld-users@vpac.org")
    parser.add_option("-c", "--scripts", dest="scripts",
                  help="filename of file containing any java scripts for the webpage, with path. No Default. Example: ['./myscript.js', 'myscript2.js']",  default="")                  
    parser.add_option("-m", "--metaType", dest="metaFlag",
                  help="type of meta file to output. options: 'dtd' (classic StGermain meta), or 'xsd' (new meta format), Default: 'dtd'", 
                  default="dtd")
    parser.add_option("-I", "--imagePath", dest="picPath",
                  help="Path for location of already existing images to include in webpages. No Default.", 
                  default="")    
    parser.add_option("-y", "--pictureDirectory", dest="pictureDirectory",
                  help="Name of picture directory to be created inside output directory for all pictures. Default: 'images'", 
                  default="images")
    parser.add_option("-x", "--indexFileName", dest="indexFileName",
                  help="Path to file that contains html content for index.html for webpages. No Default.", 
                  default="images") 
    parser.add_option("-S", "--searchBarType", dest="searchBarType",
                  help="Type of search bar to use on Index Page for searching for components on other pages. Options: 'Google', 'None'. Default: 'Google'", 
                  default="Google")      
    parser.add_option("-T", "--searchPath", dest="searchPath",
                  help="Path name to use on search bar if different from outputPath. Useful for when the output will be copied into a different directory later, or put online.  No Default.", 
                  default="")                    
    (options, args) = parser.parse_args()
    
    #Now create the projects.
    projectList = []
    projectNames =StripTextToList(options.Projects, 1)
    for projectName  in projectNames:
        # check project name in current path already
        if (string.count(str(os.path.realpath(options.path)),  str("/"+projectName)) > 0) :
           path =  os.path.realpath(options.path) + "/"
        # check if adding projectName
        else:
            path = os.path.realpath(options.path) + "/"+ str(projectName) + "/"
        print "Creating data for " + projectName
        project =Project(projectName, path,  options.metaFlag)
        project.assignMetas()
        projectList.append(project)
        
    #now create the HTML documents
    print "Now creating HTML documents in directory: " + os.path.realpath(options.documentPath)
    #turn Ext.webs into useful list
    extWeb = StripTextToList(options.externalWebpages, 2)
    stylesheetList = StripTextToList(options.stylesheetList, 1)
    
    if not(options.scripts == ""):
        scripts = StripTextToList(options.scripts, 1)
    else:
        scripts = []
    htmlDocuments = createHTMLDocuments(options.documentName, options.documentBlurb,   projectList,  options.documentPath,stylesheetList,  extWeb,  options.email, scripts,  options.pictureDirectory,  options.indexFileName,  options.searchBarType, options.searchPath)
    htmlDocuments.createHTMLPages()    
    htmlDocuments.copyPictures(options.picPath)
    htmlDocuments.copyStylesheets()

