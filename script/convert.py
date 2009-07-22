
import stgDtd
import stgMetaXsd

def dtd2MetaXsd( xml_text ):
	# Parse DTD
	try:
		dtd = stgDtd.readXML( xml_text )
	except:
		print 'Failed to parse as a StGermain DTD'
		raise

	# Convert DTD-dict to XSD-dict
	try:
		xsd = dtdDict2metaXsdDict( dtd )
	except:
		print 'Failed to convert information from a StGermain Meta DTD to a StGermain Meta XSD'
		raise

	# Write XSD
	try:
		doc = stgMetaXsd.createXML( xsd )
		return doc.toprettyxml()
	except:
		print 'Failed to write to StGermain Meta XSD'
		raise


def dtdDict2metaXsdDict( stgDtdDict ):
	meta = {}

	# Info...
	meta["info"] = {}
	meta["info"]["title"] = stgDtdDict["Name"]
	# Often Author has been left blank... use Organisation where this has occured
	try:
		if stgDtdDict["Author"] != '...':
			meta["info"]["creator"] = stgDtdDict["Author"]
		else:
			meta["info"]["creator"] = stgDtdDict["Organisation"]
	except KeyError:
		meta["info"]["creator"] = stgDtdDict["Organisation"]
	meta["info"]["publisher"] = stgDtdDict["Organisation"]
	meta["info"]["rights"] = stgDtdDict["License"]
	meta["info"]["source"] = stgDtdDict["Location"]
	meta["info"]["subject"] = stgDtdDict["Project"]
	meta["info"]["description"] = stgDtdDict["Description"]

	# Code...
	meta["code"] = {}
	#meta["code"]["example-documentation"] = 
	try:
		if stgDtdDict["Example"] != '...':
			meta["code"]["example-code"] = stgDtdDict["Example"]
	except KeyError:
		pass
	try:
		if stgDtdDict["Parent"] != '...':
			meta["code"]["inherits"] = stgDtdDict["Parent"]
	except KeyError:
		pass

	# Implements...
	meta["implements"] = {}
	try:
		if stgDtdDict["Reference"] != '...':
			meta["implements"]["reference"] = stgDtdDict["Reference"]
	except KeyError:
		pass
	try:
		if stgDtdDict["Equation"] != '...':
			meta["implements"]["equation"] = stgDtdDict["Equation"]
	except KeyError:
		pass

	# Parameters...
	meta["parameters"] = []
	try:
		for param in stgDtdDict["Params"]:
			p = {}
			p["name"] = param["Name"]
			# Convert value...
			if param["Type"].lower().strip() == "double":
				p["type"] = "xsd:double"
			elif param["Type"].lower().strip() == "int" or param["Type"].lower().strip() == "integer":
				p["type"] = "xsd:int"
			elif param["Type"].lower().strip() == "long":
				p["type"] = "xsd:long"
			elif param["Type"].lower().strip() == "char":
				p["type"] = "xsd:byte"
			elif param["Type"].lower().strip() == "bool":
				p["type"] = "xsd:boolean"
			elif param["Type"].lower().strip() == "string":
				p["type"] = "xsd:string"
			elif param["Type"].lower().strip() == "unsignedint" or param["Type"].lower().strip() == "unsigned int":
				p["type"] = "xsd:unsignedInt"
			elif param["Type"].lower().strip() == "unsignedlong" or param["Type"].lower().strip() == "unsigned long":
				p["type"] = "xsd:unsignedLong"
			elif param["Type"].lower().strip() == "float":
				p["type"] = "xsd:float"
			elif param["Type"].lower().strip() == "list":
				p["type"] = "stg:list"
			else:
				raise RuntimeError( 'Unknown parameter type' )
			try:
				p["default"] = param["Default"]
			except KeyError:
				pass
			try:
				p["documentation"] = param["Description"]
			except KeyError:
				pass
			meta["parameters"].append( p )
	except KeyError:
		pass

	# Associations...
	meta["associations"] = []
	try:
		for assoc in stgDtdDict["Dependencies"]:
			a = {}
			a["name"] = assoc["Name"]
			a["type"] = assoc["Type"]
			# Take opporunity to enforce a nillability state (even though XML XSD doesn't require it). Assume nillable=false by default.
			try:
				if assoc["Essential"].lower().strip() == "false" or assoc["Essential"].lower().strip() == "f" or assoc["Essential"].lower().strip() == "no":
					a["nillable"] = "true"
				else:
					a["nillable"] = "false"
			except KeyError:
				a["nillable"] = "false"
			try:
				a["documentation"] = assoc["Description"]
			except KeyError:
				pass
			meta["associations"].append( a )
	except KeyError:
		pass

	return meta

## The purpose of this function is to provide the C code for the includes needed by the code generated by 
# "metaXsdDict2stgDictionaryCode" and "metaXsdDict2stgStrings"
def metaXsdDict2stgCodeHeader():
	s = ''
	s += '#include <stdarg.h>\n'
	# Purposely refer to internal StGermain headers, as StGermain/StGermain.h would necessarily exist yet if a meta file exists
	# within StGermain itself.
	s += '#include "StGermain/Base/Foundation/Foundation.h"\n'
	s += '#include "StGermain/Base/IO/IO.h"\n'

	return s

## The purpose of this function is to write C code that creates a const char* for the component type name and the xml meta data
def metaXsdDict2stgStrings( xsdDict ):
	s = ''

	# The Name of the component
	s += 'const char* ' + safecvar( xsdDict["info"]["title"] ) +  '_Name = "' + safecval( xsdDict["info"]["title"] ) + '";\n'
	s += 'const char* ' + safecvar( xsdDict["info"]["title"] ) +  '_GetName() {\n'
	s += '\treturn ' + safecvar( xsdDict["info"]["title"] ) +  '_Name;\n'
	s += '}\n'
	s += '\n'

	# The xml of the meta of the component
	xsdDoc = stgMetaXsd.createXML( xsdDict )
	xsdTxt = xsdDoc.toprettyxml()
	xsdTxt = safecval( xsdTxt )	
	s += 'const char* ' + safecvar( xsdDict["info"]["title"] ) +  '_Meta = "' + xsdTxt + '";\n'
	s += 'const char* ' + safecvar( xsdDict["info"]["title"] ) +  '_GetMetadata() {\n'
	s += '\treturn ' + safecvar( xsdDict["info"]["title"] ) +  '_Meta;\n'
	s += '}\n'
	# The _Type variant exists because of macro used for ComponentRegister_Add does a stringify on the Component_Type argument
	s += 'const char* ' + safecvar( xsdDict["info"]["title"] ) +  '_Type_GetMetadata() {\n'
	s += '\treturn ' + safecvar( xsdDict["info"]["title"] ) +  '_Meta;\n'
	s += '}\n'

	return s
	

## The purpose of this function is to write C code that creates a StGermain Dictionary that is equivalent to the passed in python
# dictionary representing a StGermain meta Xsd
def metaXsdDict2stgDictionaryCode( xsdDict ):
	s = ''
	s += 'Dictionary* ' + safecval( xsdDict["info"]["title"] ) + '_MetaAsDictionary() {\n'
	s += '\tDictionary* meta;\n'
	s += '\tDictionary* info;\n'
	s += '\tDictionary* code;\n'
	s += '\tDictionary* implements;\n'
	s += '\tDictionary* parameters;\n'
	for param in xsdDict["parameters"]:
		s+= '\tDictionary* ' + safecvar( param["name"] ) + 'Param;\n'
	s += '\tDictionary* associations;\n'
	for assoc in xsdDict["associations"]:
		s+= '\tDictionary* ' + safecvar( assoc["name"] ) + 'Assoc;\n'
	s += '\n'

	s += '\tmeta = Dictionary_New();\n'
	s += '\n'

	# XML ... (requires metaXsdDict2stgStrings( xsdDict ) to have been called first)
	s += '\tDictionary_Add( meta, "xml", Dictionary_Entry_Value_FromString( ' + safecvar( xsdDict["info"]["title"] ) + '_Meta ));\n'
	s += '\n'
	
	# Info...
	s += '\tinfo = Dictionary_New();\n'
	s += '\tDictionary_Add( info, "title", Dictionary_Entry_Value_FromString( "' + safecval( xsdDict["info"]["title"] ) + '" ));\n'
	s += '\tDictionary_Add( info, "creator", Dictionary_Entry_Value_FromString( "' + safecval( xsdDict["info"]["creator"] ) + '" ));\n'
	s += '\tDictionary_Add( info, "publisher", Dictionary_Entry_Value_FromString( "' + safecval( xsdDict["info"]["publisher"] ) + '" ));\n'
	s += '\tDictionary_Add( info, "rights", Dictionary_Entry_Value_FromString( "' + safecval( xsdDict["info"]["rights"] ) + '" ));\n'
	s += '\tDictionary_Add( info, "source", Dictionary_Entry_Value_FromString( "' + safecval( xsdDict["info"]["source"] ) + '" ));\n'
	s += '\tDictionary_Add( info, "subject", Dictionary_Entry_Value_FromString( "' + safecval( xsdDict["info"]["subject"] ) + '" ));\n'
	s += '\tDictionary_Add( info, "description", Dictionary_Entry_Value_FromString( "' + safecval( xsdDict["info"]["description"] ) + '" ));\n'
	s += '\tDictionary_Add( meta, "info", Dictionary_Entry_Value_FromStruct( info ) );\n'
	s += '\n'

	# Code...
	s += '\tcode = Dictionary_New();\n'
	try:
		s += '\tDictionary_Add( code, "example-documentation", Dictionary_Entry_Value_FromString( "' + safecval( xsdDict["code"]["example-documentation"] ) + '" ));\n'
	except KeyError:
		pass
	try:
		s += '\tDictionary_Add( code, "example-code", Dictionary_Entry_Value_FromString( "' + safecval( xsdDict["code"]["example-code"] ) + '" ));\n'
	except KeyError:
		pass
	try:
		s += '\tDictionary_Add( code, "inherits", Dictionary_Entry_Value_FromString( "' + safecval( xsdDict["code"]["inherits"] ) + '" ));\n'
	except KeyError:
		pass
	s += '\tDictionary_Add( meta, "code", Dictionary_Entry_Value_FromStruct( code ) );\n'
	s += '\n'

	# Implements...
	s += '\timplements = Dictionary_New();\n'
	try:
		s += '\tDictionary_Add( implements, "reference", Dictionary_Entry_Value_FromString( "' + safecval( xsdDict["implements"]["reference"] ) + '" ));\n'
	except KeyError:
		pass
	try:
		# BIG ASSUMPTION: equation is in latex
		s += '\tDictionary_Add( implements, "equation", Dictionary_Entry_Value_FromString( "' + safecvalFromLatex( xsdDict["implements"]["equation"] ) + '" ));\n'
	except KeyError:
		pass
	s += '\tDictionary_Add( meta, "implements", Dictionary_Entry_Value_FromStruct( implements ) );\n'
	s += '\n'

	# Parameters...
	s += '\tparameters = Dictionary_New();\n'
	for param in xsdDict["parameters"]:
		s += '\t' + safecvar( param["name"] ) + 'Param = Dictionary_New();\n'
		s += '\tDictionary_Add( ' + safecvar( param["name"] ) + 'Param, "name", Dictionary_Entry_Value_FromString( "' + safecval( param["name"] ) + '" ));\n'
		s += '\tDictionary_Add( ' + safecvar( param["name"] ) + 'Param, "type", Dictionary_Entry_Value_FromString( "' + safecval( param["type"] ) + '" ));\n'
		try:
			s += '\tDictionary_Add( ' + safecvar( param["name"] ) + 'Param, "default", Dictionary_Entry_Value_FromString( "' + safecval( param["default"] ) + '" ));\n'
		except KeyError:
			pass
		try:
			s += '\tDictionary_Add( ' + safecvar( param["name"] ) + 'Param, "documentation", Dictionary_Entry_Value_FromString( "' + safecval( param["documentation"] ) + '" ));\n'
		except KeyError:
			pass
		s += '\tDictionary_Add( parameters, "' + safecval( param["name"] ) + '", Dictionary_Entry_Value_FromStruct( ' + safecvar( param["name"] ) + 'Param ) );\n'
		s += '\n'
	s += '\tDictionary_Add( meta, "parameters", Dictionary_Entry_Value_FromStruct( parameters ) );\n'
	s += '\n'

	# Associations...
	s += '\tassociations = Dictionary_New();\n'
	for assoc in xsdDict["associations"]:
		s += '\t' + safecvar( assoc["name"] ) + 'Assoc = Dictionary_New();\n'
		s += '\tDictionary_Add( ' + safecvar( assoc["name"] ) + 'Assoc, "name", Dictionary_Entry_Value_FromString( "' + safecval( assoc["name"] ) + '" ));\n'
		s += '\tDictionary_Add( ' + safecvar( assoc["name"] ) + 'Assoc, "type", Dictionary_Entry_Value_FromString( "' + safecval( assoc["type"] ) + '" ));\n'
		try:
			s += '\tDictionary_Add( ' + safecvar( assoc["name"] ) + 'Assoc, "nillable", Dictionary_Entry_Value_FromString( "' + safecval( assoc["nillable"] ) + '" ));\n'
		except KeyError:
			pass
		try:
			s += '\tDictionary_Add( ' + safecvar( assoc["name"] ) + 'Assoc, "documentation", Dictionary_Entry_Value_FromString( "' + safecval( assoc["documentation"] ) + '" ));\n'
		except KeyError:
			pass
		s += '\tDictionary_Add( associations, "' + safecval( assoc["name"] ) + '", Dictionary_Entry_Value_FromStruct( ' + safecvar( assoc["name"] ) + 'Assoc ) );\n'
		s += '\n'
	s += '\tDictionary_Add( meta, "associations", Dictionary_Entry_Value_FromStruct( associations ) );\n'
	s += '\n'
	s += '\treturn meta;\n'
	s += '}\n'

	# The _Type variant exists because of macro used for ComponentRegister_Add does a stringify on the Component_Type argument
	s += 'Dictionary* ' + safecval( xsdDict["info"]["title"] ) + '_Type_MetaAsDictionary() {\n'
	s += '\treturn ' + safecvar( xsdDict["info"]["title"] ) +  '_MetaAsDictionary();\n'
	s += '}\n'

	return s


## Convert the python unicode string 's' into an ascii string that is safe to put in a C string as a value
def safecval( s ):
	return s.replace( "\t", "\\t" ).replace( "\n", "\\n" ).replace( "\"", "\\\"" ).encode( 'ascii', 'replace' )


## Convert the python unicode string 's' into as ascii string that is safe to use as a C variable
def safecvar( s ):
	# The followign doesn't work because we use "_" as quasi namespacing for component names:)	
	#r = ''
	#for a in s.encode( 'ascii', 'replace' ):
	#	if a.isalnum():
	#		r += a
	#return r
	return s.replace( " ", "" ).replace( "-", "" ).replace( "\t", "" ).replace( "\n", "" ).replace( "\"", "" ).replace( "&", "" ).replace( "[", "" ).replace( "]", "" ).replace( "<", "" ).replace( ">", "" ).encode( 'ascii', 'ignore' )

## Convert the python unicode string 's' that is of latex into an ascii string that is safe to put in a C string as a value
def safecvalFromLatex( s ):
	return (safecval( s ).replace( "\\", "\\\\" )).encode('ascii','replace')



