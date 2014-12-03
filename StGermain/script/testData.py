
def stgMetaDtdDict():
	meta = {}

	# Info...
	meta["Name"] = "Arrhenius"
	meta["Author"] = "Monash University"
	meta["Organisation"] = "Monash University"
	meta["License"] = "http://www.opensource.org/licenses/bsd-license.php"
	meta["Location"] = "Underworld/Rheology/src"
	meta["Project"] = "Underworld"
	meta["Description"] = "A temperature dependent rheology."

	# Code...
	meta["Example"] = "<Arrhenius>\n\t<TemperatureField> temperatureField </TemperatureField>\n\t<eta0> 1.0e-6 </eta0>\n</Arrhenius>\n"
	meta["Parent"] = "Rheology"

	# Implements...
	meta["Reference"] = "..."
	meta["Equation"] = "$\\eta = \\eta_0 \exp \\left( \\frac {E + V d} {T + T_0} \\right)$"

	# Deprecated...
	meta["ProjectWeb"] = "http://www.mcc.monash.edu.au/Software/Underworld"
	meta["Copyright"] = "Copyright (c) 2005, Monash Cluster Computing"
	meta["Summary"] = "..."

	# Parameters...
	meta["Params"] = []
	p0 = {}
	p0["Name"] = "eta0"
	p0["Type"] = "Double"
	p0["Default"] = "1.0"
	p0["Description"] = "This is the $\\eta_0$ in the equation above."
	meta["Params"].append( p0 )

	# Associations...
	meta["Dependencies"] = []
	a0 = {}
	a0["Name"] = "TemperatureField"
	a0["Type"] = "FeVariable"
	a0["Essential"] = "True"
	a0["Description"] = "The TemperatureField that provides the $T$ in the equation above."
	meta["Dependencies"].append( a0 )

	return meta


def stgMetaXsdDict():
	meta = {}

	# Info...
	meta["info"] = {}
	meta["info"]["title"] = "Arrhenius"
	meta["info"]["creator"] = "Monash University"
	meta["info"]["publisher"] = "Monash University"
	meta["info"]["rights"] = "http://www.opensource.org/licenses/bsd-license.php"
	meta["info"]["source"] = "Underworld/Rheology/src"
	meta["info"]["subject"] = "Underworld"
	meta["info"]["description"] = "A temperature dependent rheology."

	# Code...
	meta["code"] = {}
	#meta["code"]["example-documentation"] = "This example demonstrates setting the 'parameters' and 'associations'"
	meta["code"]["example-code"] = "<Arrhenius>\n\t<TemperatureField> temperatureField </TemperatureField>\n\t<eta0> 1.0e-6 </eta0>\n</Arrhenius>\n"
	meta["code"]["inherits"] = "Rheology"

	# Implements...
	meta["implements"] = {}
	meta["implements"]["reference"] = "..."
	meta["implements"]["equation"] = "$\\eta = \\eta_0 \exp \\left( \\frac {E + V d} {T + T_0} \\right)$"

	# Parameters...
	meta["parameters"] = []
	p0 = {}
	p0["name"] = "eta0"
	p0["type"] = "xsd:double"
	p0["default"] = "1.0"
	p0["documentation"] = "This is the $\\eta_0$ in the equation above."
	meta["parameters"].append( p0 )

	# Associations...
	meta["associations"] = []
	a0 = {}
	a0["name"] = "TemperatureField"
	a0["type"] = "FeVariable"
	a0["nillable"] = "false"
	a0["documentation"] = "The TemperatureField that provides the $T$ in the equation above."
	meta["associations"].append( a0 )

	return meta

def stgMetaXsdXml():
	return u'<?xml version="1.0" ?>\n<meta xmlns="urn:stgermainmeta-schema" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://purl.org/dc/elements/1.1/ dc.xsd http://www.w3.org/2001/XMLSchema XMLSchema.xsd urn:stgermainmeta-schema stgermainmeta.xsd">\n\t<info>\n\t\t<dc:title>\n\t\t\tArrhenius\n\t\t</dc:title>\n\t\t<dc:creator>\n\t\t\tMonash University\n\t\t</dc:creator>\n\t\t<dc:publisher>\n\t\t\tMonash University\n\t\t</dc:publisher>\n\t\t<dc:rights>\n\t\t\thttp://www.opensource.org/licenses/bsd-license.php\n\t\t</dc:rights>\n\t\t<dc:source>\n\t\t\tUnderworld/Rheology/src\n\t\t</dc:source>\n\t\t<dc:subject>\n\t\t\tUnderworld\n\t\t</dc:subject>\n\t\t<dc:description>\n<![CDATA[A temperature dependent rheology.]]>\t\t</dc:description>\n\t</info>\n\t<code>\n\t\t<xsd:annotation>\n\t\t\t<xsd:documentation/>\n\t\t\t<xsd:appinfo>\n<![CDATA[<Arrhenius>\n\t<TemperatureField> temperatureField </TemperatureField>\n\t<eta0> 1.0e-6 </eta0>\n</Arrhenius>\n]]>\t\t\t</xsd:appinfo>\n\t\t</xsd:annotation>\n\t\t<inherits>\n\t\t\tRheology\n\t\t</inherits>\n\t</code>\n\t<implements>\n\t\t<reference>\n\t\t\t...\n\t\t</reference>\n\t\t<equation>\n<![CDATA[$\\eta = \\eta_0 \\exp \\left( \\frac {E + V d} {T + T_0} \\right)$]]>\t\t</equation>\n\t</implements>\n\t<parameters>\n\t\t<xsd:element default="1.0" name="eta0" type="xsd:double">\n\t\t\t<xsd:annotation>\n\t\t\t\t<xsd:documentation>\n<![CDATA[This is the $\\eta_0$ in the equation above.]]>\t\t\t\t</xsd:documentation>\n\t\t\t</xsd:annotation>\n\t\t</xsd:element>\n\t</parameters>\n\t<associations>\n\t\t<xsd:element name="TemperatureField" nillable="false" type="FeVariable">\n\t\t\t<xsd:annotation>\n\t\t\t\t<xsd:documentation>\n<![CDATA[The TemperatureField that provides the $T$ in the equation above.]]>\t\t\t\t</xsd:documentation>\n\t\t\t</xsd:annotation>\n\t\t</xsd:element>\n\t</associations>\n</meta>\n'

