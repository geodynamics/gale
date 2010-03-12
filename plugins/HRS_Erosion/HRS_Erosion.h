/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Victorian Partnership for Advanced Computing (VPAC) Ltd, Australia
** (C) 2003-2005 All Rights Reserved
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	David May, PhD Student Monash University, VPAC. (david.may@sci.maths.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**
**
** Role:
**	Defines any header info, such as new structures, needed by this plugin
**
** Assumptions:
**
** Comments:
**
** $Id $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Underworld_plugins_HRS_Erosion_HRS_Erosion_h__
#define __Underworld_plugins_HRS_Erosion_HRS_Erosion_h__

	Index Underworld_HRS_Erosion_Register( PluginsManager* pluginsMgr );

	void* _Underworld_HRS_Erosion_DefaultNew( Name name );

	void _Underworld_HRS_Erosion_AssignFromXML( void* codelet, Stg_ComponentFactory* cf, void* data );

	void _Underworld_HRS_Erosion_Build( void* codelet, void* data );

	void _Underworld_HRS_Erosion_Destroy( void* codelet, void* data );

#endif /* __Underworld_plugins_HRS_Erosion_HRS_Erosion_h__ */
