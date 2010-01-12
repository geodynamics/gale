/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Lesser General Public
**  License as published by the Free Software Foundation; either
**  version 2.1 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Lesser General Public License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License along with this library; if not, write to the Free Software
**  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
**
*/
/** \file
**  Role:
**
** Assumptions:
**
** Comments:
**
** $Id: Plugin.h 3192 2005-08-25 01:45:42Z AlanLo $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StGermain_Base_Extensibility_Plugin_h__
#define __StGermain_Base_Extensibility_Plugin_h__
	
	/** The prototype for the Register function in a plugin */
	typedef Index (Plugin_RegisterFunction) ( void* pluginsManager );

	/* Textual name of this class */
	extern const Type Plugin_Type;

	/* Plugins info */
	#define __Plugin \
		/* General info */ \
		__Module \
		\
		/* Virtual info */ \
		\
		/* Plugin info */ \
		Plugin_RegisterFunction*    Register;
		
	struct Plugin { __Plugin };


	/* Create a new Plugin */
	Plugin* Plugin_New( Name name, Stg_ObjectList* directories );
	Module* Plugin_Factory( Name name, Stg_ObjectList* directories );
	
	/* Creation implementation / Virtual constructor */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define PLUGIN_DEFARGS \
                MODULE_DEFARGS

	#define PLUGIN_PASSARGS \
                MODULE_PASSARGS

	Plugin* _Plugin_New(  PLUGIN_DEFARGS  );
	
	/* Initialisation implementation */
	void _Plugin_Init( Plugin* self );
	
	/* Delete implementation */
	void _Plugin_Delete( void* plugin );
	
	/* Print implementation */
	void _Plugin_Print( void* plugin, Stream* stream );

	/* MangleName implementation */
	char* _Plugin_MangleName( char* name );

	/** Get the function pointer the to the plugin's register function */
	Plugin_RegisterFunction* Plugin_GetRegisterFunc( void* plugin );
	
#endif /* __StGermain_Base_Extensibility_Plugin_h__ */

