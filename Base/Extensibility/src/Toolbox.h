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
** $Id: Toolbox.h 3192 2005-08-25 01:45:42Z AlanLo $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Base_Extensibility_Toolbox_h__
#define __Base_Extensibility_Toolbox_h__
	
	/** The prototype for the Register function in a toolbox */
	typedef Index (Toolbox_RegisterFunction)   ( void* toolboxesManager );
	typedef Index (Toolbox_InitialiseFunction) ( void* toolboxesManager, int* argc, char*** argv );
	typedef Index (Toolbox_FinaliseFunction)   ( void* toolboxesManager );

	/* Textual name of this class */
	extern const Type Toolbox_Type;

	/* Toolboxs info */
	#define __Toolbox \
		/* General info */ \
		__Module \
		\
		/* Virtual info */ \
		\
		/* Toolbox info */ \
		Toolbox_RegisterFunction*    Register; \
		Toolbox_InitialiseFunction*  Initialise; \
		Toolbox_FinaliseFunction*    Finalise;
		
	struct Toolbox { __Toolbox };


	/* Create a new Toolbox */
	Toolbox* Toolbox_New( Name name, Stg_ObjectList* directories );
	Module* Toolbox_Factory( Name name, Stg_ObjectList* directories );
	
	/* Creation implementation / Virtual constructor */
	Toolbox* _Toolbox_New( 
		SizeT                        _sizeOfSelf,
		Type                         type,
		Stg_Class_DeleteFunction*    _delete,
		Stg_Class_PrintFunction*     _print,
		Stg_Class_CopyFunction*      _copy, 
		Name                         name,
		Stg_ObjectList*              directories );
	
	/* Initialisation implementation */
	void _Toolbox_Init( Toolbox* self );
	
	/* Delete implementation */
	void _Toolbox_Delete( void* toolbox );
	
	/* Print implementation */
	void _Toolbox_Print( void* toolbox, Stream* stream );

	/** Get the function pointer the to the toolbox's register function */
	Toolbox_RegisterFunction* Toolbox_GetRegisterFunc( void* toolbox );
	
	/** Get the function pointer the to the toolbox's register function */
	Toolbox_InitialiseFunction* Toolbox_GetInitialiseFunc( void* toolbox );
	
	/** Get the function pointer the to the toolbox's register function */
	Toolbox_FinaliseFunction* Toolbox_GetFinaliseFunc( void* toolbox );
	
#endif /* __Base_Extensibility_Toolbox_h__ */
