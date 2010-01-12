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
** Invariants:
**
** Comments:
**
** $Id: PtrSet.h 2225 1970-01-02 13:48:23Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StGermain_Base_Container_PtrSet_h__
#define __StGermain_Base_Container_PtrSet_h__
	

	/* Textual name of this class */
	extern const Type PtrSet_Type;

	/* Virtual function types */
	
	/* Support structures */
	
	/** PtrSet class contents */
	#define __PtrSet \
		/* General info */ \
		__Set \
		\
		/* Virtual info */ \
		\
		/* PtrSet info ... */ \

	struct PtrSet { __PtrSet };
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	/* Create an instance with all parameters */
	PtrSet* PtrSet_New( 
		Dictionary*					dictionary );
	
	/* Creation implementation */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define PTRSET_DEFARGS \
                SET_DEFARGS

	#define PTRSET_PASSARGS \
                SET_PASSARGS

	PtrSet* _PtrSet_New(  PTRSET_DEFARGS  );
	
	
	/* Initialise an instance */
	void PtrSet_Init(
		PtrSet*						self,
		Dictionary*					dictionary );
	
	/* Initialisation implementation functions */
	void _PtrSet_Init(
		PtrSet*						self );
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/
	
	/* Stg_Class_Delete implementation */
	void _PtrSet_Delete(
		void*						ptrSet );
	
	/* Print implementation */
	void _PtrSet_Print(
		void*						ptrSet, 
		Stream*						stream );

	void* _PtrSet_Union( void* ptrSet, void* operand );

	void* _PtrSet_Intersection( void* ptrSet, void* operand );

	void* _PtrSet_Subtraction( void* ptrSet, void* operand );
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	int _PtrSet_CompareData( void* left, void* right );

	void _PtrSet_DeleteData( void* data );


#endif /* __StGermain_Base_Container_PtrSet_h__ */

