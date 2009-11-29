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
** Invariants:
**
** Comments:
**
** $Id: MGOpGenerator.h 675 2006-12-22 01:43:03Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Experimental_Solvers_MGOpGenerator_h__
#define __Experimental_Solvers_MGOpGenerator_h__

	/** Textual name of this class */
	extern const Type MGOpGenerator_Type;

	/** Virtual function types */
	typedef void (MGOpGenerator_SetNumLevelsFunc)( void* mgOpGenerator, unsigned nLevels );
	typedef Bool (MGOpGenerator_HasExpiredFunc)( void* mgOpGenerator );
	//typedef void (MGOpGenerator_GenerateFunc)( void* mgOpGenerator, Matrix*** pOps, Matrix*** rOps );
	typedef void (MGOpGenerator_GenerateFunc)( void* mgOpGenerator, Mat** pOps, Mat** rOps );

	/** MGOpGenerator class contents */
	#define __MGOpGenerator							\
		/* General info */						\
		__Stg_Component							\
										\
		/* Virtual info */						\
		MGOpGenerator_SetNumLevelsFunc*		setNumLevelsFunc;	\
		MGOpGenerator_HasExpiredFunc*		hasExpiredFunc;		\
		MGOpGenerator_GenerateFunc*		generateFunc;		\
										\
		/* MGOpGenerator info */					\
		MGSolver_PETScData*	solver;					\
		unsigned		nLevels;

	struct MGOpGenerator { __MGOpGenerator };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/





	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define MGOPGENERATOR_DEFARGS \
                STG_COMPONENT_DEFARGS, \
                MGOpGenerator_SetNumLevelsFunc*  setNumLevelsFunc, \
                MGOpGenerator_HasExpiredFunc*      hasExpiredFunc, \
                MGOpGenerator_GenerateFunc*          generateFunc

	#define MGOPGENERATOR_PASSARGS \
                STG_COMPONENT_PASSARGS, \
	        setNumLevelsFunc, \
	        hasExpiredFunc,   \
	        generateFunc    

	MGOpGenerator* _MGOpGenerator_New(  MGOPGENERATOR_DEFARGS  );
	void _MGOpGenerator_Init( MGOpGenerator* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _MGOpGenerator_Delete( void* mgOpGenerator );
	void _MGOpGenerator_Print( void* mgOpGenerator, Stream* stream );
	void _MGOpGenerator_AssignFromXML( void* mgOpGenerator, Stg_ComponentFactory* cf, void* data );
	void _MGOpGenerator_Build( void* mgOpGenerator, void* data );
	void _MGOpGenerator_Initialise( void* mgOpGenerator, void* data );
	void _MGOpGenerator_Execute( void* mgOpGenerator, void* data );
	void _MGOpGenerator_Destroy( void* mgOpGenerator, void* data );

	void _MGOpGenerator_SetNumLevels( void* mgOpGenerator, unsigned nLevels );

	#define MGOpGenerator_SetNumLevels( self, nLevels )		\
		VirtualCall( self, setNumLevelsFunc, self, nLevels )

	#define MGOpGenerator_HasExpired( self )			\
		VirtualCall( self, hasExpiredFunc, self )

	#define MGOpGenerator_Generate( self, pOps, rOps )		\
		VirtualCall( self, generateFunc, self, pOps, rOps )

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void MGOpGenerator_SetMatrixSolver( void* mgOpGenerator, void* solver );
	unsigned MGOpGenerator_GetNumLevels( void* mgOpGenerator );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

#endif /* __Experimental_Solvers_MGOpGenerator_h__ */

