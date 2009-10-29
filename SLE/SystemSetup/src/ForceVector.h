/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**	Computational Infrastructure for Geodynamics - http://www.geodynamics.org
**
** Contributors:
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
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
**	Finite Element force vector - holds references to the data relevant to a RHS "force vector", and provides
**		an entry point that the user fills in to build it.
**
** Assumptions:
**
** Comments:
**
** $Id: ForceVector.h 1210 2008-08-25 01:17:12Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_SLE_SystemSetup_ForceVector_h__
#define __StgFEM_SLE_SystemSetup_ForceVector_h__
	
	
	/* Textual name of this class */
	extern const Type ForceVector_Type;
	
	/* StiffnessMatrix information */
	#define __ForceVector  \
		/* General info */ \
		__SolutionVector \
		\
		/* Virtual info */ \
		\
		/* StiffnessMatrix info */ \
		Index				localSize;  \
		Dimension_Index			dim;  \
		EntryPoint_Register*		entryPoint_Register;  \
		FeEntryPoint*			assembleForceVector;  \
		Name				_assembleForceVectorEPName;  \
		Stg_ObjectList*			forceTermList;  \
		Stg_Component*			applicationDepExtraInfo; /**< Default is NULL: passed to elForceVec during assembly */\
		Assembler*			bcAsm;  \
		IArray*				inc;  \
                int nModifyCBs;                       \
                Callback* modifyCBs;
	
	struct ForceVector { __ForceVector };
	
	/* Creation implementation / Virtual constructor */
	
	ForceVector* ForceVector_New(
		Name                                      name,
		FeVariable*                               feVariable,
		Dimension_Index                           dim,
		void*                                     entryPoint_Register,
		MPI_Comm                                  comm );

	ForceVector* _ForceVector_New( 
		SizeT                                     _sizeOfSelf,
		Type                                      type,
		Stg_Class_DeleteFunction*                 _delete,
		Stg_Class_PrintFunction*                  _print,
		Stg_Class_CopyFunction*                   _copy, 
		Stg_Component_DefaultConstructorFunction* _defaultConstructor,
		Stg_Component_ConstructFunction*          _construct,
		Stg_Component_BuildFunction*              _build,
		Stg_Component_InitialiseFunction*         _initialise,
		Stg_Component_ExecuteFunction*            _execute,
		Stg_Component_DestroyFunction*            _destroy,
		Name                                      name,
		Bool                                      initFlag,
		FeVariable*                               feVariable,
		Dimension_Index                           dim,
		void*                                     entryPoint_Register,
		MPI_Comm                                  comm );
	
	void _ForceVector_Init(
		void*                                     forceVector,
		Dimension_Index                           dim,
		EntryPoint_Register*                      entryPoint_Register );
	
	/* 'Stg_Class' Virtual Functions */
	void _ForceVector_Delete( void* forceVector );
	void _ForceVector_Print( void* forceVector, Stream* stream );
	#define ForceVector_Copy( self ) \
		(ForceVector*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define ForceVector_DeepCopy( self ) \
		(ForceVector*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _ForceVector_Copy( void* forceVector, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	/* 'Stg_Component' Virtual Functions */
	void* _ForceVector_DefaultNew( Name name );
	void _ForceVector_AssignFromXML( void* forceVector, Stg_ComponentFactory* cf, void* data );
	void _ForceVector_Build( void* forceVector, void* data );
	void _ForceVector_Initialise( void* forceVector, void* data );
	void _ForceVector_Execute( void* forceVector, void* data );
	void _ForceVector_Destroy( void* forceVector, void* data );
	
	/** Interface to assemble this Force Vector. Calls an entry point, meaning the user can specify if, and then how,
	it should be assembled. */
	void ForceVector_Assemble( void* forceVector );
	
	/** Prints the contents of a single element's force vector */
	void ForceVector_PrintElementForceVector(
		ForceVector* self,
		Element_LocalIndex element_lI,
		Dof_EquationNumber** elementLM,
		double* elForceVecToAdd );
	
	void ForceVector_GlobalAssembly_General( void* forceVector ) ;
	void ForceVector_AssembleElement( void* forceVector, Element_LocalIndex element_lI, double* elForceVecToAdd ) ;
	void ForceVector_AddForceTerm( void* forceVector, ForceTerm* forceTerm ) ;

void ForceVector_AddModifyCallback( ForceVector* self, void* callback, void* object );

#endif /* __StgFEM_SLE_SystemSetup_ForceVector_h__ */
