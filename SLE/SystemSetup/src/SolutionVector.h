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
** Role:
**	Holds the solution to a FEM problem
**
** Assumptions:
**
** Comments:
**
** $Id: SolutionVector.h 656 2006-10-18 06:45:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_SLE_SystemSetup_SolutionVector_h__
#define __StgFEM_SLE_SystemSetup_SolutionVector_h__
	
	/* Textual name of this class */
	extern const Type SolutionVector_Type;
	
	/* StiffnessMatrix information */
	#define __SolutionVector  \
		/* General info */ \
		__Stg_Component \
		\
		/* Virtual info */ \
		\
		/* StiffnessMatrix info */ \
		Stream*						debug; \
		Vector*						vector; \
		MPI_Comm					comm; \
		FeVariable*					feVariable; /** need to get # of global unconstrained dofs */\

	struct SolutionVector { __SolutionVector };
	
	/* Creation implementation / Virtual constructor */
	void* SolutionVector_DefaultNew( Name name );

	SolutionVector* SolutionVector_New(
		Name						name,
		MPI_Comm					comm,
		FeVariable*					feVariable );

	SolutionVector* SolutionVector_New_FromArray(
		Name						name,
		MPI_Comm					comm,
		FeVariable*					feVariable );

	void SolutionVector_Init(
		SolutionVector*					self,
		Name						name,
		MPI_Comm					comm,
		FeVariable*					feVariable );

	void SolutionVector_Init_FromArray(
		SolutionVector*					self,
		Name						name,
		MPI_Comm					comm,
		FeVariable*					feVariable );


	SolutionVector* _SolutionVector_New(
		SizeT						_sizeOfSelf,
		Type						type,
		Stg_Class_DeleteFunction*				_delete,
		Stg_Class_PrintFunction*				_print,
		Stg_Class_CopyFunction*				_copy, 
		Stg_Component_DefaultConstructorFunction*			_defaultConstructor,
		Stg_Component_ConstructFunction*			_construct,
		Stg_Component_BuildFunction*			_build,
		Stg_Component_InitialiseFunction*			_initialise,
		Stg_Component_ExecuteFunction*			_execute,
		Stg_Component_DestroyFunction*			_destroy,
		Name						name,
		Bool						initFlag,
		MPI_Comm					comm,
		FeVariable*					feVariable );
		
	void SolutionVector_LoadFromDict( void* solutionVector, Dictionary* subDict, Dictionary* dictionary, Stg_ObjectList* objList );

	/* Initialise implementation */
	void _SolutionVector_Init( 
		SolutionVector*					self,
		MPI_Comm					comm,
		FeVariable*					feVariable );
	
	
	/* Stg_Class_Delete a ElementType construst */
	void _SolutionVector_Delete( void* solutionVector );
	
	/* Print the contents of an ElementType construct */
	void _SolutionVector_Print( void* solutionVector, Stream* stream );
	
	/* Copy */
	#define SolutionVector_Copy( self ) \
		(ForceVector*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define SolutionVector_DeepCopy( self ) \
		(ForceVector*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
	
	void* _SolutionVector_Copy( void* solutionVector, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	void _SolutionVector_Build( void* solutionVector, void* data );
	
	void _SolutionVector_Construct( void* solutionVector, Stg_ComponentFactory* cf, void* data );
	
	/* Initialisation implementation */
	void _SolutionVector_Initialise( void* solutionVector, void* data );
	
	/* Execution implementation */
	void _SolutionVector_Execute( void* solutionVector, void* data );
	
	void _SolutionVector_Destroy( void* solutionVector, void* data );
	
	void SolutionVector_ApplyBCsToVariables( void* solutionVector, void* data );

	void SolutionVector_UpdateSolutionOntoNodes( void* solutionVector );

	/** Loads the current value at each dof of the feVariable related to this solution vector onto the vector itself */
	void SolutionVector_LoadCurrentFeVariableValuesOntoVector( void* solutionVector );
	
#endif /* __StgFEM_SLE_SystemSetup_SolutionVector_h__ */
