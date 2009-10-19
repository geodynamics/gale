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
**	Class allowing the user to construct and solve a FEM problem
**
** Assumptions:
**	The context which inherits from this one, or extensions which add to it, will take care of deleting the base level
**	objects (e.g. meshes, swarms etc), as this context can't tell which FE variables share which objects.
**
** Comments:
**	There's an issue inside FE_Context_Build() see inside the code.
**
** $Id: Context.h 1212 2008-08-31 13:57:54Z LukeHodkinson $
*
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_SLE_SystemSetup_Context_h__
#define __StgFEM_SLE_SystemSetup_Context_h__
	
	/* Textual name of this class */
	extern const Type FiniteElementContext_Type;

	extern const Name FiniteElementContext_EP_CalcDt;
	
	#define __FiniteElementContext \
		/* General info */ \
		__DomainContext \
		\
		/* Virtual info */ \
		\
		/* FiniteElementContext info */ \
		SystemLinearEquationList*           slEquations;             \
		double                              prevTimestepDt;          \
		Bool                                limitTimeStepIncreaseRate; \
		double                              maxTimeStepIncreasePercentage; \
                double                              maxTimeStepSize;\
		EntryPoint*                         calcDtEP;                 \
		
	struct FiniteElementContext { __FiniteElementContext };
	
	/* Constructors ----------------------------------------------------------------------------------------------------------*/
	
	/** Constructor */
	void* FiniteElementContext_DefaultNew( Name name );

	FiniteElementContext* FiniteElementContext_New( 
		Name                                            name,
		double						start,
		double						stop,
		MPI_Comm					communicator,
		Dictionary*					dictionary );

	
	/** Creation implementation / Virtual constructor */
	FiniteElementContext* _FiniteElementContext_New( 
		SizeT						sizeOfSelf,
		Type						type,
		Stg_Class_DeleteFunction*			_delete,
		Stg_Class_PrintFunction*			_print,
		Stg_Class_CopyFunction*				_copy, 
		Stg_Component_DefaultConstructorFunction*	_defaultConstructor,
		Stg_Component_ConstructFunction* 		_construct,
		Stg_Component_BuildFunction*			_build,
		Stg_Component_InitialiseFunction*		_initialise,
		Stg_Component_ExecuteFunction* 			_execute,
		Stg_Component_DestroyFunction*			_destroy,
		Name 						name,
		Bool						initFlag,
		AbstractContext_SetDt*				_setDt,
		double						start,
		double						stop,
		MPI_Comm					communicator,
		Dictionary*					dictionary );
	
	/** Initialisation implementation */
	void  _FiniteElementContext_Init( FiniteElementContext* self );

	/* Virtual Functions -----------------------------------------------------------------------------------------------------*/
	
	/* Stg_Class_Delete implementation */
	void _FiniteElementContext_Delete( void* context );
	
	/* Print implementation */
	void _FiniteElementContext_Print( void* context, Stream* stream );

	/* Set the dt */
	void _FiniteElementContext_CalcNewDt( void* context, double dt );
	
	/* Construct EntryPoint hook */
	void _FiniteElementContext_Construct( void* context, Stg_ComponentFactory* cf, void* data );
	
	/* Build EntryPoint hook */
	void _FiniteElementContext_Build( void* context );
	
	/* Initialise EntryPoint hook */
	void _FiniteElementContext_Initialise( void* context );
	
	/* Solve EntryPoint hook */
	void _FiniteElementContext_Solve( void* context );

	/* PostSolve EntryPoint hook: for things such as updating the Dt */
	void _FiniteElementContext_PostSolve( void* context ) ;
	
	/* Dt related functions */
	/** Function to assign the dt mediated by the abstract context to use for time integration
		in the current step - usually just the last timestep's dt */
	void _FiniteElementContext_SetDt( void* context, double dt );
	/** Function to pass the last-calculated dt back to the context. If loading from checkpoint,
		will load the dt from timeInfo file */
	double _FiniteElementContext_GetDt( void* context );

	/** Function to calculate the new dt based on the solution just obtained, which
	  * will be used next timestep. Calls as entry point of a corresponding name,
	  * which each SLE can add a hook to based on its own criterion */
	double FiniteElementContext_CalcNewDt( void* context ) ;

	/* Public functions ------------------------------------------------------------------------------------------------------*/
	
	void FiniteElementContext_AddSLE_Func( void* context, void* sle );

	#define FiniteElementContext_AddSLE_Macro( self, sle ) \
		Stg_ObjectList_Append( (self)->slEquations, sle )

	#ifdef MACRO_AS_FUNC
		#define FiniteElementContext_AddSLE FiniteElementContext_AddSLE_Func
	#else
		#define FiniteElementContext_AddSLE FiniteElementContext_AddSLE_Macro
	#endif
	#define AddSLE FiniteElementContext_AddSLE
	
	SystemLinearEquations* FiniteElementContext_GetSLE_Func( void* context, Name sleName );
	#define FiniteElementContext_GetSLE_Macro( self, sleName ) \
		((SystemLinearEquations*)Stg_ObjectList_Get( (self)->slEquations, sleName ))
	#ifdef MACRO_AS_FUNC
		#define FiniteElementContext_GetSLE FiniteElementContext_GetSLE_Func
	#else
		#define FiniteElementContext_GetSLE FiniteElementContext_GetSLE_Macro
	#endif
	#define GetSLE FiniteElementContext_GetSLE

	/* Private functions -----------------------------------------------------------------------------------------------------*/

	/** Saves all FeVariables known about by the context to file */
	void _FiniteElementContext_SaveFeVariables( void* context );

	/** Saves all Swarms known about by the context to file (necessary to put this here since 
	gauss integration swarms may need to be saved and reloaded) */
	void _FiniteElementContext_SaveSwarms( void* context );

   void _FiniteElementContext_SaveMesh( void* context );
   
   #ifdef WRITE_HDF5
	void _FiniteElementContext_DumpMeshHDF5( void* context, FeMesh* mesh );
   #else
	void _FiniteElementContext_DumpMeshAscii( void* context, FeMesh* mesh );
   #endif

#endif /* __StgFEM_SLE_SystemSetup_Context_h__ */
