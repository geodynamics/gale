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
**	Abtstract class to manage the set up, building, initialisation etc of a System of
**	Linear Equations. Subclasses should define exactly what Matrices and Vectors
**	are required. The actual Solve() function also needs to be provided for each
**	subclass.
**
** Assumptions:
**
** Comments:
**	Note - up to 1 September 2004 (rev 1994), much of the functionality of this class was in
**	the SLE_Solver class. The 2 classes have been reorganised: This class is responsible
**	for storing and managing the matrices and vectors that make up a system, but uses
**	the SLE_Solver class to actually implement a solution mechanism for the given eqn.
**
** $Id: SystemLinearEquations.h 1230 2008-09-15 01:44:43Z DavidLee $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_SLE_SystemSetup_SystemLinearEquations_h__
#define __StgFEM_SLE_SystemSetup_SystemLinearEquations_h__

	/** Textual name of this class */
	extern const Type SystemLinearEquations_Type;
	
	/* virtual function interface */
	typedef void (SystemLinearEquations_LM_SetupFunction) ( void* sle, void* data );
	typedef void (SystemLinearEquations_MatrixSetupFunction) ( void* sle, void* data );
	typedef void (SystemLinearEquations_VectorSetupFunction) ( void* sle, void* data );
	typedef void (SystemLinearEquations_UpdateSolutionOntoNodesFunc) ( void* sle, void* data );
	typedef void (SystemLinearEquations_MG_SelectStiffMatsFunc) ( void* _sle, unsigned* nSMs, StiffnessMatrix*** sms );

	/* for solving non linear systems using Newton's method */
	//typedef int (SystemLinearEquations_BuildFFunc) ( void* nls, Vector* x, Vector* f, void* context );
	//typedef int (SystemLinearEquations_BuildJFunc) ( void* nls, Vector* x, Matrix** A, Matrix** B, void* matStruct, void* context );	
	//typedef void (SystemLinearEquations_SetFFunc) ( Vector** F, void* context );
	typedef int (SystemLinearEquations_BuildFFunc) ( SNES nls, Vec x, Vec f, void* context );
	typedef int (SystemLinearEquations_BuildJFunc) ( SNES nls, Vec x, Mat* A, Mat* B, MatStructure* matStruct, void* context );	
	typedef void (SystemLinearEquations_SetFFunc) ( Vec* F, void* context );
	typedef void (SystemLinearEquations_ConfigureNonlinearSolver) ( void* nls, void* data );

	//typedef void (SLE_FormFunctionFunc)( void *someSLE, Vector *X, Vector *F, void *ctx );
	typedef void (SLE_FormFunctionFunc)( void *someSLE, Vec X, Vec F, void *ctx );

	
	/*
	** SystemLinearEquations class contents.
	*/
	
	#define __SystemLinearEquations \
		/* General info */ \
		__Stg_Component \
		ExtensionManager*                                   extensionManager; \
		\
		/* Virtual info */ \
		SystemLinearEquations_LM_SetupFunction*             _LM_Setup;                 \
		SystemLinearEquations_MatrixSetupFunction*          _matrixSetup;              \
		SystemLinearEquations_VectorSetupFunction*          _vectorSetup;              \
		SystemLinearEquations_UpdateSolutionOntoNodesFunc*  _updateSolutionOntoNodes;  \
		SystemLinearEquations_MG_SelectStiffMatsFunc*       _mgSelectStiffMats;        \
		\
		/* SystemLinearEquations info */ \
		Stream*                                             debug;                     \
		Stream*                                             info;                      \
		Stream*                                             convergenceStream;         \
		Bool                                                makeConvergenceFile;       \
		MPI_Comm                                            comm;                      \
		StiffnessMatrixList*                                stiffnessMatrices;         \
		ForceVectorList*                                    forceVectors;              \
		SolutionVectorList*                                 solutionVectors;           \
		SLE_Solver*	                                    solver;                    \
		FiniteElementContext*                               context;                   \
		Name                                                executeEPName;             \
		EntryPoint*                                         executeEP;                 \
		Name                                                integrationSetupEPName;    \
		EntryPoint*                                         integrationSetupEP;        \
		EntryPoint_Register*                                entryPoint_Register;       \
		\
		Bool                                                bcRemoveQuery;             \
		\
		/* Non-linear info */ \
		Bool                                                isNonLinear;               \
		Stg_Component_ExecuteFunction*                      linearExecute;             \
		double                                              nonLinearTolerance;        \
		Iteration_Index                                     nonLinearMaxIterations;    \
		Iteration_Index                                     nonLinearIteration_I;      \
		Bool                                                killNonConvergent;         \
		Iteration_Index                                     nonLinearMinIterations;    \
		double			                            curResidual;               \
		double                                              curSolveTime;              \
		/* BEGIN LUKE'S FRICTIONAL BCS BIT */					       \
		char*						    nlSetupEPName;	       \
		EntryPoint*					    nlSetupEP;		       \
		char*						    nlEPName;	       	       \
		EntryPoint*					    nlEP;		       \
		char*						    postNlEPName;	       \
		EntryPoint*					    postNlEP;		       \
		char*						    nlConvergedEPName;	       \
		EntryPoint*					    nlConvergedEP;	       \
                Bool                                                nlFormJacobian; \
                Vec                                                 nlCurIterate; \
		/* END LUKE'S FRICTIONAL BCS BIT */					       \
		/* Multi-grid data. */ \
		Bool                                                mgEnabled;                 \
		Bool                                                mgUpdate; \
		unsigned                                            nMGHandles;                \
		unsigned*                                           mgHandles; /* one per MG'd 'StiffnessMatrix' */ \
		/* for solving non linear systems using Newton's method */\
		Name						    nonLinearSolutionType;     \
		SystemLinearEquations_BuildFFunc*		    _buildF;		       \
		SystemLinearEquations_BuildJFunc*		    _buildJ;		       \
		void*						    buildFContext;             \
		void*						    buildJContext;             \
		SNES				    		    nlSolver;		       \
		Bool						    linearSolveInitGuess;      \
		Vec    						    F;			       \
		Vec    						    X;		       	       \
		Mat    						    J;			       \
		Mat    						    P;			       \
		SystemLinearEquations_SetFFunc*		    	    _setFFunc;		       \
		SystemLinearEquations_ConfigureNonlinearSolver*	    _configureNLSolverFunc;    \
		SystemLinearEquations_SetFFunc*		    	    _updateOldFields;	       \
		Name						    optionsPrefix;	       \
		/* parameter and methods relevant to PICARD */ \
		Name                   picard_form_function_type;  \
		double                 alpha;                      \
		double                 rtol,ttol,xtol,abstol;      \
		Bool                   picard_monitor;             \
		SLE_FormFunctionFunc* _sleFormFunction;
		
		
	/** Abstract class to manage the set up, building, initialisation etc of a System of
	Linear Equations - see SystemLinearEquations.h */
	struct SystemLinearEquations { __SystemLinearEquations };

	/** Constructor */
	SystemLinearEquations* SystemLinearEquations_New(
		Name                                               name,
		SLE_Solver*                                        solver,
		void*				   		   nlSolver,
		FiniteElementContext*                              context,
		Bool                                               isNonLinear,
		double                                             nonLinearTolerance,
		Iteration_Index                                    nonLinearMaxIterations,
		Bool                                               killNonConvergent,		
		EntryPoint_Register*                               entryPoint_Register,
		MPI_Comm                                           comm );

	/** Creation implementation / Virtual constructor */
	SystemLinearEquations* _SystemLinearEquations_New( 
		SizeT                                              sizeOfSelf,
		Type                                               type,
		Stg_Class_DeleteFunction*                          _delete,
		Stg_Class_PrintFunction*                           _print,
		Stg_Class_CopyFunction*                            _copy, 
		Stg_Component_DefaultConstructorFunction*          _defaultConstructor,
		Stg_Component_ConstructFunction*                   _construct,
		Stg_Component_BuildFunction*                       _build,
		Stg_Component_InitialiseFunction*                  _initialise,
		Stg_Component_ExecuteFunction*                     _execute,
		Stg_Component_DestroyFunction*                     _destroy,
		SystemLinearEquations_LM_SetupFunction*            _LM_Setup,
		SystemLinearEquations_MatrixSetupFunction*         _matrixSetup,
		SystemLinearEquations_VectorSetupFunction*         _vectorSetup,
		SystemLinearEquations_UpdateSolutionOntoNodesFunc* _updateSolutionOntoNodes,  
		SystemLinearEquations_MG_SelectStiffMatsFunc*      _mgSelectStiffMats, 
		Name                                               name );

	void SystemLinearEquations_InitAll( 
		void*                                              sle, 
		SLE_Solver*                                        solver, 
		void*				   		   nlSolver,
		FiniteElementContext*                              context, 
		Bool                                               isNonLinear,
		double                                             nonLinearTolerance,
		Iteration_Index                                    nonLinearMaxIterations,
		Bool                                               killNonConvergent,
		EntryPoint_Register*                               entryPoint_Register,
		MPI_Comm                                           comm );

	/* Stg_Class_Delete() implementations */
	void _SystemLinearEquations_Delete( void* sle );

	/* Stg_Class_Print() implementation */
	void _SystemLinearEquations_Print( void* sle, Stream* stream );
	
	/* Copy implementation */
	#define SystemLinearEquations_Copy( self ) \
		(SystemLinearEquations*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define SystemLinearEquations_DeepCopy( self ) \
		(SystemLinearEquations*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
		
	void* _SystemLinearEquations_Copy( void* sle, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	/* +++ Virtual Functions +++ */
	void* _SystemLinearEquations_DefaultNew( Name name ) ;
	void _SystemLinearEquations_Construct( void* sle, Stg_ComponentFactory* cf, void* data );
	void _SystemLinearEquations_Build( void* sle, void* data );
	void _SystemLinearEquations_Initialise( void* sle, void* data );

	/** Stg_Component_Execute() implementation: Assembles the correct values into all matrices of the
	system, then executes the system's solver. */
	void _SystemLinearEquations_Execute( void* sle, void* data );
	
	void _SystemLinearEquations_Destroy( void* sle, void* data );

	void SystemLinearEquations_ExecuteSolver( void* sle, void* data ) ;

	void SystemLinearEquations_BC_Setup( void* sle, void* data );

	/* LM Setup */
	void SystemLinearEquations_LM_Setup( void* sle, void* data );
	void _SystemLinearEquations_LM_Setup( void* sle, void* data );

	void SystemLinearEquations_IntegrationSetup( void* sle, void* data );
	
	/* Matrix Setup */
	void SystemLinearEquations_MatrixSetup( void* sle, void* data );
	void _SystemLinearEquations_MatrixSetup( void* sle, void* data );

	/* Vector Setup */
	void SystemLinearEquations_VectorSetup( void* sle, void* data );
	void _SystemLinearEquations_VectorSetup( void* sle, void* data );

	/* +++ Public Functions / Macros +++ */
	
	/* AddStiffnessMatrix macro */
	#define SystemLinearEquations_AddStiffnessMatrix( sle, stiffnessMatrix ) \
		Stg_ObjectList_Append( ((sle)->stiffnessMatrices), stiffnessMatrix )
	Index _SystemLinearEquations_AddStiffnessMatrix( void* sle, StiffnessMatrix* stiffnessMatrix );

	/* GetStiffnessMatrix macro */
	#define SystemLinearEquations_GetStiffnessMatrix( sle, stiffnessMatrixName ) \
		((StiffnessMatrix*)Stg_ObjectList_Get( ((sle)->stiffnessMatrices), stiffnessMatrixName ))
	StiffnessMatrix* _SystemLinearEquations_GetStiffnessMatrix( void* sle, Name stiffnessMatrixName );

	#define SystemLinearEquations_GetStiffnessMatrixAt( sle, stiffnessMatrixIndex ) \
		((StiffnessMatrix*)Stg_ObjectList_At( ((sle)->stiffnessMatrices), stiffnessMatrixIndex ))

	/* AddForceVector macro */
	#define SystemLinearEquations_AddForceVector( sle, forceVector ) \
		Stg_ObjectList_Append( ((sle)->forceVectors), forceVector )
	Index _SystemLinearEquations_AddForceVector( void* sle, ForceVector* forceVector );

	/* GetForceVector macro */
	#define SystemLinearEquations_GetForceVector( sle, forceVectorName ) \
		((ForceVector*)Stg_ObjectList_Get( ((sle)->forceVectors), forceVectorName ))
	ForceVector* _SystemLinearEquations_GetForceVector( void* sle, Name forceVectorName );

	#define SystemLinearEquations_GetForceVectorAt( sle, forceVectorIndex ) \
		((ForceVector*)Stg_ObjectList_At( ((sle)->forceVectors), forceVectorIndex ))

	/* AddSolutionVector macro */
	#define SystemLinearEquations_AddSolutionVector( sle, solutionVector ) \
		Stg_ObjectList_Append( ((sle)->solutionVectors), solutionVector )
	Index _SystemLinearEquations_AddSolutionVector( void* sle, SolutionVector* solutionVector );

	/* GetSolutionVector macro */
	#define SystemLinearEquations_GetSolutionVector( sle, solutionVectorName ) \
		((SolutionVector*)Stg_ObjectList_Get( ((sle)->solutionVectors), solutionVectorName ))
	SolutionVector* _SystemLinearEquations_GetSolutionVector( void* sle, Name solutionVectorName );

	#define SystemLinearEquations_GetSolutionVectorAt( sle, solutionVectorIndex ) \
		((SolutionVector*)Stg_ObjectList_At( ((sle)->solutionVectors), solutionVectorIndex ))

	/* Update all solution vectors back onto mesh nodes */
	void SystemLinearEquations_UpdateSolutionOntoNodes( void* sle, void* data ) ;
	void _SystemLinearEquations_UpdateSolutionOntoNodes( void* sle, void* data );

	void SystemLinearEquations_ZeroAllVectors( void* sle, void* data ) ;

	/* Non-linear stuff */
	/* matrix free finite difference newton's method non linear solve */
	void SystemLinearEquations_NewtonMFFDExecute( void* sle, void* data );
	/* solitary waves model with hand rolled J */
	void SystemLinearEquations_NewtonInitialise( void* sle, void* data );
	void SystemLinearEquations_NewtonExecute( void* sle, void* data );
	void SystemLinearEquations_NewtonFinalise( void* sle, void* data );
	void SystemLinearEquations_NonLinearExecute( void* sle, void* data ) ;
	void SystemLinearEquations_AddNonLinearSetupEP( void* sle, const char* name,
                                                        EntryPoint_2VoidPtr_Cast func );
	void SystemLinearEquations_AddNonLinearEP( void* sle, const char* name,
                                                   EntryPoint_2VoidPtr_Cast func );
	void SystemLinearEquations_AddPostNonLinearEP( void* sle, const char* name,
                                                       EntryPoint_2VoidPtr_Cast func );
	void SystemLinearEquations_SetToNonLinear( void* sle ) ;
	void SystemLinearEquations_CheckIfNonLinear( void* sle ) ;
	
	/*
	** All the multi-grid virtual functions and their general implementations.
	*/
	
	void SystemLinearEquations_MG_Enable( void* _sle );
	
	void SystemLinearEquations_MG_SelectStiffMats( void* _sle, unsigned* nSMs, StiffnessMatrix*** sms );
	void _SystemLinearEquations_MG_SelectStiffMats( void* _sle, unsigned* nSMs, StiffnessMatrix*** sms );
	
#endif
