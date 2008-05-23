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
**	Brings together and manages the life cycle of all the components required by the
**	Finite Element Method about a variable to be solved for.
**
** Assumptions:
**
** Comments:
**
** $Id $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_Discretisation_FeVariable_h__
#define __StgFEM_Discretisation_FeVariable_h__
	
	/** Textual name of this class */
	extern const Type FeVariable_Type;

	typedef void (FeVariable_InterpolateWithinElementFunction) (void* fieldVariable, Element_DomainIndex dElement_I, double* xi, double* value );
	typedef void (FeVariable_GetValueAtNodeFunction) (void* feVariable, Node_DomainIndex dNode_I, double* value );
	typedef void (FeVariable_SyncShadowValuesFunc)( void* feVariable );
	
	/* Function prototypes for import / export */
	typedef void (FeVariable_ReadNodalValuesFromFile_Function) (void* feVariable, const char* prefixStr, unsigned int timeStep );
	typedef void (FeVariable_SaveNodalValuesToFile_Function) (void* feVariable, const char* prefixStr, unsigned int timeStep );

	/** Struct containing info of how to read from / export to a certain file format */
	/* We expect several of these guys to be registered in the Discretisation Init phase, and possibly by later
	 * plugins */
	typedef struct FeVariable_ImportExportInfo {
		FeVariable_ReadNodalValuesFromFile_Function*  readNodalValuesFromFile;   
		FeVariable_SaveNodalValuesToFile_Function*    saveNodalValuesToFile;   
	} FeVariable_ImportExportInfo;

	void FeVariable_ImportExportInfo_Delete( void* ptr );
	void FeVariable_ImportExportInfo_Print ( void* ptr, struct Stream* stream );
	void* FeVariable_ImportExportInfo_Copy( 
		void*					ptr, 
		void*					dest,
		Bool					deep,
		Name					nameExt, 
		struct PtrMap*				ptrMap );

	/* A global list of import/export info objects - can be added to later by plugins. Needs to be initialised in
	 * FeDiscretisation_Init() */
	/* TODO: maybe should move this to the StgFEM_Context later???  */
	extern Stg_ObjectList*   FeVariable_FileFormatImportExportList;

	extern const char*       StgFEM_Native_ImportExportType;
	
	/** FeVariable class contents */
	#define __FeVariable \
		/* General info */ \
		__FieldVariable \
		\
		/* Virtual info */ \
		FeVariable_InterpolateWithinElementFunction*      _interpolateWithinElement; \
		FeVariable_GetValueAtNodeFunction*                _getValueAtNode;          \
		FeVariable_SyncShadowValuesFunc*		_syncShadowValues; \
		\
		/* FeVariable info */ \
		\
		Stream*                                           debug; \
		/** Mesh that this variable is discretised over */ \
		FeMesh*						feMesh; \
		/** The mesh that this variable can create the determinant of the jacobian from */ \
		FeMesh*						geometryMesh; \
		/** DofLayout for this variable: relates each mesh node to the Variable's */ \
		DofLayout*                                        dofLayout; \
		/** Boundary conditions applied to this variable - Compulsory, so the eq num table can be worked out*/ \
		VariableCondition*                                bcs; \
		DynamicVC*					dynamicBCs[2]; /* Temporary hack */	\
		Bool						removeBCs;	\
		/** Boundary conditions applied to this variable - Optional, may be NULL */ \
		VariableCondition*                                ics; \
		/** Info on which dofs are linked together: optional, may be NULL */ \
		LinkedDofInfo*                                    linkedDofInfo; \
		/** Equation number array: maps where each dof of this Variable goes in any matrices
		based off it. */ \
		FeEquationNumber*                                 eqNum;  \
		/** Records whether the user has sync'd shadow values yet. */ \
		Bool                                              shadowValuesSynchronised;  \
		/** A "template" feVariable this one is based on - ie this one's mesh and BCs is based off that one */ \
		FeVariable*                                       templateFeVariable; \
		/** A type recording what import/export system for loading and saving should be used */ \
		char*                                             importFormatType; \
		char*                                             exportFormatType; \
		char*                                             customInputPath; \
		char*                                             customOutputPath; \
		/** Records whether this FeVariable is a reference solution - and should be loaded from a file, regardless
		 * of checkpointing status */ \
		Bool                                              isReferenceSolution; \
		/* if self->isReferenceSolution is true, this param determines if it's loaded once at the start, or every
		 * timestep. */ \
		Bool                                              loadReferenceEachTimestep; \
									\
		IArray* inc;


	/* Brings together and manages the life cycle of all the components required by the
	Finite Element Method about a variable to be solved for - see FeVariable.h */
	struct FeVariable { __FeVariable };
	
	/* --- Contstructors / Destructors --- */
	
	/** Create a new FeVariable and initialises it. The default one - no template. */
	
	void* FeVariable_DefaultNew( Name name );
	
	FeVariable* FeVariable_New(
		Name                                            name,
		void*                                           feMesh,
		void*                                           geometryMesh,
		DofLayout*                                      dofLayout, 
		void*                                          	bcs,
		void*                                          	ics,
		void*                                           linkedDofInfo,
		Dimension_Index                                 dim,
		Bool                                            isCheckpointedAndReloaded,
		const char* const                               importFormatType,
		const char* const                               exportFormatType,
		const char* const                               customInputPath,
		const char* const                               customOutputPath,
		Bool                                            referenceSoulution,
		Bool                                            loadReferenceEachTimestep,
		FieldVariable_Register*                         fV_Register );

	/** Create a new FeVariable and initialises it. Mesh, bcs and eqNum table is based off a template one.
	User is required to provide new ICs (we figured this would be the case in the vast majority of
	implementations. */
	FeVariable* FeVariable_New_FromTemplate(
		Name                                            name,
		void*                                          	templateFeVariable,
		DofLayout*                                      dofLayout, 
		void*                                          	ics,
		const char* const                               importFormatType,
		const char* const                               exportFormatType,
		const char* const                               customInputPath,
		const char* const                               customOutputPath,
		Bool                                            isReferenceSolution,
		Bool                                            loadReferenceEachTimestep,
		FieldVariable_Register*                         fV_Register );
	
	/** Create a new FeVariable and initialises it. User chooses whether to pass a template or not. */
	FeVariable* FeVariable_New_Full(
		Name                                            name,
		void*                                           feMesh,
		void*                                           geometryMesh,
		DofLayout*                                      dofLayout, 
		void*                                          	bcs,
		void*                                          	ics,
		void*                                           linkedDofInfo,
		void*                                          	templateFeVariable,
		Index                                           fieldComponentCount,
		Dimension_Index                                 dim,
		Bool                                            isCheckpointedAndReloaded,
		const char* const                               importFormatType,
		const char* const                               exportFormatType,
		const char* const                               customInputPath,
		const char* const                               customOutputPath,
		Bool                                            referenceSoulution,
		Bool                                            loadReferenceEachTimestep,
		MPI_Comm                                        communicator,
		FieldVariable_Register*                         fV_Register );
	
	/* Creation implementation / Virtual constructor */
	FeVariable* _FeVariable_New(
 		SizeT                                           _sizeOfSelf,
		Type                                            type,
		Stg_Class_DeleteFunction*                       _delete,
		Stg_Class_PrintFunction*                        _print,
		Stg_Class_CopyFunction*                         _copy, 
		Stg_Component_DefaultConstructorFunction*       _defaultConstructor,
		Stg_Component_ConstructFunction*                _construct,
		Stg_Component_BuildFunction*                    _build,
		Stg_Component_InitialiseFunction*               _initialise,
		Stg_Component_ExecuteFunction*                  _execute,
		Stg_Component_DestroyFunction*                  _destroy,
		Name                                            name,
		Bool                                            initFlag,
		FieldVariable_InterpolateValueAtFunction*       _interpolateValueAt,
		FieldVariable_GetValueFunction*                 _getMinGlobalFieldMagnitude,
		FieldVariable_GetValueFunction*                 _getMaxGlobalFieldMagnitude,
		FieldVariable_GetCoordFunction*                 _getMinAndMaxLocalCoords,
		FieldVariable_GetCoordFunction*                 _getMinAndMaxGlobalCoords,		
		FeVariable_InterpolateWithinElementFunction*    _interpolateWithinElement,	
		FeVariable_GetValueAtNodeFunction*              _getValueAtNode,
		FeVariable_SyncShadowValuesFunc*		_syncShadowValues, 
		void*                                           feMesh,
		void*                                           geometryMesh,
		DofLayout*                                      dofLayout, 
		void*                                          	bcs,
		void*                                          	ics,
		void*                                           linkedDofInfo,
		void*                                          	templateFeVariable,
		Index                                           fieldComponentCount,
		Dimension_Index                                 dim,
		Bool                                            isCheckpointedAndReloaded,
		const char* const                               importFormatType,
		const char* const                               exportFormatType,
		const char* const                               customInputPath,
		const char* const                               customOutputPath,
		Bool                                            referenceSoulution,
		Bool                                            loadReferenceEachTimestep,
		MPI_Comm                                        communicator,
		FieldVariable_Register*                         fV_Register );
	
	/* Initialise implementation */
	void _FeVariable_Init( 
		FeVariable*                                     self,
		void*                                           feMesh,
		void*                                           geometryMesh,
		DofLayout*                                      dofLayout, 
		void*                                           bcs,
		void*                                           ics,
		void*                                           linkedDofInfo,
		void*                                           templateFeVariable,
		const char* const                               importFormatType,
		const char* const                               exportFormatType,
		const char* const                               customInputPath,
		const char* const                               customOutputPath,
		Bool                                            referenceSoulution,
		Bool                                            loadReferenceEachTimestep );
	
	/** Stg_Class_Delete a FeVariable construst */
	void _FeVariable_Delete( void* variable );
	
	/* --- Virtual Function Implementations --- */
	#define FeVariable_InterpolateWithinElement( feVariable, dElement_I, xi, value ) \
		( ((FeVariable*) feVariable)->_interpolateWithinElement( feVariable, dElement_I, xi, value ) )

	#define FeVariable_GetValueAtNode( feVariable, dNode_I, value ) \
		( ((FeVariable*) feVariable)->_getValueAtNode( feVariable, dNode_I, value ) )

	#define FeVariable_SyncShadowValues( feVariable ) \
		( ((FeVariable*) feVariable)->_syncShadowValues( feVariable ) )

	/** Print the contents of an FeVariable construct */
	void _FeVariable_Print( void* variable, Stream* stream );
	
	/** Copy */
	#define FeVariable_Copy( self ) \
		(FeVariable*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
	#define FeVariable_DeepCopy( self ) \
		(FeVariable*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	
	void* _FeVariable_Copy( void* feVariable, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	/** Stg_Component_Build() implementation */
	void _FeVariable_Build( void* variable, void* data );
	
	/** Stg_Component_Construct() implementation */
	void _FeVariable_Construct( void* variable, Stg_ComponentFactory* cf, void* data );
	
	/** Stg_Component_Initialise() implementation */
	void _FeVariable_Initialise( void* variable, void* data );
	
	/** Stg_Component_Execute() implementation */
	void _FeVariable_Execute( void* variable, void* data );
	
	/** Stg_Component_Destroy() implementation */
	void _FeVariable_Destroy( void* variable, void* data );
	
	/** Apply BCs for this variable */
	void FeVariable_ApplyBCs( void* variable, void* data );
	Bool FeVariable_IsBC( void* variable, int node, int dof );
	
	/** Interpolate the value of the FE variable at a particular coord **/
	InterpolationResult _FeVariable_InterpolateValueAt( void* variable, double* coord, double* value );
	
	/* Implementations of the min and max val functions */
	double _FeVariable_GetMinGlobalFieldMagnitude( void* feVariable );

	double _FeVariable_GetMaxGlobalFieldMagnitude( void* feVariable );
	
	/* Implementations of the coord-getting functions */
	void _FeVariable_GetMinAndMaxLocalCoords( void* feVariable, double* min, double* max ) ;

	void _FeVariable_GetMinAndMaxGlobalCoords( void* feVariable, double* min, double* max ) ;

	/** Prints out the value at each DOF for this FeVariable */
	void FeVariable_PrintLocalDiscreteValues( void* variable, Stream* stream );

	void _FeVariable_GetValueAtNode( void* feVariable, Node_DomainIndex dNode_I, double* value ) ;


	/* --- Public Functions --- */

	/** Finds the value of the field at the node and broadcasts it to the rest of the processors 
	 *  It calls MPI_Allreduce - so this function must be called by all processors for it to work */
	void FeVariable_GetValueAtNodeGlobal( void* feVariable, Node_GlobalIndex gNode_I, double* value ) ;

	/** Finds the coordinate of the node and broadcasts it to the rest of the processors
	 *  It calls MPI_Allreduce - so this function must be called by all processors for it to work */
	void FeVariable_GetCoordAtNodeGlobal( void* feVariable, Node_GlobalIndex gNode_I, double* coord ) ;

	/** Zeros the value of the field at every nodal position */
	void FeVariable_ZeroField( void* feVariable ) ;

	/** Calculates the domain element & element local coord that a particular global coord lives in.
		Same return status conventions as for the InterpolateValueAt function. */
	InterpolationResult FeVariable_GetElementLocalCoordAtGlobalCoord( void* feVariable, double* globalCoord,
		double* elLocalCoord, Element_DomainIndex* elementCoordInPtr );

	/** Updates a single component of the value at a certain node */
	#define FeVariable_SetComponentAtNode( feVariable, dNode_I, dof_I, componentVal ) \
		DofLayout_SetValueDouble( (feVariable)->dofLayout, dNode_I, dof_I, componentVal );
		
	/** Updates all the components at a given node that user passes in with the componentValues array */
	void FeVariable_SetValueAtNode( void* feVariable, Node_DomainIndex dNode_I, double* componentValues );
	
	/** Gets the value of the FeVariable at a given node.
		if a scalar, just returns the value.
		if a vector, gives the magnitude.
	Probably should split these into separate functions */	
	double FeVariable_GetScalarAtNode( void* feVariable, Node_LocalIndex lNode_I );	

	/** Prints the values at the nodes in a manner that's easier to interpret if the geometry
	used is a 2D box. */
	void FeVariable_PrintLocalDiscreteValues_2dBox( void* variable, Stream* stream );

	/** Saves the current mesh coordinates, and value of each dof in the feVariable, to file */
	void FeVariable_SaveToFile( void* feVariable, const char* prefixStr, unsigned int timeStep );
	void FeVariable_SaveNodalValuesToFile_StgFEM_Native( void* feVariable, const char* prefixStr, unsigned int timeStep );

	/** Reads in everything to initialise a built FeVariable from a file */
	void FeVariable_ReadFromFile( void* feVariable, const char* prefixStr, unsigned int timeStep );
	void FeVariable_ReadNodalValuesFromFile_StgFEM_Native( void* feVariable, const char* prefixStr, unsigned int timeStep );

	/** Evaluates Spatial Derivatives using shape functions */
	Bool FeVariable_InterpolateDerivativesAt( void* variable, double* globalCoord, double* value ) ;
	
	void FeVariable_InterpolateDerivativesToElLocalCoord( void* _feVariable, Element_LocalIndex lElement_I, double* elLocalCoord, double* value ) ;
	
	void FeVariable_InterpolateDerivatives_WithGNx( void* _feVariable, Element_LocalIndex lElement_I, double** GNx, double* value ) ;
	
	void FeVariable_GetMinimumSeparation( void* feVariable, double* minSeparationPtr, double minSeparationEachDim[3] );

	/** Synchronises each processor's shadow dof values to be the same as the values on their "home" processors.
	 * Collective. */
	void _FeVariable_SyncShadowValues( void* feVariable );

	/** Perhaps should be moved into feVariable interface? */
	void FeVariable_PrintDomainDiscreteValues( void* feVariable, Stream* stream );
	void FeVariable_PrintCoordsAndValues( void* _feVariable, Stream* stream ) ;

	/** Use this function when you want the value of the field but the local coordinates you are using may not be appropriate for the mesh of this FeVariable */
	void FeVariable_InterpolateFromMeshLocalCoord( void* feVariable, FeMesh* mesh, Element_DomainIndex dElement_I, double* localCoord, double* value ) ;
	
	#define FeVariable_IntegrateElement( feVariable, swarm, dElement_I ) \
		FeVariable_IntegrateElement_AxisIndependent( \
				feVariable, swarm, dElement_I, ((FeVariable*) feVariable)->dim, \
				I_AXIS, J_AXIS, K_AXIS ) 

	double FeVariable_IntegrateElement_AxisIndependent( 
			void* feVariable, void* _swarm,
			Element_DomainIndex dElement_I, Dimension_Index dim, 
			Axis axis0, Axis axis1, Axis axis2 ) ;
	double FeVariable_Integrate( void* feVariable, void* _swarm ) ;

	/** Functions assumes IJK Topology Elements */
	double FeVariable_AverageTopLayer( void* feVariable, void* swarm, Axis layerAxis ) ;
	double FeVariable_AverageBottomLayer( void* feVariable, void* swarm, Axis layerAxis ) ;
	double FeVariable_AverageLayer( void* feVariable, void* swarm, Axis layerAxis, Index layerIndex ) ;

	#define FeVariable_IntegrateLayer( feVariable, swarm, layerAxis, layerIndex ) \
		FeVariable_IntegrateLayer_AxisIndependent( \
				feVariable, swarm, layerAxis, layerIndex, ((FeVariable*) feVariable)->dim,\
				I_AXIS, J_AXIS, K_AXIS ) 

	double FeVariable_IntegrateLayer_AxisIndependent( 
		void* feVariable, void* _swarm,
		Axis layerAxis, Index layerIndex, Dimension_Index dim, 
		Axis axis0, Axis axis1, Axis axis2 ) ;
	double FeVariable_AveragePlane( void* feVariable, Axis planeAxis, double planeHeight ) ;
	double FeVariable_IntegratePlane( void* feVariable, Axis planeAxis, double planeHeight ) ;

	/* --- Private Functions --- */

	/** Evaluates the value at a point within a given element, based on the current values of the nodes in that element */
	void _FeVariable_InterpolateNodeValuesToElLocalCoord( void* self, Element_DomainIndex element_lI, double* elLocalCoord, double* value );

	/** Utility function:- saves duplication of the print local and print domain functions */
	void _FeVariable_PrintLocalOrDomainValues( void* variable, Index localOrDomainCount, Stream* stream );


#endif /* __StgFEM_Discretisation_FeVariable_h__ */
