
#ifndef __Underworld_Utils_SPR_StrainRate_h__
#define __Underworld_Utils_SPR_StrainRate_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type SPR_StrainRate_Type;

	/** SPR_StrainRate class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */

	#define __SPR_StrainRate \
 	/* Macro defining parent goes here - This means you can cast this class as its parent */ \
 		__BaseRecoveryFeVar \
		/* Virtual functions go here */ \

	struct SPR_StrainRate { __SPR_StrainRate };

  /* public constructor for c code */
  SPR_StrainRate* SPR_StrainRate_New(
		Name							name,
		DomainContext*				context,
		void*                   feMesh,
 		void*                   geometryMesh,
		DofLayout*              dofLayout,                                                                                 
		void*                   bcs,
		void*                   ics,
		void*                   linkedDofInfo,
		void*                   templateFeVariable,    
		Index							fieldComponentCount,
		Dimension_Index			dim,
		Bool                    isCheckpointedAndReloaded, 
		Bool                    isReferenceSolution,
		Bool                    loadReferenceEachTimestep,
		MPI_Comm						communicator,
		FieldVariable_Register*	fV_Register,
		Variable_Register*		vr, 
		FeVariable*					rawField, 
		int							rawOrderOfInterpolation,
		Bool							coeffInterpolation );

 	void* _SPR_StrainRate_DefaultNew( Name name );

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define SPR_STRAINRATE_DEFARGS \
                BASERECOVERYFEVAR_DEFARGS

	#define SPR_STRAINRATE_PASSARGS \
                BASERECOVERYFEVAR_PASSARGS

	SPR_StrainRate* _SPR_StrainRate_New(  SPR_STRAINRATE_DEFARGS  );

	void _SPR_StrainRate_Execute(void* sprVar, void* data);

	void _SPR_StrainRate_Destroy(void* sprVar, void* data);

	void _SPR_StrainRate_Delete( void* sprVar );

	void _SPR_StrainRate_Print( void* sprVar, Stream* stream );

	void* _SPR_StrainRate_Copy( void* sprVar, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

	void _SPR_StrainRate_Init( SPR_StrainRate* self );

	void _SPR_StrainRate_AssignFromXML( void* sprVar, Stg_ComponentFactory* cf, void* data );

	void _SPR_StrainRate_Build( void* sprVar, void* data );

	void _SPR_StrainRate_Initialise( void* sprVar, void* data );

	int FeMesh_PopulateBoundaryNodesInfo( FeMesh* mesh, BoundaryNodesInfo* bninfo );

	void _SPR_StrainRate_AssembleSolveLocalPatchs( void* sprVar );

	void _SPR_StrainRate_AssemblePatch( SPR_StrainRate* self, int node_I, double** AMatrix, double** bVector);

	void _SPR_StrainRate_Boundaries( void* spr );

	void _SPR_StrainRate_GetCoeffAtNode( void* feVariable, Node_DomainIndex dNode_I, double* coeff );

 	void _SPR_StrainRate_GetValueInElement( void* feVariable, Element_Index lEl_I, double* xi, double* value );

	void _SPR_StrainRate_GetValueAtNode( void* feVariable, Node_DomainIndex dNode_I, double* value );

 	double _SPR_StrainRate_ApplyCoeff( SPR_StrainRate* self, double* coeff, double* coord, int order, double* pVec );

#endif

