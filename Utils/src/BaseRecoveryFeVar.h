#ifndef __Underworld_Utils_BaseRecoveryFeVar_h__
#define __Underworld_Utils_BaseRecoveryFeVar_h__

   typedef void (GetCoeffAtNode)( void* feVariable, Node_DomainIndex dNode_I, double* coeff );

   typedef void (Build_pVec)(  double* globalCoord, double* p_pVec );

   extern const Type BaseRecoveryFeVar_Type;
   extern Stg_ObjectList* repRequiredRawFields_Reg;

   #define __BaseRecoveryFeVar \
      /* parent class */ \
      __FeVariable \
      /* virtual methods */ \
      GetCoeffAtNode*		_getCoeffAtNode; \
      Build_pVec*				_makePoly; \
      /* member attributes */ \
      FeVariable*				rawField; \
      Variable*            dataVariable; \
      BoundaryNodesInfo*	bninfo; \
      Variable_Register*	vr;\
      double					*pVec; \
      double					orderOfInterpolation;

   struct BaseRecoveryFeVar { __BaseRecoveryFeVar };

	void* _BaseRecoveryFeVar_DefaultNew( Name name );

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define BASERECOVERYFEVAR_DEFARGS \
                FEVARIABLE_DEFARGS

	#define BASERECOVERYFEVAR_PASSARGS \
                FEVARIABLE_PASSARGS

	BaseRecoveryFeVar* _BaseRecoveryFeVar_New(  BASERECOVERYFEVAR_DEFARGS  );

	/* The virtual 5 functions */
	void _BaseRecoveryFeVar_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data );

	void _BaseRecoveryFeVar_Build( void* _self, void* data );

	void _BaseRecoveryFeVar_Initialise( void* _self, void* data );

	void _BaseRecoveryFeVar_Execute( void* _self, void* data );

	void _BaseRecoveryFeVar_Delete(void* _self);

	void _BaseRecoveryFeVar_Init(
      BaseRecoveryFeVar*	self, 
      Variable_Register*	vr, 
      FeVariable*				rawField, 
      int						rawOrderOfInterpolation,
      Bool						coeffInterpolation );

void _BaseRecoveryFeVar_Destroy(void* _self, void* data);

	void _BaseRecoveryFeVar_Boundaries( void* _self );

	void _BaseRecoveryFeVar_GetCoeffAtNode( void* feVariable, Node_DomainIndex dNode_I, double* coeff );

	void _BaseRecoveryFeVar_GetValueAtNode( void* feVariable, Node_DomainIndex dNode_I, double* value );

	void _BaseRecoveryFeVar_GetValueInElementWithCoeffInterpolation( void* feVariable, Element_Index lEl_I, double* xi, double* value );

	void _BaseRecoveryFeVar_GetValueInElementWithStdInterpolation( void* feVariable, Element_Index lEl_I, double* xi, double* value );

	double _BaseRecoveryFeVar_ApplyCoeff( BaseRecoveryFeVar* self, double* coeff, double* coord, int order, double* pVec );

	void _BaseRecoveryFeVar_pVec_2Dorder1( double *global, double *pVec );

	void _BaseRecoveryFeVar_pVec_3Dorder1( double *global, double *pVec );

   void BaseUtils_PopulateBoundaryNodesInfo( FeMesh* mesh, BoundaryNodesInfo* bninfo );

#endif

