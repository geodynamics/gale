

#ifndef __Underworld_Rheology_DepthDependentViscosity_h__
#define __Underworld_Rheology_DepthDependentViscosity_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type DepthDependentViscosity_Type;

	/** Rheology class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __DepthDependentViscosity \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__Rheology \
		/* Virtual functions go here */ \
		/* Material Parameters */\
		FeMesh*    		                            feMesh;                             \
		double                                              eta0;                               \
		double                                              gamma;                              \
		Axis                                                variationAxis;                      \
		double                                              referencePoint;                     

	struct DepthDependentViscosity { __DepthDependentViscosity };

	/** Public Constructor */
   DepthDependentViscosity* DepthDependentViscosity_New(
      Name                  name,
      AbstractContext*      context,
      FeMesh*               feMesh,
      double                eta0,
      double                gamma,
      Axis                  variationAxis,
      double                referencePoint );
	
	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define DEPTHDEPENDENTVISCOSITY_DEFARGS \
                RHEOLOGY_DEFARGS

	#define DEPTHDEPENDENTVISCOSITY_PASSARGS \
                RHEOLOGY_PASSARGS

	DepthDependentViscosity* _DepthDependentViscosity_New(  DEPTHDEPENDENTVISCOSITY_DEFARGS  );

	
	/* 'Stg_Component' implementations */
	void* _DepthDependentViscosity_DefaultNew( Name name ) ;
	void _DepthDependentViscosity_AssignFromXML( void* rheology, Stg_ComponentFactory* cf, void* data );
   void _DepthDependentViscosity_Init( DepthDependentViscosity* self, FeMesh* feMesh, double eta0, double gamma, Axis variationAxis, double referencePoint );
   void _DepthDependentViscosity_Destroy( void* rheology, void* data ) ;
   
	void _DepthDependentViscosity_ModifyConstitutiveMatrix( 
		void*                                              rheology, 
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               swarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi );
#endif

