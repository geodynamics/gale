
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Rheology/Rheology.h>

#include <assert.h>

typedef struct {
	__AnalyticSolution
	MaterialViscosity*               materialViscosity;
	NonNewtonian*                    nonNewtonianRheology;
	BlockGeometry*                   geometry;
	double                           velocityTopOfBox;
} NonNewtonianShearSolution;

const Type NonNewtonianShearSolution_Type = "NonNewtonianShearSolution";

void NonNewtonianShearSolution_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity ) {
	NonNewtonianShearSolution*   self = (NonNewtonianShearSolution*)analyticSolution;
	double                       height;
	
	/* Get Parameters */
	height = self->geometry->max[J_AXIS] - self->geometry->min[J_AXIS];

	velocity[ I_AXIS ] = coord[J_AXIS] / height * self->velocityTopOfBox;
	velocity[ J_AXIS ] = 0.0;
}

void NonNewtonianShearSolution_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate ) {
	NonNewtonianShearSolution*   self = (NonNewtonianShearSolution*)analyticSolution;
	double                       height;
	
	/* Get Parameters */
	height = self->geometry->max[J_AXIS] - self->geometry->min[J_AXIS];

	strainRate[ 0 ] = 0.0;
	strainRate[ 1 ] = 0.0;
	strainRate[ 2 ] = 0.5 * self->velocityTopOfBox / height;
}

void NonNewtonianShearSolution_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress ) {
	NonNewtonianShearSolution*   self = (NonNewtonianShearSolution*)analyticSolution;
	double                         eta;
	double                         height;
	double                         tau;
	double                         n;
	
	/* Get Parameters */
	eta = self->materialViscosity->eta0;
	n = self->nonNewtonianRheology->stressExponent;
	height = self->geometry->max[J_AXIS] - self->geometry->min[J_AXIS];

	/* Calculate stress - without considering cohesion */
	tau = pow( eta * self->velocityTopOfBox / height, 1/n);

	/* Calculate Analytic Solution for the stress */
	stress[0] = 0.0;
	stress[1] = 0.0;
	stress[2] = tau;
}


void NonNewtonianShearSolution_ViscosityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* viscosity ) {
	NonNewtonianShearSolution*     self = (NonNewtonianShearSolution*)analyticSolution;
	double                         eta0;
	double                         height;
	double                         n;
	
	/* Get Parameters */
	eta0 = self->materialViscosity->eta0;
	n = self->nonNewtonianRheology->stressExponent;
	height = self->geometry->max[J_AXIS] - self->geometry->min[J_AXIS];

	/* Calculate stress - without considering cohesion */
	*viscosity = eta0 * pow( eta0 * self->velocityTopOfBox / height, 1/n-1.0);
}

void NonNewtonianShearSolution_UpdateVelocityBC( NonNewtonianShearSolution* self, FiniteElementContext* context ) {
	self->velocityTopOfBox += context->dt;
}

void NonNewtonianShearSolution_VelocityBC( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	FiniteElementContext *	         context       = (FiniteElementContext*)_context;
	NonNewtonianShearSolution*   self          = (NonNewtonianShearSolution*) LiveComponentRegister_Get( context->CF->LCRegister, NonNewtonianShearSolution_Type );
	double*                          result        = (double*) _result;

	*result = self->velocityTopOfBox;
}

void _NonNewtonianShearSolution_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	NonNewtonianShearSolution* self = (NonNewtonianShearSolution*)analyticSolution;
	FiniteElementContext*          context;
	FeVariable*                    velocityField;
	FeVariable*                    strainRateField;
	FeVariable*                    stressField;
	FeVariable*                    viscosityField;

	/* Construct parent */
	_AnalyticSolution_Construct( self, cf, data );
	context = Stg_CheckType( self->context, FiniteElementContext );
	
	ConditionFunction_Register_Add( context->condFunc_Register, 
			ConditionFunction_New( NonNewtonianShearSolution_VelocityBC, "ShearTrigger") );	

	/* Create Analytic Velocity Field */
	velocityField = Stg_ComponentFactory_ConstructByName( cf, "VelocityField", FeVariable, True, data );
	AnalyticSolution_CreateAnalyticVectorField( self, velocityField, NonNewtonianShearSolution_VelocityFunction );

	/* Create Analytic Strain Rate Field */
	strainRateField = Stg_ComponentFactory_ConstructByName( cf, "StrainRateField", FeVariable, True, data );
	AnalyticSolution_CreateAnalyticSymmetricTensorField( 
			self, strainRateField, NonNewtonianShearSolution_StrainRateFunction );

	/* Create Analytic Stress Field */
	stressField = Stg_ComponentFactory_ConstructByName( cf, "StressField", FeVariable, True, data );
	AnalyticSolution_CreateAnalyticSymmetricTensorField( self, stressField, NonNewtonianShearSolution_StressFunction );

	/* Create Analytic Viscosity Field */
	viscosityField = Stg_ComponentFactory_ConstructByName( cf, "ViscosityField", FeVariable, True, data );
	AnalyticSolution_CreateAnalyticField( self, viscosityField, NonNewtonianShearSolution_ViscosityFunction );

	self->materialViscosity = Stg_ComponentFactory_ConstructByName( cf, "layerViscosity", MaterialViscosity, True, data );
	self->nonNewtonianRheology = 
		Stg_ComponentFactory_ConstructByName( cf, "nonNewtonianRheology", NonNewtonian, True, data );
	self->geometry = Stg_ComponentFactory_ConstructByName( cf, "geometry", BlockGeometry, True, data );

	/* Set Velocity Stuff */
	EP_AppendClassHook( Context_GetEntryPoint( context, AbstractContext_EP_UpdateClass ),
								NonNewtonianShearSolution_UpdateVelocityBC, self );
	self->velocityTopOfBox = Stg_ComponentFactory_GetRootDictDouble( cf, "velocityTopOfBox", 0.5 );
}


/* This function will provide StGermain the abilty to instantiate (create) this codelet on demand. */
void* _NonNewtonianShearSolution_DefaultNew( Name name ) {
	return _AnalyticSolution_New(
			sizeof(NonNewtonianShearSolution),
			NonNewtonianShearSolution_Type,
			_AnalyticSolution_Delete,
			_AnalyticSolution_Print,
			_AnalyticSolution_Copy,
			_NonNewtonianShearSolution_DefaultNew,
			_NonNewtonianShearSolution_Construct,
			_AnalyticSolution_Build,
			_AnalyticSolution_Initialise,
			_AnalyticSolution_Execute,
			_AnalyticSolution_Destroy,
			name );	
}
	
/* This function is automatically run by StGermain when this plugin is loaded. The name must be "<plugin-name>_Register". */
Index Underworld_NonNewtonianShearSolution_Register( PluginsManager* pluginsManager ) {
	/* A plugin is only properly registered once it returns the handle provided when submitting a codelet to StGermain. */
	return PluginsManager_Submit( pluginsManager, NonNewtonianShearSolution_Type, "0",
		_NonNewtonianShearSolution_DefaultNew );
}

