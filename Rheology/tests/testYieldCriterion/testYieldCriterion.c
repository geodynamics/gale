
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Rheology/Rheology.h>

#include <assert.h>

const Type Underworld_testYieldCriterion_Type = "Underworld_testYieldCriterion";

/* Define plugin structure */
typedef struct {
	__Codelet
	YieldRheology_HasYieldedFunction* realHasYieldedFunction;
	FeMesh* 		          mesh;
	XYZ                               min;
	XYZ                               max;
	Bool                              hasYielded;
} Underworld_testYieldCriterion;
	
FiniteElementContext*   context;

void testYieldCriterion_HasYielded( 
		void*                            yieldRheology,
		ConstitutiveMatrix*              constitutiveMatrix,
		MaterialPointsSwarm*             materialPointsSwarm,
		Element_LocalIndex               lElement_I,
		MaterialPoint*                   materialPoint,
		double                           yieldCriterion,
		double                           yieldIndicator ) 
{
	Underworld_testYieldCriterion* self;
	Dimension_Index dim_I;

	/* Get pointer to plugin struct */
	self = (Underworld_testYieldCriterion*) LiveComponentRegister_Get(
					context->CF->LCRegister,
					Underworld_testYieldCriterion_Type );

	/* Call real 'HasYielded' function */
	self->realHasYieldedFunction( 
			yieldRheology, constitutiveMatrix, materialPointsSwarm, lElement_I, materialPoint, yieldCriterion, yieldIndicator );

	/* Don't output information if this is the first non-linear iteration */
	if ( constitutiveMatrix->sleNonLinearIteration_I == 0 ) {
		return;
	}

	/* Store information */
	self->hasYielded = True;
	for ( dim_I = 0 ; dim_I < context->dim ; dim_I++ ) {
		if ( materialPoint->coord[ dim_I ] < self->min[ dim_I ] )
			self->min[ dim_I ] = materialPoint->coord[ dim_I ];
		if ( materialPoint->coord[ dim_I ] > self->max[ dim_I ] )
			self->max[ dim_I ] = materialPoint->coord[ dim_I ];
	}
}

double Underworld_testYieldCriterion_dt( FiniteElementContext* context ) {
	if ( context->currentTime >= 0.65 ) {
		return 0.01;
	}
	return HUGE_VAL;
}

void Underworld_testYieldCriterion_Check( FiniteElementContext* context ) {
	Stream* stream = Journal_Register( Dump_Type, Underworld_testYieldCriterion_Type );
	Underworld_testYieldCriterion* self;

	/* Get pointer to plugin struct */
	self = (Underworld_testYieldCriterion*) LiveComponentRegister_Get(
					context->CF->LCRegister,
					Underworld_testYieldCriterion_Type );
	assert( self );

	/* Don't do anything if nothing has yielded yet */
	if ( !self->hasYielded ) {
		return;
	}

	/* Get Calculation to stop */
	context->maxTimeSteps = context->timeStep;

	/* Set the stream to point to our output file (so we can do a diff on it later) */
	Stream_Enable( stream, True );
	Stream_RedirectFile_WithPrependedPath( stream, context->outputPath, "testYieldCriterion.dat" );

	Journal_Printf( stream, "Material yielded at time %4g (step %u) within:\n", context->currentTime, context->timeStep ); 

	/* Output information */
	Journal_Printf( stream, "\tx: %12.4g - %12.4g\n", self->min[ I_AXIS ], self->max[ I_AXIS ] );
	Journal_Printf( stream, "\ty: %12.4g - %12.4g\n", self->min[ J_AXIS ], self->max[ J_AXIS ] );
	if ( context->dim == 3 ) {
		Journal_Printf( stream, "\tz: %12.4g - %12.4g\n", self->min[ K_AXIS ], self->max[ K_AXIS ] );
	}
}

void _Underworld_testYieldCriterion_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
	YieldRheology*          yieldRheology;
	Underworld_testYieldCriterion* self;

	context = Stg_ComponentFactory_ConstructByName( cf, "context", FiniteElementContext, True, data ); 

	/* Get pointer to plugin struct */
	self = (Underworld_testYieldCriterion*) LiveComponentRegister_Get(
					context->CF->LCRegister,
					Underworld_testYieldCriterion_Type );

	/* get pointer to the mesh */
	self->mesh = Stg_ComponentFactory_ConstructByName( cf, "linearMesh", FeMesh, True, data ); 
	
	/* Get a pointer the yield rheology that we are trying to test */
	yieldRheology = (YieldRheology*) LiveComponentRegister_Get( context->CF->LCRegister, "yieldRheology" );
	
	/* Store the pointer to the original 'HasYielded' function */
	self->realHasYieldedFunction = yieldRheology->_hasYielded;

	/* Reset this function pointer with our own */
	yieldRheology->_hasYielded = testYieldCriterion_HasYielded;

	ContextEP_Append( context, AbstractContext_EP_Step, Underworld_testYieldCriterion_Check );
	EP_AppendClassHook( Context_GetEntryPoint( context, FiniteElementContext_EP_CalcDt ), 
			Underworld_testYieldCriterion_dt, context );

	self->min[ I_AXIS ] = HUGE_VAL;
	self->min[ J_AXIS ] = HUGE_VAL;
	self->min[ K_AXIS ] = HUGE_VAL;
	self->max[ I_AXIS ] = -HUGE_VAL;
	self->max[ J_AXIS ] = -HUGE_VAL;
	self->max[ K_AXIS ] = -HUGE_VAL;
}

void* _Underworld_testYieldCriterion_DefaultNew( Name name ) {
	return _Codelet_New(
			sizeof( Underworld_testYieldCriterion ),
			Underworld_testYieldCriterion_Type,
			_Codelet_Delete,
			_Codelet_Print,
			_Codelet_Copy,
			_Underworld_testYieldCriterion_DefaultNew,
			_Underworld_testYieldCriterion_AssignFromXML,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
	}

Index Underworld_testYieldCriterion_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Underworld_testYieldCriterion_Type, "0", _Underworld_testYieldCriterion_DefaultNew );
}
