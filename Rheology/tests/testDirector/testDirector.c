
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Rheology/Rheology.h>

#include <assert.h>

#define TOLERANCE 0.01

double Underworld_testDirector_dt( FiniteElementContext* context ) {
	return 0.05;
}

void Underworld_testDirector_Function( FiniteElementContext* context ) {
	AlignmentSwarmVariable* alignment              = (AlignmentSwarmVariable*)LiveComponentRegister_Get( 
								context->CF->LCRegister, 
								"alignment" );
	FeVariable*             velocityField          = Stg_ComponentFactory_ConstructByName( 
								context->CF, 
								"VelocityField", 
								FeVariable, 
								True,
								0 /* dummy */ );
	Director*               director;
	Particle_Index          lParticle_I;
	GlobalParticle*         particle;
	double                  time                   = context->currentTime + context->dt;
	Swarm*                  swarm;
	double                  error                  = 0.0;
	double                  alignmentValue;
	XYZ                     normal;
	double                  angle                  = 0.5 * M_PI - atan(1.0/(2.0*time) );
	XYZ                     velocity;
	double                  analyticAlignmentValue = 1.0 - cos( angle );
	Stream*                 stream                 = Journal_Register( Info_Type, CURR_MODULE_NAME );

	assert( alignment );

	swarm = alignment->swarm;
	director = alignment->director;

	for ( lParticle_I = 0 ; lParticle_I < swarm->particleLocalCount ; lParticle_I++ ) {
		particle = (GlobalParticle*)Swarm_ParticleAt( swarm, lParticle_I );
		SwarmVariable_ValueAt( alignment, lParticle_I, &alignmentValue );
		SwarmVariable_ValueAt( director->directorSwarmVariable, lParticle_I, normal );

		FieldVariable_InterpolateValueAt( velocityField, particle->coord, velocity );

		error += fabs( alignmentValue - analyticAlignmentValue );
	}
	error /= (double) swarm->particleLocalCount;

	Journal_Printf( stream, "Test %s at timestep %u\n", 
			error < TOLERANCE ? "passed" : "failed", context->timeStep );
	if ( error >= TOLERANCE ) {
		Journal_Printf( stream, "\terror = %g analyticAlignmentValue = %g last alignmentValue = %g angle = %g\n",
				error, analyticAlignmentValue, alignmentValue, angle * 180.0 / M_PI );
	}
	Stream_Flush( stream );
}

const Type Underworld_testDirector_Type = "Underworld_testDirector";
typedef struct {
	__Codelet
} Underworld_testDirector;

void _Underworld_testDirector_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
	FiniteElementContext*   context;
	Stream*                 stream                 = Journal_Register( Info_Type, CURR_MODULE_NAME );

	context = Stg_ComponentFactory_ConstructByName( cf, "context", FiniteElementContext, True, data ); 

	Stream_RedirectFile_WithPrependedPath( stream, context->outputPath, "testDirector.dat" );
	
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput, Underworld_testDirector_Function );
	EP_AppendClassHook( Context_GetEntryPoint( context, FiniteElementContext_EP_CalcDt ), 
			Underworld_testDirector_dt, context );
}

void* _Underworld_testDirector_DefaultNew( Name name ) {
	return Codelet_New(
			Underworld_testDirector_Type,
			_Underworld_testDirector_DefaultNew,
			_Underworld_testDirector_AssignFromXML,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
	}

Index Underworld_testDirector_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Underworld_testDirector_Type, "0", _Underworld_testDirector_DefaultNew );
}
