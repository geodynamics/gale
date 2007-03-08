
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Rheology/Rheology.h>

#include <assert.h>
#include <math.h>

#define TEST_RAND_DIR_ERROR 1.0

void Underworld_testDirectorRandomness_Function( FiniteElementContext* context ) {
	AlignmentSwarmVariable* alignment              = (AlignmentSwarmVariable*) LiveComponentRegister_Get( context->CF->LCRegister, "alignment" );
	Director*               director;
	Particle_Index          lParticle_I;
	GlobalParticle*         particle;
	Swarm*                  swarm;
	XYZ                     normal;
	Stream*                 stream                 = Journal_Register( Info_Type, CURR_MODULE_NAME );
	int                     ii;
	double                  dotProduct;
	double                  angleBetween;
	double                  circleAngle;
	int                     circleAngleCounts[36];
	Bool                    circleErrorFlag = False;
	int                     circleAngleSum;
	double                  circleAngleAverage;
	double                  circleAngleStdDev;
	int                  circleAngleLowerBound;
	int                  circleAngleUpperBound;
	
	assert( alignment );

	swarm = alignment->swarm;
	director = alignment->director;

	for ( ii =0; ii < 36; ii++ ) {	
		circleAngleCounts[ ii ] = 0;
	}

	for ( lParticle_I = 0 ; lParticle_I < swarm->particleLocalCount ; lParticle_I++ ) {
		particle = (GlobalParticle*)Swarm_ParticleAt( swarm, lParticle_I );
		SwarmVariable_ValueAt( director->directorSwarmVariable, lParticle_I, normal );
		/* Calculate dot product between normal and (0,1), then get an angle */
		dotProduct = normal[0] * 0 + normal[1] * 1;
		angleBetween = acos( dotProduct ) * 180 / M_PI; 
		if ( normal[0] > 0 ) {
			circleAngle = angleBetween;
		}
		else {
			circleAngle = 360 - angleBetween;
		}
		circleAngleCounts[ (int)(circleAngle+0.5) / 10 ] += 1;
	}

	Journal_Printf( stream, "Swarm has %d local particles, therefore expect average of %.1f director "
		"normals in each 10 degrees of arc.\n\n", swarm->particleLocalCount,
		(double) swarm->particleLocalCount / 36 );
	circleAngleAverage = (double)swarm->particleLocalCount / 36;
	/*NB. This definition is determined based on a set no. of particlesPerCell.
		  Currently this value is = 20 */
	#define TheoreticalStandardDeviation 13.64
	circleAngleLowerBound = (int)(floor(circleAngleAverage - 3*TheoreticalStandardDeviation));
	if (circleAngleLowerBound < 0 ) {
		circleAngleLowerBound = 0;
	}
	circleAngleUpperBound = (int)(floor(circleAngleAverage + 3*TheoreticalStandardDeviation));
	
	Journal_Printf( stream, "For confidence interval of 99 percent:\n");
    circleAngleSum = 0;
	for ( ii =0; ii < 36; ii++ ) {	
		if (( circleAngleCounts[ii] < circleAngleLowerBound )|| (circleAngleCounts[ii] > circleAngleUpperBound)) {
			Journal_Printf( stream, "Error: Circle angle %d-%d count = %d is outside range [%d, %d].\n"
				, ii*10, ii*10+9, circleAngleCounts[ii], circleAngleLowerBound, circleAngleUpperBound );
			circleErrorFlag = True;			
		}
		circleAngleSum = ((circleAngleAverage - circleAngleCounts[ii]) * 
						(circleAngleAverage - circleAngleCounts[ii])) + circleAngleSum;
	}
	/* Calculate standard deviation */
	circleAngleStdDev = sqrt(circleAngleSum / (36-1));
	//Journal_Printf( stream, "Standard Deviation = %.2f\n", circleAngleStdDev);
	
	
	
	if (fabs(circleAngleStdDev - TheoreticalStandardDeviation) < TEST_RAND_DIR_ERROR ) {
		Journal_Printf( stream, "The standard deviation of director normals is within "
		"tolerance, %.1f, of expected standard deviation, %.1f \n\n", 
		TEST_RAND_DIR_ERROR, TheoreticalStandardDeviation);
	}
	else {
		Journal_Printf( stream, "The standard deviation of director normals is not within"
		" tolerance, %.1f, of expected standard deviation, %.1f, instead is = %.1f \n\n", 
		TEST_RAND_DIR_ERROR, TheoreticalStandardDeviation, circleAngleStdDev);
	}
	
	if (circleErrorFlag == False) {
		Journal_Printf( stream, "All circle angle counts are within range.\n");
	}
	else {
		Journal_Printf( stream, "Some circle angle counts not within range \n");
	}
	
	Stream_Flush( stream );
}

const Type Underworld_testDirectorRandomness_Type = "Underworld_testDirectorRandomness";
typedef struct {
	__Codelet
} Underworld_testDirectorRandomness;

void _Underworld_testDirectorRandomness_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	FiniteElementContext*   context;
	Stream*                 stream                 = Journal_Register( Info_Type, CURR_MODULE_NAME );

	context = Stg_ComponentFactory_ConstructByName( cf, "context", FiniteElementContext, True, data ); 

	Stream_RedirectFile_WithPrependedPath( stream, context->outputPath, "testDirectorRandomness.dat" );
	
	ContextEP_Append_AlwaysLast( context, AbstractContext_EP_FrequentOutput,
		Underworld_testDirectorRandomness_Function );
}

void* _Underworld_testDirectorRandomness_DefaultNew( Name name ) {
	return Codelet_New(
			Underworld_testDirectorRandomness_Type,
			_Underworld_testDirectorRandomness_DefaultNew,
			_Underworld_testDirectorRandomness_Construct,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
	}

Index Underworld_testDirectorRandomness_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Underworld_testDirectorRandomness_Type, "0", _Underworld_testDirectorRandomness_DefaultNew );
}
