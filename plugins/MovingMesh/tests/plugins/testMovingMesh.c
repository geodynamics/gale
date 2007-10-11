
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include <math.h>
#include <string.h>
#include <assert.h>

ExtensionInfo_Index handle;

double dt( void* class, PICelleratorContext* context ) {
	return Dictionary_GetDouble_WithDefault( context->dictionary, "dt", 0.01 );
}

void MovingMeshTestVelBCs_IncreasingWithY( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	FeVariable*             feVariable         = NULL;
	FeMesh*     		mesh               = NULL;
	double*                 result             = (double*) _result;
	double			min[3], max[3];
	double*                 coord;
	
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh       = feVariable->feMesh;

	/* Find coordinate of node */
	Mesh_GetGlobalCoordRange( mesh, min, max );
	coord = Mesh_GetVertex( mesh, node_lI );
	
	*result = ( coord[J_AXIS] - min[J_AXIS] ) / ( max[J_AXIS] - min[J_AXIS] );
}


void MovingMeshTestVelBCs_ToCentreY( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	FeVariable*             feVariable         = NULL;
	FeMesh* 		mesh               = NULL;
	double*                 result             = (double*) _result;
	double*                 coord;
	double                  toCentreJ_Max = 0;
	double                  boxHeight = 0;
	double                  centreInJ = 0;
	double                  coordJ_RelativeToCentreJ = 0;
	double			min[3], max[3];
	
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh       = feVariable->feMesh;

	/* Find coordinate of node */
	Mesh_GetGlobalCoordRange( mesh, min, max );
	coord = Mesh_GetVertex( mesh, node_lI );
	
	boxHeight = max[J_AXIS] - min[J_AXIS];
	centreInJ = ( max[J_AXIS] + min[J_AXIS] ) / 2;
	coordJ_RelativeToCentreJ = coord[J_AXIS] - centreInJ;

	toCentreJ_Max  = Dictionary_GetDouble_WithDefault( context->dictionary, "toCentreJ_Max", 1.0 );
	/* Need the -1 factor to "flip" so both above and below coords go towards the centre */
	*result = (-1) * coordJ_RelativeToCentreJ / (boxHeight / 2) * toCentreJ_Max;
}


void construct( PICelleratorContext* context ) {
	MaterialPointsSwarm* materialSwarm;

	/* Add original pos to particle */
	materialSwarm = (MaterialPointsSwarm*) LiveComponentRegister_Get( context->CF->LCRegister, "materialSwarm" );
	handle = ExtensionManager_Add( materialSwarm->particleExtensionMgr, CURR_MODULE_NAME, sizeof( Coord ) );

	ConditionFunction_Register_Add( context->condFunc_Register,
		ConditionFunction_New( MovingMeshTestVelBCs_IncreasingWithY, "MovingMeshTestVelBCs_IncreasingWithY" ) );
	ConditionFunction_Register_Add( context->condFunc_Register,
		ConditionFunction_New( MovingMeshTestVelBCs_ToCentreY, "MovingMeshTestVelBCs_ToCentreY" ) );

	/* Need to do this so we can test the remesher handles incorrect BCs properly */	
	stJournal->firewallProducesAssert = False; 
}

void storeOriginalPos( PICelleratorContext* context ) {
	MaterialPointsSwarm*   materialSwarm;
	GlobalParticle*   particle;
	double*           originalCoord;
	Particle_Index    lParticle_I;

	materialSwarm = (MaterialPointsSwarm*) LiveComponentRegister_Get( context->CF->LCRegister, "materialSwarm" );

	for ( lParticle_I = 0 ; lParticle_I < materialSwarm->particleLocalCount ; lParticle_I++ ) {
		particle = (GlobalParticle*) Swarm_ParticleAt( materialSwarm, lParticle_I );
		originalCoord = ExtensionManager_Get( materialSwarm->particleExtensionMgr, particle, handle );

		memcpy( originalCoord, particle->coord, sizeof(Coord) );
	}
}

void check( PICelleratorContext* context ) {
	MaterialPointsSwarm*   materialSwarm;
	GlobalParticle*   particle;
	double*           originalCoord;
	double*           coord;
	Particle_Index    lParticle_I;
	double            maxError             = 0.0;
	double            errorGlobal;
	double            currLeft;
	double            currRight;
	double            startLeft;
	double            startRight;
	double            velocityLeft;
	double            velocityRight;
	double            analyticCoord;
	/* Since this test is occuring _after_ time integration using the latest dt, have to add dt on to the
		total time up to the start of this timestep */
	double            time                 = context->currentTime + context->dt;
	Dictionary*       dictionary           = context->dictionary;
	Stream*           stream               = Journal_Register( Info_Type, CURR_MODULE_NAME );
	double            tolerance;
	Dimension_Index   dim;

	/* BIG HACK - SHOULD DO PROPERLEY */
	startLeft  = Dictionary_GetDouble( dictionary, "minX" );
	startRight = Dictionary_GetDouble( dictionary, "maxX" );
	velocityLeft  = -0.5;
	velocityRight = +0.5;
	dim = 2;
	
	/* Add original pos to particle */
	materialSwarm = (MaterialPointsSwarm*) LiveComponentRegister_Get( context->CF->LCRegister, "materialSwarm" );

	tolerance = Dictionary_GetDouble_WithDefault( dictionary, "tolerance", 1.0e-6 );

	Journal_Printf( stream, "Timestep = %u: ", context->timeStep );

	for ( lParticle_I = 0 ; lParticle_I < materialSwarm->particleLocalCount ; lParticle_I++ ) {
		particle      = (GlobalParticle*)Swarm_ParticleAt( materialSwarm, lParticle_I );

		coord         = particle->coord;
		originalCoord = ExtensionManager_Get( materialSwarm->particleExtensionMgr, particle, handle );

		currLeft  = startLeft + velocityLeft * time;
		currRight = startRight + velocityRight * time;

		analyticCoord = 
			currLeft + ( originalCoord[ I_AXIS ] - startLeft ) * ( currRight - currLeft ) / ( startRight - startLeft );

		maxError = MAX( maxError, fabs( analyticCoord - coord[ I_AXIS ] ) );

		Journal_Firewall( maxError < tolerance, stream, 
			"Failed - for local particle %u, orig coord was (%.3f,%.3f,%.3f), thus given "
			"startLeft of %.3f, startRight of %.3f, velocityLeft of %.3f, velocityRight of "
			"%.3f, and currTime of %.2f (calc. currLeft is %.3f and currRight is %.3f):\n"
			"\tanalytic coord in I axis is %.3f, but actual coord in I is %.3f.\n",
			lParticle_I, originalCoord[0], originalCoord[1], originalCoord[2],
			startLeft, startRight, velocityLeft, velocityRight, time, currLeft, currRight,
			analyticCoord, coord[I_AXIS] );

		if ( maxError > 0.45 )
			abort();
	}

	MPI_Allreduce( &maxError,  &errorGlobal,  1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );

	if ( errorGlobal < tolerance )
		Journal_Printf( stream, "Passed\n" );
	else 
		Journal_Printf( stream, "Failed - Error = %g - Tolerance = %g\n", errorGlobal, tolerance );
}

const Type testMovingMesh_Type = "testMovingMesh";
typedef struct {
	__Codelet
} testMovingMesh;

void _testMovingMesh_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	DomainContext* context;
	Stream*                stream               = Journal_Register( Info_Type, CURR_MODULE_NAME );

	context = Stg_ComponentFactory_ConstructByName( cf, "context", DomainContext, True, data ); 

	ContextEP_Append( context, AbstractContext_EP_ConstructExtensions, construct );
	ContextEP_Append( context, AbstractContext_EP_Initialise, storeOriginalPos );
	ContextEP_Prepend( context, AbstractContext_EP_Step, check );
	EP_AppendClassHook( Context_GetEntryPoint( context, FiniteElementContext_EP_CalcDt ), dt, context );

	Stream_RedirectFile_WithPrependedPath( stream, context->outputPath, "output.dat" );
	Stream_SetPrintingRank( stream, 0 );
	
}

void* _testMovingMesh_DefaultNew( Name name ) {
	return Codelet_New(
			testMovingMesh_Type,
			_testMovingMesh_DefaultNew,
			_testMovingMesh_Construct,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
}

Index testMovingMesh_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, testMovingMesh_Type, "0", _testMovingMesh_DefaultNew );
}
