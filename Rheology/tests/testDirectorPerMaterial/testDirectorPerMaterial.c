
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Rheology/Rheology.h>

#include <assert.h>

#define DIR_TEST_NUM_MAT 3
#define DIR_X_AXIS 0
#define DIR_Y_AXIS 1
#define DIR_Z_AXIS 2
#define ANGLE_ERROR 1e-15

void Underworld_testDirectorPerMaterial_Function( PICelleratorContext* context ) {
	AlignmentSwarmVariable* alignment              = (AlignmentSwarmVariable*) LiveComponentRegister_Get( context->CF->LCRegister, "alignment" );
	Materials_Register*     materials_Register     = context->materials_Register;
	Director*               director;
	Particle_Index          lParticle_I;
	Swarm*                  swarm;
	XYZ                     testVector;
	Stream*                 stream                 = Journal_Register( Info_Type, CURR_MODULE_NAME );
	Material_Index          materialsCount;
	int                     material_I;
	XYZ*                    matDirectionVectors;
	double                  angle;
	XYZI                    materialErrorCounts;
	XYZ                     accumulatedErrors;
	
	assert( alignment );

	swarm = alignment->swarm;
	director = alignment->director;
	materialsCount = Materials_Register_GetCount( materials_Register);
	
	/*  construct test for testDirectorPerMaterial.xml  */
	/* assume a direction for each material and check that */
	/* each particle has the same direction for it's material type. */
	/* and that it's equal to the default value for each material. */
	/* get materials and their initial directions */
	
	/* check that there are at least 3 materials */
	/* check that the initialDirections are == certain defaults */
	if (materialsCount != DIR_TEST_NUM_MAT ) {
		Journal_Printf(stream, "Amount of materials does not match default of %d\n", DIR_TEST_NUM_MAT);
	}
	matDirectionVectors = Memory_Alloc_Array(XYZ, DIR_TEST_NUM_MAT, "materialDirectionVectors");

	/* Define vectors */
	
	matDirectionVectors[0][DIR_X_AXIS] = 1.0;
	matDirectionVectors[0][DIR_Y_AXIS] = 0.0;
	matDirectionVectors[0][DIR_Z_AXIS] = 0.0;
	
	matDirectionVectors[1][DIR_X_AXIS] = 0.0;
	matDirectionVectors[1][DIR_Y_AXIS] = 1.0;
	matDirectionVectors[1][DIR_Z_AXIS] = 0.0;
	
	matDirectionVectors[2][DIR_X_AXIS] = 1.0;
	matDirectionVectors[2][DIR_Y_AXIS] = 1.0;
	matDirectionVectors[2][DIR_Z_AXIS] = 0.0;
	
	/* set error counters */
	for (material_I=0; material_I < DIR_TEST_NUM_MAT; material_I++) {
		materialErrorCounts[material_I] = 0;
		accumulatedErrors[material_I] = 0;
	}		

	for ( lParticle_I = 0 ; lParticle_I < swarm->particleLocalCount ; lParticle_I++ ) {

		/* get particle's material type */
		material_I = MaterialPointsSwarm_GetMaterialIndexAt(
						swarm, 
						lParticle_I );
		/* get particle's initialDirection */
		SwarmVariable_ValueAt( director->directorSwarmVariable, lParticle_I, testVector );
	
		/* check that angle is within error bar of correct value */
		angle = StGermain_AngleBetweenVectors(matDirectionVectors[material_I],testVector, swarm->dim);
		
		if (fabs(angle) > ANGLE_ERROR) {
			materialErrorCounts[material_I] = materialErrorCounts[material_I] + 1;
			accumulatedErrors[material_I] = accumulatedErrors[material_I] + 1;
		}	
		
	}
	
	/* if all particles within tolerance, then print out a 'pass' message */
	for (material_I = 0; material_I < DIR_TEST_NUM_MAT; material_I++) {
		if (materialErrorCounts[material_I] == 0 ) {		
			Journal_Printf(stream, "Angles for material %d within tolerance %g\n", 
				material_I, ANGLE_ERROR);
		}
		else {
			Journal_Printf(stream, "Angles for material %d not within tolerance %g\n", 
				material_I, ANGLE_ERROR);
			Journal_Printf(stream, "%d particles outside tolerance, with accumulated error = %f\n", 
				materialErrorCounts[material_I], accumulatedErrors[material_I]);
		}
	}
	Stream_Flush( stream );
}

const Type Underworld_testDirectorPerMaterial_Type = "Underworld_testDirectorPerMaterial";
typedef struct {
	__Codelet
} Underworld_testDirectorPerMaterial;

void _Underworld_testDirectorPerMaterial_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	FiniteElementContext*   context;
	Stream*                 stream                 = Journal_Register( Info_Type, CURR_MODULE_NAME );

	context = Stg_ComponentFactory_ConstructByName( cf, "context", FiniteElementContext, True, data );

	Stream_RedirectFile_WithPrependedPath( stream, context->outputPath, "testDirectorPerMaterial.dat" );
	
	ContextEP_Append_AlwaysLast( context, AbstractContext_EP_FrequentOutput,
		Underworld_testDirectorPerMaterial_Function );
}

void* _Underworld_testDirectorPerMaterial_DefaultNew( Name name ) {
	return Codelet_New(
			Underworld_testDirectorPerMaterial_Type,
			_Underworld_testDirectorPerMaterial_DefaultNew,
			_Underworld_testDirectorPerMaterial_Construct,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
	}

Index Underworld_testDirectorPerMaterial_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Underworld_testDirectorPerMaterial_Type, "0", _Underworld_testDirectorPerMaterial_DefaultNew );
}
