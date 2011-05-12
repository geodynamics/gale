/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2005, Monash Cluster Computing 
** All rights reserved.
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
** 		* Redistributions of source code must retain the above copyright notice, 
** 			this list of conditions and the following disclaimer.
** 		* Redistributions in binary form must reproduce the above copyright 
**			notice, this list of conditions and the following disclaimer in the 
**			documentation and/or other materials provided with the distribution.
** 		* Neither the name of the Monash University nor the names of its contributors 
**			may be used to endorse or promote products derived from this software 
**			without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
** THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
** PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS 
** BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
** CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
** SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
** HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
** LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
** OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
**
** Contact:
*%		Louis Moresi - Louis.Moresi@sci.monash.edu.au
*%
** Contributors:
*+		Robert Turnbull
*+		Vincent Lemiale
*+		Louis Moresi
*+		David May
*+		David Stegman
*+		Mirko Velic
*+		Patrick Sunter
*+		Julian Giordani
*+		Kathleen Humble
*+
** $Id: Director.c 629 2007-11-14 05:47:33Z BelindaMay $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "RheologyClass.h"
#include "YieldRheology.h"
#include "VonMises.h"
#include "Director.h"
#include "ConstitutiveMatrix.h"

#include <assert.h>
#include <string.h>

/* Define the default initial direction for particles */ 
#define DIRECTOR_DEFAULT_DIR_X 0.0
#define DIRECTOR_DEFAULT_DIR_Y 1.0
#define DIRECTOR_DEFAULT_DIR_Z 0.0


/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type Director_Type = "Director";

/* Public Constructor */
Director* Director_New( 
		Name                   name,
		DomainContext*         context,
		TimeIntegrator*        timeIntegrator, 
		Variable*              variable,
		Index                  dataCount, 
		Stg_Component**        data,
		Bool                   allowFallbackToFirstOrder,
		FeVariable*            velGradField,
		MaterialPointsSwarm*   materialPointsSwarm,
		InitialDirectionType   initialDirectionType,
		double                 globalInitialDirectionX,
		double                 globalInitialDirectionY,
		double                 globalInitialDirectionZ,
		int                    randomInitialDirectionSeed,
		Bool                   dontUpdate )
{
	Director*	self;

	self = (Director*) _Director_DefaultNew( name );
	
	_TimeIntegrand_Init( self, context, timeIntegrator, variable, dataCount, data, allowFallbackToFirstOrder );
   _Director_Init(
		self,
		velGradField,
		materialPointsSwarm,
		initialDirectionType,
		globalInitialDirectionX,
		globalInitialDirectionY,
		globalInitialDirectionZ,
		randomInitialDirectionSeed,
		dontUpdate );

	self->isConstructed = True;
	return self;
}

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
Director* _Director_New(  DIRECTOR_DEFARGS  ) 
{
	Director*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(Director) );
	self = (Director*) _TimeIntegrand_New(  TIMEINTEGRAND_PASSARGS  );
	
	/* Function pointers for this class that are not on the parent class should be set here */
	
	return self;
}

void _Director_Init(
		Director*                                          self,
		FeVariable*                                        velGradField,
		MaterialPointsSwarm*                               materialPointsSwarm,
		InitialDirectionType                               initialDirectionType,
		double                                             globalInitialDirectionX,
		double                                             globalInitialDirectionY,
		double                                             globalInitialDirectionZ,
		int                                                randomInitialDirectionSeed,
		Bool                                               dontUpdate )
{
	/* Assign Values */
	self->velGradField               = velGradField;
	self->materialPointsSwarm        = materialPointsSwarm;
	self->dontUpdate                 = dontUpdate;
	self->globalInitialDirection[ I_AXIS ] = globalInitialDirectionX;
	self->globalInitialDirection[ J_AXIS ] = globalInitialDirectionY;
	self->globalInitialDirection[ K_AXIS ] = globalInitialDirectionZ;
	self->randomInitialDirectionSeed = randomInitialDirectionSeed;
	StGermain_VectorNormalise( self->globalInitialDirection, materialPointsSwarm->dim );
	self->initialDirectionType = initialDirectionType;
	/****** Setup Variables *****/

	/* First check to see if a particle extension has already been created for this swarm */
	self->particleExtHandle     = ExtensionManager_GetHandle( materialPointsSwarm->particleExtensionMgr, (Name)Director_Type );

	/* If there isn't one then create the particle extension - otherwise just use the one already there*/
	if ( self->particleExtHandle == (ExtensionInfo_Index) -1  ) {
		StandardParticle      particle;
		Director_ParticleExt* particleExt;

		/* Add particle extension */
		self->particleExtHandle = 
			ExtensionManager_Add( materialPointsSwarm->particleExtensionMgr, (Name)Director_Type, sizeof(Director_ParticleExt)  );	

		particleExt = (Director_ParticleExt*)ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, &particle, self->particleExtHandle );

		self->directorSwarmVariable = Swarm_NewVectorVariable( materialPointsSwarm, (Name)"Director", (ArithPointer) &particleExt->director - (ArithPointer) &particle,
			Variable_DataType_Double,
			materialPointsSwarm->dim,
			"DirectorX",
			"DirectorY",
			"DirectorZ" );
		self->dontUpdateParticle = Swarm_NewScalarVariable( materialPointsSwarm, (Name)"dontUpdateParticle", (ArithPointer) &particleExt->dontUpdateParticle - (ArithPointer) &particle, Variable_DataType_Int  );
	}
	else {
		char* variableName;

		/* Get Variables already created */
		variableName = Stg_Object_AppendSuffix( materialPointsSwarm, (Name)"Director"  );
		self->directorSwarmVariable = SwarmVariable_Register_GetByName( materialPointsSwarm->swarmVariable_Register, variableName );
		assert( self->directorSwarmVariable );
		Memory_Free( variableName );
		/* Get Variables already created */
		variableName = Stg_Object_AppendSuffix( materialPointsSwarm, (Name)"dontUpdateParticle"  );
		self->dontUpdateParticle = SwarmVariable_Register_GetByName( materialPointsSwarm->swarmVariable_Register, variableName );
		assert( self->dontUpdateParticle );
		Memory_Free( variableName );
	}
	
	self->variable = self->directorSwarmVariable->variable;

	TimeIntegrator_AppendSetupEP( self->timeIntegrator,
                                      "Director_UpdateVariables", (void*)Director_UpdateVariables,  self->name, self );
}

void* _Director_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                               _sizeOfSelf = sizeof(Director);
	Type                                                       type = Director_Type;
	Stg_Class_DeleteFunction*                               _delete = _TimeIntegrand_Delete;
	Stg_Class_PrintFunction*                                 _print = _TimeIntegrand_Print;
	Stg_Class_CopyFunction*                                   _copy = _TimeIntegrand_Copy;
	Stg_Component_DefaultConstructorFunction*   _defaultConstructor = _Director_DefaultNew;
	Stg_Component_ConstructFunction*                     _construct = _Director_AssignFromXML;
	Stg_Component_BuildFunction*                             _build = _Director_Build;
	Stg_Component_InitialiseFunction*                   _initialise = _Director_Initialise;
	Stg_Component_ExecuteFunction*                         _execute = _TimeIntegrand_Execute;
	Stg_Component_DestroyFunction*                         _destroy = _Director_Destroy;
	TimeIntegrand_CalculateTimeDerivFunction*  _calculateTimeDeriv = _Director_TimeDerivative;
	TimeIntegrand_IntermediateFunction*              _intermediate = _Director_Intermediate;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _Director_New(  DIRECTOR_PASSARGS  );
}

void _Director_AssignFromXML( void* director, Stg_ComponentFactory* cf, void* data ){
	Director*               self           = (Director*)director;
	MaterialPointsSwarm*    materialPointsSwarm;
	FeVariable*             velGradField;
	char*                   initialDirectionTypeName;
	InitialDirectionType    initialDirectionType;
	
	/* Construct Parent */
	_TimeIntegrand_AssignFromXML( self, cf, data );
	
	/* Construct 'Director' stuff */
	/* TODO: 'KeyFallback' soon to be deprecated/updated */
        velGradField   = Stg_ComponentFactory_ConstructByNameWithKeyFallback( cf, self->name, (Name)"VelocityGradientsField", (Dictionary_Entry_Key)"VelocityGradientsField", FeVariable, True, data  );
	/*
	velGradField   = Stg_ComponentFactory_ConstructByKey( 
			cf, self->name, "VelocityGradientsField", FeVariable, True );
	*/		
	materialPointsSwarm = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"MaterialPointsSwarm", MaterialPointsSwarm, True, data  );
	
	
	/* Define the initial Direction type for the problem
		The options are: global definition of direction,
		Random direction,
		a different direction per material.
	*/
	initialDirectionTypeName = Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"initialDirectionType", "global" );
	if (0 == strcasecmp(initialDirectionTypeName, "global")) {
		initialDirectionType = INIT_DIR_GLOBAL;
	}
	else if (0 == strcasecmp(initialDirectionTypeName, "random")) {
		initialDirectionType = INIT_DIR_RANDOM;
	}
	else if (0 == strcasecmp(initialDirectionTypeName,"perMaterial")) {
		initialDirectionType = INIT_DIR_PER_MAT;
	}
	else {
			Journal_Firewall(False, Journal_Register( Error_Type, (Name)"Director"  ), 
			"Error in '%s', do not understand initialDirectionType = %s\n."
			" Options are: \"global\", \"random\" or \"perMaterial\" ", __func__, 
			initialDirectionTypeName);		
	}
	
	_Director_Init(
			self, 
			velGradField,
			materialPointsSwarm,
			initialDirectionType,
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"initialDirectionX", DIRECTOR_DEFAULT_DIR_X  ),
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"initialDirectionY", DIRECTOR_DEFAULT_DIR_Y  ),
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"initialDirectionZ", DIRECTOR_DEFAULT_DIR_Z  ),
			Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"randomInitialDirectionSeed", 1  ),
			Stg_ComponentFactory_GetBool( cf, self->name, (Dictionary_Entry_Key)"dontUpdate", False )  );

}

void _Director_Build( void* director, void* data ) {
	Director*                       self               = (Director*) director;

	/* Build parent */
	_TimeIntegrand_Build( self, data );

	Stg_Component_Build( self->directorSwarmVariable, data, False );
	Stg_Component_Build( self->dontUpdateParticle, data, False );
}

void _Director_Initialise( void* director, void* data ) {
	Director*                       self               = (Director*) director;
	Particle_Index                  lParticle_I;
	Particle_Index                  particleLocalCount = self->variable->arraySize;
	double*                         normal;
	Dimension_Index                 dim_I;

	/* Initialise Parent */
	_TimeIntegrand_Initialise( self, data );

	Stg_Component_Initialise( self->materialPointsSwarm, data, False );
	/* We should only set initial directors if in regular non-restart mode. If in restart mode, then
	the directors will be set correctly when we re-load the Swarm. */
	if ( self->context->loadFromCheckPoint == False ) {

      Stg_Component_Initialise( self->directorSwarmVariable, data, False );
      Stg_Component_Initialise( self->dontUpdateParticle, data, True );
   
      /* Update variables */
      Variable_Update( self->variable );

		particleLocalCount = self->variable->arraySize;
		if ( self->initialDirectionType == INIT_DIR_GLOBAL ) {
			for ( lParticle_I = 0 ; lParticle_I < particleLocalCount ; lParticle_I++ ) {
				/* Initialise the norm of each director */
				normal = Variable_GetPtrDouble( self->variable, lParticle_I );

				memcpy( normal, self->globalInitialDirection, sizeof(XYZ) );
			}	
		}
		else if (self->initialDirectionType == INIT_DIR_RANDOM){
			Particle_Index	gParticle_I;
			unsigned	approxGlobalParticleCount = particleLocalCount * self->materialPointsSwarm->nProc;
			unsigned	startIndex = particleLocalCount * self->materialPointsSwarm->myRank;
			double		norm[3];

			/* create random directions for each particle */
			srand( self->randomInitialDirectionSeed );

			lParticle_I = 0;
			for( gParticle_I = 0; gParticle_I < approxGlobalParticleCount; gParticle_I++ ) {
				for ( dim_I = 0; dim_I < self->materialPointsSwarm->dim; dim_I++ ) {
					norm[dim_I] = ( (float) rand() - RAND_MAX/2 ) / RAND_MAX;
				}	
				if( gParticle_I >= startIndex ) {
					normal = Variable_GetPtrDouble( self->variable, lParticle_I );
					for ( dim_I = 0; dim_I < self->materialPointsSwarm->dim; dim_I++ ) {
						normal[dim_I] = norm[dim_I];
					}
					lParticle_I++;
				}
				if( lParticle_I >= particleLocalCount )
					break;
			}
/*
			for ( lParticle_I = 0 ; lParticle_I < particleLocalCount ; lParticle_I++ ) {
				normal = Variable_GetPtrDouble( self->variable, lParticle_I );
				for ( dim_I = 0; dim_I < self->materialPointsSwarm->dim; dim_I++ ) {
					normal[dim_I] = ( (float) rand() - RAND_MAX/2 ) / RAND_MAX;
				}	
				StGermain_VectorNormalise( normal, self->materialPointsSwarm->dim );
			}
*/
		}
		else if (self->initialDirectionType == INIT_DIR_PER_MAT) {
			/* Assign initial direction based on material
			  and check first is material is defined as random.*/
			Material_Index	materialsCount = Materials_Register_GetCount( self->materialPointsSwarm->materials_Register);
			XYZ*				materialDirectionVectors;
			int				material_I;
			Material*		material;
			Bool*				randomInitialDirections;
			int*				randomInitialDirectionSeeds;
			/* int				materialOfParticle; */
			
			materialDirectionVectors = Memory_Alloc_Array(XYZ, materialsCount, "materialDirectionVectors");
			randomInitialDirectionSeeds = Memory_Alloc_Array(int, materialsCount,
			"randomInitialDirectionSeeds");
			randomInitialDirections = Memory_Alloc_Array(Bool, materialsCount,
			"randomInitialDirections");

			/* Loop over materials and get material properties from dictionary */
			for ( material_I = 0 ; material_I < materialsCount ; material_I++ ) {
				material = Materials_Register_GetByIndex( 
						self->materialPointsSwarm->materials_Register, 
						material_I );
				randomInitialDirections[material_I] = Dictionary_GetBool_WithDefault( material->dictionary, (Dictionary_Entry_Key)"randomInitialDirection", False );
				/* If no value is set, use default value from this file,
				and not the (not-used) global value */ 
				if ( randomInitialDirections[material_I] == False) {
					materialDirectionVectors[ material_I ][0] = 
						Dictionary_GetDouble_WithDefault( material->dictionary, (Dictionary_Entry_Key)"initialDirectionX", DIRECTOR_DEFAULT_DIR_X  );
				
					materialDirectionVectors[ material_I ][1] = 
						Dictionary_GetDouble_WithDefault( material->dictionary, (Dictionary_Entry_Key)"initialDirectionY", DIRECTOR_DEFAULT_DIR_Y  );				   		
				
					materialDirectionVectors[ material_I ][2] = 
						Dictionary_GetDouble_WithDefault( material->dictionary, (Dictionary_Entry_Key)"initialDirectionZ", DIRECTOR_DEFAULT_DIR_Z  );
				}
				else {
					/* If no random seed set, then use global value ( which if 
					not set will go to default value) */
					randomInitialDirectionSeeds[material_I] = Dictionary_GetUnsignedInt_WithDefault(
							material->dictionary,
							"randomInitialDirectionSeed", 
							self->randomInitialDirectionSeed);
				}
			}
			/* Assign value to particle */
			
			/* If material is random, set the local srand, 
			locate all random particles, and set their director */
			for (material_I = 0; material_I < materialsCount; material_I++) {
				if (randomInitialDirections[material_I] == True) {
					Particle_Index	gParticle_I;
					unsigned	approxGlobalParticleCount = particleLocalCount * self->materialPointsSwarm->nProc;
					unsigned	startIndex = particleLocalCount * self->materialPointsSwarm->myRank;
					double		norm[3];

					srand( self->randomInitialDirectionSeed );
					lParticle_I = 0;
					for( gParticle_I = 0; gParticle_I < approxGlobalParticleCount; gParticle_I++ ) {
						for ( dim_I = 0; dim_I < self->materialPointsSwarm->dim; dim_I++ ) {
							norm[dim_I] = ( (float) rand() - RAND_MAX/2 ) / RAND_MAX;
						}	
						if( gParticle_I >= startIndex ) {
							normal = Variable_GetPtrDouble( self->variable, lParticle_I );
							for ( dim_I = 0; dim_I < self->materialPointsSwarm->dim; dim_I++ ) {
								normal[dim_I] = norm[dim_I];
							}
							lParticle_I++;
						}
						if( lParticle_I >= particleLocalCount )
							break;
					}
#if 0
					/* create random directions for each particle */
					srand(randomInitialDirectionSeeds[material_I]);
					for ( lParticle_I = 0 ; lParticle_I < particleLocalCount ; lParticle_I++ ){
						materialOfParticle = MaterialPointsSwarm_GetMaterialIndexAt( self->materialPointsSwarm, lParticle_I );
						if (materialOfParticle == material_I){
							normal = Variable_GetPtrDouble( self->variable, lParticle_I );
							for ( dim_I = 0; dim_I < self->materialPointsSwarm->dim; dim_I++ ) {
								normal[dim_I] = ( (float) rand() - RAND_MAX/2 ) / RAND_MAX;
							}	
							/* Normalising the direction vector */
							StGermain_VectorNormalise( normal, self->materialPointsSwarm->dim );
						}
					}
#endif
				}
			}
		    /* For each non-random particle, set the initial direction */
			for ( lParticle_I = 0 ; lParticle_I < particleLocalCount ; lParticle_I++ ) {
				/* Initialise the norm of each director */
				material_I = MaterialPointsSwarm_GetMaterialIndexAt(
						self->materialPointsSwarm, 
						lParticle_I );
				if (randomInitialDirections[material_I] == False) {
					normal = Variable_GetPtrDouble(self->variable,lParticle_I);
					memcpy(normal, materialDirectionVectors[material_I], sizeof(XYZ));
					StGermain_VectorNormalise( normal, self->materialPointsSwarm->dim );			
				}
			}	
			Memory_Free(materialDirectionVectors);		
			Memory_Free(randomInitialDirections);
			Memory_Free(randomInitialDirectionSeeds);
		}
		else {
			Journal_Firewall(False, Journal_Register( Error_Type, (Name)"Director"  ), 
			"Error in '%s', do not understand self->initialDirectionType = %u\n"
			"initialDirectionType must match enumerated type, 'InitialDirectionType'"
			"in Director.h\n", __func__, 
			self->initialDirectionType);
		}
      /* Initialise the dontUpdate flag on each particle */
      for ( lParticle_I = 0 ; lParticle_I < particleLocalCount ; lParticle_I++ ) { 
         Variable_SetValueInt(self->dontUpdateParticle->variable, lParticle_I, 0);
      }	
	}

}

void _Director_Destroy( void* _self, void* data ) {
	Director*          self               = (Director*) _self;

   Stg_Component_Destroy( self->velGradField, data, False );
   Stg_Component_Destroy( self->materialPointsSwarm, data, False );
   Stg_Component_Destroy( self->directorSwarmVariable, data, False );
   Stg_Component_Destroy( self->dontUpdateParticle, data, False );
	/* Destroy parent */
	_TimeIntegrand_Destroy( self, data );

}

Bool _Director_TimeDerivative( void* director, Index lParticle_I, double* timeDeriv ) {
	Director*                self              = (Director*) director;
	MaterialPointsSwarm*     materialPointsSwarm = self->materialPointsSwarm;
	TensorArray              velGrad;
	double*                  normal;
	Element_LocalIndex       lElement_I;
	MaterialPoint*           materialPoint = (MaterialPoint*) Swarm_ParticleAt( materialPointsSwarm, lParticle_I );
	Director_ParticleExt*    particleExt;
	InterpolationResult      result;

	/* Get particle extension info - this will contain vector with director */
	particleExt = (Director_ParticleExt*) ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr,
		materialPoint, self->particleExtHandle );
	
	/* Check if this particle should be updated or not */
	if ( self->dontUpdate || particleExt->dontUpdateParticle ) {
		memset( timeDeriv, 0, sizeof(double) * materialPointsSwarm->dim );
		return True;
	}

	normal     = particleExt->director;

	lElement_I = materialPoint->owningCell;
	
	result = FieldVariable_InterpolateValueAt( self->velGradField, materialPoint->coord, velGrad );
	/* if in debug mode, perform some tests */
#ifdef DEBUG
	if ( result == OTHER_PROC || result == OUTSIDE_GLOBAL || isinf(velGrad[0]) || isinf(velGrad[1]) || 
		( self->materialPointsSwarm->dim == 3 && isinf(velGrad[2]) ) ) 
	{
		Journal_Printf( Journal_Register( Error_Type, (Name)self->type  ),
			"Error in func '%s' for particle with index %u.\n\tPosition (%g, %g, %g)\n\tVelGrad here is (%g, %g, %g)."
			"\n\tInterpolation result is %s.\n",
			__func__, lParticle_I, materialPoint->coord[0], materialPoint->coord[1], materialPoint->coord[2], 
			velGrad[0], velGrad[1], ( self->materialPointsSwarm->dim == 3 ? velGrad[2] : 0.0 ),
			InterpolationResultToStringMap[result]  );
		return False;		
	}
#endif
	if ( materialPointsSwarm->dim == 2 ) {
		timeDeriv[0] = -velGrad[MAP_2D_TENSOR(0,0)] * normal[0] - velGrad[MAP_2D_TENSOR(1,0)] * normal[1];
		timeDeriv[1] = -velGrad[MAP_2D_TENSOR(0,1)] * normal[0] - velGrad[MAP_2D_TENSOR(1,1)] * normal[1];
	}
	else {
		timeDeriv[0] = -velGrad[MAP_3D_TENSOR(0,0)] * normal[0] - velGrad[MAP_3D_TENSOR(1,0)] * normal[1] - velGrad[MAP_3D_TENSOR(2,0)] * normal[2];
		timeDeriv[1] = -velGrad[MAP_3D_TENSOR(0,1)] * normal[0] - velGrad[MAP_3D_TENSOR(1,1)] * normal[1] - velGrad[MAP_3D_TENSOR(2,1)] * normal[2];
		timeDeriv[2] = -velGrad[MAP_3D_TENSOR(0,2)] * normal[0] - velGrad[MAP_3D_TENSOR(1,2)] * normal[1] - velGrad[MAP_3D_TENSOR(2,2)] * normal[2] ;
	}

	return True;
}

/* This function is called after each of the time integration steps - 
 * here we just want to normalise the vector */
void _Director_Intermediate( void* director, Index lParticle_I ) {
	Director*          self              = (Director*) director;
	double*            normal;
	
	normal = Variable_GetPtrDouble( self->variable, lParticle_I );
	StGermain_VectorNormalise( normal, self->materialPointsSwarm->dim );
}

/* Update these variables in case the swarm has been reallocated */
void Director_UpdateVariables( void* timeIntegrator, Director* self ) {
	Variable_Update( self->variable );
	Variable_Update( self->dontUpdateParticle->variable );
	Variable_Update( self->materialPointsSwarm->particleCoordVariable->variable );
	Variable_Update( self->materialPointsSwarm->owningCellVariable->variable );
}
	
void Director_GetNormal( void* director, void* particle, XYZ normal ) {
	Director*             self              = (Director*) director;
	Director_ParticleExt* particleExt;

	particleExt = (Director_ParticleExt*) 
		ExtensionManager_Get( self->materialPointsSwarm->particleExtensionMgr, particle, self->particleExtHandle );

	memcpy( normal, particleExt->director, sizeof(XYZ) );
}

void Director_SetNormal( void* director, void* particle, XYZ normal ) {
	Director*             self              = (Director*) director;
	Director_ParticleExt* particleExt;
	
	particleExt = (Director_ParticleExt*) 
			ExtensionManager_Get( self->materialPointsSwarm->particleExtensionMgr, particle, self->particleExtHandle );
	
	memcpy( particleExt->director, normal , sizeof(XYZ) );
}

double* Director_GetNormalPtr( void* director, void* particle) {
	Director*             self              = (Director*) director;
	Director_ParticleExt* particleExt;
	
	particleExt = (Director_ParticleExt*) 
			ExtensionManager_Get( self->materialPointsSwarm->particleExtensionMgr, particle, self->particleExtHandle );
	
	return particleExt->director;
}

void Director_SetDontUpdateParticleFlag( void* director, void* particle, Particle_Bool dontUpdateFlag ){
	Director*             self              = (Director*) director;
	Director_ParticleExt* particleExt;
	
	particleExt = (Director_ParticleExt*) 
			ExtensionManager_Get( self->materialPointsSwarm->particleExtensionMgr, particle, self->particleExtHandle );
	
	particleExt->dontUpdateParticle = dontUpdateFlag;
}


