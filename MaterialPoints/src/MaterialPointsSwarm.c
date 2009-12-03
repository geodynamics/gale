/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
** Copyright (c) 2005-2006, Monash Cluster Computing, Building 28, Monash University Clayton Campus,
**	Victoria, 3800, Australia
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**
** Contributors:
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Patrick D. Sunter, Software Engineer, VPAC. (patrick@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	David Stegman, Postdoctoral Fellow, Monash University. (david.stegman@sci.monash.edu.au)
**	Wendy Sharples, PhD Student, Monash University (wendy.sharples@sci.monash.edu.au)
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Lesser General Public
**  License as published by the Free Software Foundation; either
**  version 2.1 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Lesser General Public License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License along with this library; if not, write to the Free Software
**  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
**
** $Id:  $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include <PICellerator/PopulationControl/PopulationControl.h>
#include <PICellerator/Weights/Weights.h>

#include "MaterialPoints.h"
#include "MaterialPointsSwarm.h"
#include <assert.h>
#include <string.h>
#include <math.h>


/* Textual name of this class */
const Type MaterialPointsSwarm_Type = "MaterialPointsSwarm";

MaterialPointsSwarm* MaterialPointsSwarm_New(
      Name                                  name,
      AbstractContext*                      context,
      void*                                 cellLayout,
      void*                                 particleLayout,
      Dimension_Index                       dim,
      SizeT                                 particleSize,
      FeMesh*                               mesh,
      EscapedRoutine*                       escapedRoutine, 
      Material*                             material,
      Variable_Register*                    swarmVariable_Register,
      ExtensionManager_Register*            extensionMgr_Register,
      Materials_Register*                   materials_Register,		
      MPI_Comm                              comm,
      void*                                 ics_dummy ) 
{
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(MaterialPointsSwarm);
	Type                                                      type = MaterialPointsSwarm_Type;
	Stg_Class_DeleteFunction*                              _delete = _MaterialPointsSwarm_Delete;
	Stg_Class_PrintFunction*                                _print = _MaterialPointsSwarm_Print;
	Stg_Class_CopyFunction*                                  _copy = _MaterialPointsSwarm_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _MaterialPointsSwarm_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _MaterialPointsSwarm_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _MaterialPointsSwarm_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _MaterialPointsSwarm_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _MaterialPointsSwarm_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _MaterialPointsSwarm_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;
	void*                          ics = ZERO;

	/* The following terms are parameters that have been passed into or defined in this function but are being set before being passed onto the parent */
	Particle_InCellIndex  cellParticleTblDelta = DEFAULT_CELL_PARTICLE_TBL_DELTA;
	double                extraParticlesFactor = DEFAULT_EXTRA_PARTICLES_FACTOR;

  MaterialPointsSwarm* self = _MaterialPointsSwarm_New(  MATERIALPOINTSSWARM_PASSARGS  );

   _Swarm_Init( 
         (Swarm*)self, context,
         cellLayout,
         particleLayout,
         dim,
         particleSize,
         DEFAULT_CELL_PARTICLE_TBL_DELTA,
         DEFAULT_EXTRA_PARTICLES_FACTOR,
         extensionMgr_Register,
         swarmVariable_Register,
         comm, 
         ics_dummy );

   _MaterialPointsSwarm_Init( 
      self, 
      mesh,
      escapedRoutine, 
      material,
      materials_Register );

   return self;
}


MaterialPointsSwarm* _MaterialPointsSwarm_New(  MATERIALPOINTSSWARM_DEFARGS  )
{
	MaterialPointsSwarm* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(MaterialPointsSwarm) );
	/* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
	/* This means that any values of these parameters that are passed into this function are not passed onto the parent function
	   and so should be set to ZERO in any children of this class. */
	ics = NULL;

	self = (MaterialPointsSwarm*)_Swarm_New(  SWARM_PASSARGS  );

	return self;
}


void _MaterialPointsSwarm_Init( 
		void*                                 swarm,
		FeMesh*                               mesh,
		EscapedRoutine*                       escapedRoutine, 
		Material*                             material,
		Materials_Register*                   materials_Register )
{
	MaterialPointsSwarm*    self = (MaterialPointsSwarm*)swarm;
	MaterialPoint           particle;
	GlobalParticle          globalParticle;

	self->mesh               = mesh;
	self->swarmAdvector      = NULL;		/* If we're using a SwarmAdvector, it will 'attach' itself later on. */
	self->escapedRoutine     = escapedRoutine;
	self->material           = material;
	self->materials_Register = materials_Register;
	
	self->particleCoordVariable = Swarm_NewVectorVariable(
		self,
		"Position",
		GetOffsetOfMember( globalParticle, coord ),
		Variable_DataType_Double,
		self->dim,
		"PositionX",
		"PositionY",
		"PositionZ" );

	self->materialIndexVariable = Swarm_NewScalarVariable( 
			self,
			"MaterialIndex",
			GetOffsetOfMember( particle , materialIndex ), 
			Variable_DataType_Int ); /* Should be unsigned int */

	/* If we have an escaped routine, clear the defensive flag. */
#if 0
	if( self->escapedRoutine )
		self->particleCommunicationHandler->defensive = False;
#endif
}

/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _MaterialPointsSwarm_Delete( void* swarm ) {
	MaterialPointsSwarm* self = (MaterialPointsSwarm*)swarm;

	_Swarm_Delete( self );
}


void _MaterialPointsSwarm_Print( void* swarm, Stream* stream ) {
	MaterialPointsSwarm* self = (MaterialPointsSwarm*)swarm;
	
	_Swarm_Print( self, stream );
}

void* _MaterialPointsSwarm_Copy( void* swarm, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	MaterialPointsSwarm*	self = (MaterialPointsSwarm*)swarm;
	MaterialPointsSwarm*	newMaterialPointsSwarm;
	
	newMaterialPointsSwarm = (MaterialPointsSwarm*)_Swarm_Copy( self, dest, deep, nameExt, ptrMap );

	newMaterialPointsSwarm->materialIndexVariable = (SwarmVariable*)Stg_Class_Copy(
				self->materialIndexVariable,
				NULL,
				deep,
				nameExt,
				ptrMap );
	
	return (void*)newMaterialPointsSwarm;
}

void* _MaterialPointsSwarm_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                 _sizeOfSelf = sizeof(MaterialPointsSwarm);
	Type                                                         type = MaterialPointsSwarm_Type;
	Stg_Class_DeleteFunction*                                 _delete = _MaterialPointsSwarm_Delete;
	Stg_Class_PrintFunction*                                   _print = _MaterialPointsSwarm_Print;
	Stg_Class_CopyFunction*                                     _copy = _MaterialPointsSwarm_Copy;
	Stg_Component_DefaultConstructorFunction*     _defaultConstructor = _MaterialPointsSwarm_DefaultNew;
	Stg_Component_ConstructFunction*                       _construct = _MaterialPointsSwarm_AssignFromXML;
	Stg_Component_BuildFunction*                               _build = _MaterialPointsSwarm_Build;
	Stg_Component_InitialiseFunction*                     _initialise = _MaterialPointsSwarm_Initialise;
	Stg_Component_ExecuteFunction*                           _execute = _MaterialPointsSwarm_Execute;
	Stg_Component_DestroyFunction*                           _destroy = _MaterialPointsSwarm_Destroy;
	SizeT                                                particleSize = sizeof(MaterialPoint);

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;
	void*                          ics = ZERO;

	return _MaterialPointsSwarm_New(  MATERIALPOINTSSWARM_PASSARGS  );
}


void _MaterialPointsSwarm_AssignFromXML( void* swarm, Stg_ComponentFactory* cf, void* data ) {
	MaterialPointsSwarm*	        self          = (MaterialPointsSwarm*) swarm;
	FeMesh*             		mesh;
	EscapedRoutine*                 escapedRoutine;
	Material*                       material;
	Materials_Register*             materials_Register;
	PICelleratorContext*		context;

	_Swarm_AssignFromXML( self, cf, data );

	mesh             = Stg_ComponentFactory_ConstructByKey( cf, self->name, "FeMesh", FeMesh, True, data );
	escapedRoutine   = Stg_ComponentFactory_ConstructByKey( cf, self->name, "EscapedRoutine",     EscapedRoutine,     False, data );
	material         = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Material",           Material,           False, data );

	context = (PICelleratorContext*)self->context;
	assert( Stg_CheckType( context, PICelleratorContext ) );
	materials_Register = context->materials_Register; 
	assert( materials_Register );

	_MaterialPointsSwarm_Init(
			self,
			mesh, 
			escapedRoutine,
			material,
			materials_Register );

	self->geomodHack = Dictionary_GetBool_WithDefault( cf->rootDict, "geomodHacks", False );
}

void _MaterialPointsSwarm_Build( void* swarm, void* data ) {
	MaterialPointsSwarm*	self = (MaterialPointsSwarm*) swarm;
	int			commHandler_I;
	Bool                    movementCommHandlerFound = False;
	Stream*                 errorStream = Journal_Register( Error_Type, self->type );
	int var_I;

	_Swarm_Build( self, data );

	/* Since this swarm is being set up to advect a PICellerator material, it should make sure
	 * at least one ParticleMovementHandler-type ParticleCommHandler has been added to the base
	 * Swarm. */
	for( commHandler_I=0; commHandler_I < self->commHandlerList->count; commHandler_I++ ){
		ParticleCommHandler *pComm = NULL;

		pComm = (ParticleCommHandler*)(Stg_ObjectList_At( self->commHandlerList, commHandler_I ));
		if( pComm->type == ParticleMovementHandler_Type ) {
			movementCommHandlerFound = True;
			break;
		}
	}

	Journal_Firewall( (Stg_ObjectList_Count(self->commHandlerList) >= 1) && (movementCommHandlerFound == True),
		errorStream, "Error: for MaterialPointsSwarm Swarms, at least one ParticleMovementHandler"
			" commHandler must be registered. Please rectify this in your XML / code.\n" );

	for( var_I = 0 ; var_I < self->nSwarmVars ; var_I++ ) {
		Stg_Component_Build( self->swarmVars[var_I], data , False );
	}
}
void _MaterialPointsSwarm_Initialise( void* swarm, void* data ) {
	MaterialPointsSwarm*	self 	= (MaterialPointsSwarm*) swarm;
	AbstractContext* 	context = (AbstractContext*)self->context;
	Index            	var_I	= 0;
	Particle_Index          lParticle_I=0;
	MaterialPoint*		matPoint=NULL;

	_Swarm_Initialise( self, data );

	for( var_I = 0 ; var_I < self->nSwarmVars ; var_I++ ) {
		Stg_Component_Initialise( self->swarmVars[var_I], data , False );
	}

	/* Now setup the material properties */
   if(  False == context->loadFromCheckPoint ) {

      /* Beforehand, set each particle to have UNDEFINED_MATERIAL */
      for ( lParticle_I = 0; lParticle_I < self->particleLocalCount; lParticle_I++ ) {
         matPoint = (MaterialPoint*)Swarm_ParticleAt( self, lParticle_I );
         matPoint->materialIndex = UNDEFINED_MATERIAL;
      }
		if( self->material == NULL ) {
			/* Do it by the layout of all known materials */
			Materials_Register_SetupSwarm( self->materials_Register, self );
		}
		else {
			Material_Layout( self->material, self );
			Materials_Register_AssignParticleProperties( 
					self->materials_Register, 
					self, 
					self->swarmVariable_Register->variable_Register );
		}
	}

	/** if loading from checkpoint, particle materials etc have already been loaded in Swarm_Build() - */ 
	/** possibly need to check for empty cells (and populate) if performing a interpolation restart */
   if ( True == context->loadFromCheckPoint ){
      if ( (True == self->isSwarmTypeToCheckPointAndReload) && (True == context->interpolateRestart) ) {	   
         Particle_InCellIndex cParticle_I         = 0;
         Particle_InCellIndex particle_I          = 0;
         GlobalParticle*      particle            = NULL;
         double               minDistance         = HUGE_VAL;
         double               distanceToParticle;
         Dimension_Index      dim = self->dim;
         Cell_DomainIndex     dCell_I;
         Cell_LocalIndex      lCell_I;
         Cell_DomainIndex     belongsToCell_I = 0;
         MaterialPoint*       matNewParticle;
         MaterialPoint*       matParticleToSplit;
         Particle_Index       matNewParticle_IndexOnCPU;
         Coord                xi;
         Coord                coord;
         unsigned             nEmptyCells = 0;
         Index*				   cellID = NULL;
         Index*				   particleCPUID = NULL;
         unsigned             count;
         unsigned             ii;
         
         /** first determine how many local cells are empty */
         for( lCell_I = 0 ; lCell_I < self->cellLocalCount ; lCell_I++ )
            if (self->cellParticleCountTbl[lCell_I] == 0) nEmptyCells++;

         /** create arrays which will be later used to populate cells */
         cellID        = Memory_Alloc_Array( Index, nEmptyCells, "Cell ID for cell to be populated" );
         particleCPUID = Memory_Alloc_Array( Index, nEmptyCells, "particle ID for particle to populate cell" );

         count = 0;
         for( lCell_I = 0 ; lCell_I < self->cellLocalCount ; lCell_I++ ) {
            minDistance = HUGE_VAL;
            if (self->cellParticleCountTbl[lCell_I] == 0) {
               /** select the centre of the current cell */
               xi[0] = 0;  xi[1] = 0;  xi[2] = 0;
               /** get global coord */
               FeMesh_CoordLocalToGlobal( self->mesh, lCell_I, xi, coord );
      
               for( dCell_I = 0 ; dCell_I < self->cellDomainCount ; dCell_I++ ) {
                  /** Loop over particles find closest to cell centre */
                  for( cParticle_I = 0 ; cParticle_I < self->cellParticleCountTbl[dCell_I] ; cParticle_I++ ) {
                     particle = (GlobalParticle*)Swarm_ParticleInCellAt( self, dCell_I, cParticle_I );
               
                     /** Calculate distance to particle */
                     distanceToParticle = 
                        (particle->coord[ I_AXIS ] - coord[ I_AXIS ]) * 
                        (particle->coord[ I_AXIS ] - coord[ I_AXIS ]) +
                        (particle->coord[ J_AXIS ] - coord[ J_AXIS ]) * 
                        (particle->coord[ J_AXIS ] - coord[ J_AXIS ]) ;
               
                     if (dim == 3) {
                        distanceToParticle += 
                           (particle->coord[ K_AXIS ] - coord[ K_AXIS ]) * 
                           (particle->coord[ K_AXIS ] - coord[ K_AXIS ]) ;
                     }
                     /** Don't do square root here because it is unnessesary: i.e. a < b <=> sqrt(a) < sqrt(b) */
                        
                     /** Check if this is the closest particle */
                     if (minDistance > distanceToParticle) {
                        particle_I       = cParticle_I;
                        minDistance      = distanceToParticle;
                        belongsToCell_I  = dCell_I;
                     }
                  }
               }

               /** create new particle which will be placed at centre of empty cell */
               matNewParticle      = (MaterialPoint*) Swarm_CreateNewParticle( self, &matNewParticle_IndexOnCPU );
               /** grab closest particle, which we will copy */
               matParticleToSplit  = (MaterialPoint*) Swarm_ParticleInCellAt( self, belongsToCell_I, particle_I );
               /** copy */
               memcpy( matNewParticle, matParticleToSplit, self->particleExtensionMgr->finalSize );
               /** set the owningCell to the cellDomainCount so that its owningCell will be reinitialised (not sure if necessary) */ 
               matNewParticle->owningCell = self->cellDomainCount;
               /** Copy new global position (cell centre) to coord on new mat particle */
               memcpy( matNewParticle->coord, coord, sizeof(Coord) );
               
               /** we now store the required information to populate empty cells */
               /** note that cells are not populated at this point, as this may interfere
                   with the 0th order interpolation we are performing */ 
               cellID[count]        = lCell_I;
               particleCPUID[count] = matNewParticle_IndexOnCPU;
               count++;
      
            }
         }
         /** populate empty cells */ 
         for(ii = 0 ; ii < count ; ii++)
            Swarm_AddParticleToCell( self, cellID[ii], particleCPUID[ii] );
         Memory_Free( cellID );
         Memory_Free( particleCPUID );
      }
	/* TODO: print info / debug message */
	}
#if 0
	else {
		Particle_Index          lParticle_I=0;
		MaterialPoint*		matPoint=NULL;

		/* Beforehand, set each particle to have UNDEFINED_MATERIAL */
		for ( lParticle_I = 0; lParticle_I < self->particleLocalCount; lParticle_I++ ) {
			matPoint = (MaterialPoint*)Swarm_ParticleAt( self, lParticle_I );
			matPoint->materialIndex = UNDEFINED_MATERIAL;
		}
		/* Now setup the material properties */
		if ( self->material == NULL ) {
			/* Do it by the layout of all known materials */
			Materials_Register_SetupSwarm( self->materials_Register, self );
		}
		else {
			Material_Layout( self->material, self );
			Materials_Register_AssignParticleProperties( 
					self->materials_Register, 
					self, 
					self->swarmVariable_Register->variable_Register );
		}
	}
#endif
}
void _MaterialPointsSwarm_Execute( void* swarm, void* data ) {
	MaterialPointsSwarm*	self = (MaterialPointsSwarm*)swarm;
	
	_Swarm_Execute( self, data );
}
void _MaterialPointsSwarm_Destroy( void* swarm, void* data ) {
	MaterialPointsSwarm*	self = (MaterialPointsSwarm*)swarm;
	
	_Swarm_Destroy( self, data );
}

void _MaterialPointsSwarm_UpdateHook( void* timeIntegrator, void* swarm ) {
	MaterialPointsSwarm* self               = (MaterialPointsSwarm*)swarm;
	FeMesh*  mesh               = self->mesh;
	Index                cell;
	Index                point_I;
	MaterialPoint*       materialPoint;

	/* Need to check for escaped particles before the next block. */
	if ( self->escapedRoutine ) {
		Stg_Component_Execute( self->escapedRoutine, self, True );
	}

	/* Check that particles have not exited the box after advection */
	if ( self->swarmAdvector  ) {
		for ( point_I = 0; point_I < self->particleLocalCount; ++point_I ) {
			materialPoint = (MaterialPoint*)Swarm_ParticleAt( self, point_I );
			cell = materialPoint->owningCell;
			Journal_Firewall(
					 cell < FeMesh_GetElementDomainSize( mesh ), 
				Journal_MyStream( Error_Type, self ),
				"In func %s: MaterialPoint '%d' outside element. Coord = {%g, %g, %g}\n",
				__func__,
				point_I,
				materialPoint->coord[ I_AXIS ],
				materialPoint->coord[ J_AXIS ],
				materialPoint->coord[ K_AXIS ] );
		}
	}

}

void MaterialPointsSwarm_SetMaterialAt( void* swarm, Index point_I, Index materialIndex ) {
	MaterialPointsSwarm* self  = (MaterialPointsSwarm*)swarm;
	MaterialPoint*       point;

	point = (MaterialPoint*)Swarm_ParticleAt( self, point_I );
	point->materialIndex = materialIndex;
}


Material* MaterialPointsSwarm_GetMaterialOn( void* swarm, void* particle ) {
	MaterialPointsSwarm* self  = (MaterialPointsSwarm*)swarm;
	MaterialPoint*       materialPoint = (MaterialPoint*)particle;

	return Materials_Register_GetByIndex( self->materials_Register, materialPoint->materialIndex );
}


Material* MaterialPointsSwarm_GetMaterialAt( void* swarm, Index point_I ) {
	MaterialPointsSwarm* self  = (MaterialPointsSwarm*)swarm;
	MaterialPoint*       point;

	point = (MaterialPoint*)Swarm_ParticleAt( self, point_I );
	return Materials_Register_GetByIndex( self->materials_Register, point->materialIndex );
}

Index MaterialPointsSwarm_GetMaterialIndexAt( void* swarm, Index point_I ) {
	MaterialPointsSwarm* self  = (MaterialPointsSwarm*)swarm;
	MaterialPoint*       point;

	point = (MaterialPoint*)Swarm_ParticleAt( self, point_I );
	return point->materialIndex;
}

void* MaterialPointsSwarm_GetMaterialExtensionOn( void* swarm, void *matPoint, Index extHandle ) {
	MaterialPointsSwarm *self  = (MaterialPointsSwarm*)swarm;
	Material *mat;

	mat = MaterialPointsSwarm_GetMaterialOn(self, matPoint);
	return ExtensionManager_Get(mat->extensionMgr, mat, extHandle);
}

void* MaterialPointsSwarm_GetMaterialExtensionAt( void* swarm, int matPointInd, Index extHandle ) {
	MaterialPointsSwarm *self  = (MaterialPointsSwarm*)swarm;
	Material *mat;

	mat = MaterialPointsSwarm_GetMaterialAt(self, matPointInd);
	return ExtensionManager_Get(mat->extensionMgr, mat, extHandle);
}

void* MaterialPointsSwarm_GetExtensionAt( void* swarm, Index point_I, Index extHandle ) {
	MaterialPointsSwarm* self  = (MaterialPointsSwarm*)swarm;
	MaterialPoint*       point;

	point = (MaterialPoint*)Swarm_ParticleAt( self, point_I );
	return ExtensionManager_Get( self->particleExtensionMgr, point, extHandle );
	
}


