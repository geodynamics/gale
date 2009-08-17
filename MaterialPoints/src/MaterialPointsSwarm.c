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
		MPI_Comm                              comm) 
{
	MaterialPointsSwarm* self = (MaterialPointsSwarm*) _MaterialPointsSwarm_DefaultNew( name );
	
	self->particleSize = particleSize;

	/* 	MaterialPointsSwarm_InitAll */
	_Swarm_Init(  	(Swarm*)self, 
			cellLayout, 
			particleLayout, 
			dim, 
			DEFAULT_CELL_PARTICLE_TBL_DELTA,
			DEFAULT_EXTRA_PARTICLES_FACTOR,
			extensionMgr_Register,
			swarmVariable_Register, 
			comm,
			NULL );

	_MaterialPointsSwarm_Init(	self, 
					mesh, 
					escapedRoutine, 
					material, 
					materials_Register);

	return self;
}


MaterialPointsSwarm* _MaterialPointsSwarm_New(
		SizeT                                           _sizeOfSelf,
		Type                                            type,
		Stg_Class_DeleteFunction*                       _delete,
		Stg_Class_PrintFunction*                        _print,
		Stg_Class_CopyFunction*                         _copy,
		Stg_Component_DefaultConstructorFunction*       _defaultConstructor,
		Stg_Component_ConstructFunction*                _construct,
		Stg_Component_BuildFunction*                    _build,
		Stg_Component_InitialiseFunction*               _initialise,
		Stg_Component_ExecuteFunction*                  _execute,
		Stg_Component_DestroyFunction*                  _destroy,
		Name                                            name,
		Bool                                            initFlag,
		CellLayout*                                     cellLayout,
		ParticleLayout*                                 particleLayout,
		Dimension_Index                                 dim,
		SizeT                                           particleSize,
		Particle_InCellIndex                            cellParticleTblDelta,
		double                                          extraParticlesFactor,
		FeMesh*                   	                mesh,
		EscapedRoutine*                                 escapedRoutine, 
		Material*                                       material,
		Variable_Register*                              swarmVariable_Register,
		ExtensionManager_Register*                      extensionMgr_Register,
		Materials_Register*                             materials_Register,
		MPI_Comm                                        comm )
{
	MaterialPointsSwarm* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(MaterialPointsSwarm) );
	self = (MaterialPointsSwarm*)_Swarm_New( 
			_sizeOfSelf,
			type,
			_delete,
			_print,
			_copy,
			_defaultConstructor,
			_construct,
			_build,
			_initialise,
			_execute,
			_destroy,		
			name,
			initFlag,
			cellLayout,
			particleLayout,
			dim,
			particleSize,
			cellParticleTblDelta,
			extraParticlesFactor,
			extensionMgr_Register,
			swarmVariable_Register,
			comm,
		        NULL	);

	if (initFlag) {
	 	_MaterialPointsSwarm_Init( 
			self,
			mesh,
			escapedRoutine, 
			material,
			materials_Register );
	}	

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
	return _MaterialPointsSwarm_New(
			sizeof(MaterialPointsSwarm),
			MaterialPointsSwarm_Type,
			_MaterialPointsSwarm_Delete,
			_MaterialPointsSwarm_Print,
			_MaterialPointsSwarm_Copy,
			_MaterialPointsSwarm_DefaultNew,
			_MaterialPointsSwarm_Construct,
			_MaterialPointsSwarm_Build,
			_MaterialPointsSwarm_Initialise,
			_MaterialPointsSwarm_Execute,
			_MaterialPointsSwarm_Destroy,
			name,
			False,
			NULL,			/* cellLayout */
			NULL,                   /* particleLayout */
			0,                      /* dim */
			sizeof(MaterialPoint),  /* particleSize */
			0,                      /* cellParticleTblDelta */
			0.0,                    /* extraParticlesFactor */
			NULL,                   /* mesh */
			NULL,                   /* escapedRoutine */
			NULL,                   /* material */
			NULL,                   /* swarmVariable_Register */
			NULL,                   /* extensionMgr_Register */
			NULL,                   /* materials_Register */
			0                       /* comm */
			);
}


void _MaterialPointsSwarm_Construct( void* swarm, Stg_ComponentFactory* cf, void* data ) {
	MaterialPointsSwarm*	        self          = (MaterialPointsSwarm*) swarm;
	FeMesh*             mesh;
	EscapedRoutine*                 escapedRoutine;
	Material*                       material;
	Materials_Register*             materials_Register;

	_Swarm_Construct( self, cf, data );

	mesh             = Stg_ComponentFactory_ConstructByKey( cf, self->name, "FeMesh", FeMesh, True, data );
	escapedRoutine   = Stg_ComponentFactory_ConstructByKey( cf, self->name, "EscapedRoutine",     EscapedRoutine,     False, data );
	material         = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Material",           Material,           False, data );

	materials_Register = Stg_ObjectList_Get( cf->registerRegister, "Materials_Register" );
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
	MaterialPointsSwarm*	self = (MaterialPointsSwarm*) swarm;
	AbstractContext* context = (AbstractContext*)data;
	int var_I;
	
	_Swarm_Initialise( self, data );
	for( var_I = 0 ; var_I < self->nSwarmVars ; var_I++ ) {
		Stg_Component_Initialise( self->swarmVars[var_I], data , False );
	}

	/* if loading from checkpoint, particle materials etc have already been loaded in Swarm_Build() - 
	 * thus nothing to do here */
	if ( context && (True == context->loadFromCheckPoint) ) {
		/* TODO: print info / debug message */
	}
	/* else for a normal run, lay out the particle->material mappings based on user defined materials */
	else {
		/* Setup the material properties */
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

void* MaterialPointsSwarm_GetExtensionAt( void* swarm, Index point_I, Index extHandle ) {
	MaterialPointsSwarm* self  = (MaterialPointsSwarm*)swarm;
	MaterialPoint*       point;

	point = (MaterialPoint*)Swarm_ParticleAt( self, point_I );
	return ExtensionManager_Get( self->particleExtensionMgr, point, extHandle );
	
}
