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
** $Id: IntegrationPointsSwarm.c 518 2007-10-11 08:07:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include <PICellerator/PopulationControl/PopulationControl.h>
#include <PICellerator/Weights/Weights.h>

#include "MaterialPoints.h"
#include "Material.h"
#include "PICelleratorContext.h"

#include <assert.h>
#include <string.h>
#include <math.h>


/* Textual name of this class */
const Type IntegrationPointsSwarm_Type = "IntegrationPointsSwarm";

IntegrationPointsSwarm* IntegrationPointsSwarm_New(
		Name                                  name,
      AbstractContext*                      context,
		void*                                 cellLayout,
		void*                                 particleLayout,
		Dimension_Index                       dim,
		SizeT                                 particleSize,
		Particle_InCellIndex                  cellParticleTblDelta,
		double                                extraParticlesFactor,
		FeMesh*                   	      mesh,
		TimeIntegrator*                       timeIntegrator,
		WeightsCalculator*                    weights,
		IntegrationPointMapper*               mapper,
		Bool                                  recalculateWeights,
		ExtensionManager_Register*            extensionMgr_Register,
		Variable_Register*                    swarmVariable_Register,
		Materials_Register*                   materials_Register,
		MPI_Comm                              comm,
      void*                                 ics_dummy)
{
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(IntegrationPointsSwarm);
	Type                                                      type = IntegrationPointsSwarm_Type;
	Stg_Class_DeleteFunction*                              _delete = _IntegrationPointsSwarm_Delete;
	Stg_Class_PrintFunction*                                _print = _IntegrationPointsSwarm_Print;
	Stg_Class_CopyFunction*                                  _copy = _IntegrationPointsSwarm_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _IntegrationPointsSwarm_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _IntegrationPointsSwarm_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _IntegrationPointsSwarm_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _IntegrationPointsSwarm_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _IntegrationPointsSwarm_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _IntegrationPointsSwarm_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;
	void*                          ics = ZERO;

    IntegrationPointsSwarm* self = _IntegrationPointsSwarm_New(  INTEGRATIONPOINTSSWARM_PASSARGS  );

    _Swarm_Init( 
         (Swarm*)self, context,
         cellLayout,
         particleLayout,
         dim,
         particleSize,
         cellParticleTblDelta,
         extraParticlesFactor,
         extensionMgr_Register,
         swarmVariable_Register,
         comm, 
         ics_dummy );
   _IntegrationPointsSwarm_Init( 
      self,
      mesh, 
      timeIntegrator,
      weights,
      mapper,
      materials_Register,
      recalculateWeights );
}


void* _IntegrationPointsSwarm_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                 _sizeOfSelf = sizeof(IntegrationPointsSwarm);
	Type                                                         type = IntegrationPointsSwarm_Type;
	Stg_Class_DeleteFunction*                                 _delete = _IntegrationPointsSwarm_Delete;
	Stg_Class_PrintFunction*                                   _print = _IntegrationPointsSwarm_Print;
	Stg_Class_CopyFunction*                                     _copy = _IntegrationPointsSwarm_Copy;
	Stg_Component_DefaultConstructorFunction*     _defaultConstructor = _IntegrationPointsSwarm_DefaultNew;
	Stg_Component_ConstructFunction*                       _construct = _IntegrationPointsSwarm_AssignFromXML;
	Stg_Component_BuildFunction*                               _build = _IntegrationPointsSwarm_Build;
	Stg_Component_InitialiseFunction*                     _initialise = _IntegrationPointsSwarm_Initialise;
	Stg_Component_ExecuteFunction*                           _execute = _IntegrationPointsSwarm_Execute;
	Stg_Component_DestroyFunction*                           _destroy = _IntegrationPointsSwarm_Destroy;
	SizeT                                                particleSize = sizeof(IntegrationPoint);

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;
	void*                          ics = ZERO;

	return (void*) _IntegrationPointsSwarm_New(  INTEGRATIONPOINTSSWARM_PASSARGS  );
}


IntegrationPointsSwarm* _IntegrationPointsSwarm_New(  INTEGRATIONPOINTSSWARM_DEFARGS  )
{
	IntegrationPointsSwarm* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(IntegrationPointsSwarm) );
	/* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
	/* This means that any values of these parameters that are passed into this function are not passed onto the parent function
	   and so should be set to ZERO in any children of this class. */
	ics = NULL;

	self = (IntegrationPointsSwarm*)_Swarm_New(  SWARM_PASSARGS  );

	return self;
}


void _IntegrationPointsSwarm_AssignFromXML( void* integrationPoints, Stg_ComponentFactory* cf, void* data ) {
	IntegrationPointsSwarm*	        self          = (IntegrationPointsSwarm*) integrationPoints;
	FeMesh*             		mesh;
	TimeIntegrator*                 timeIntegrator;
	WeightsCalculator*              weights;
	IntegrationPointMapper*         mapper;
	Materials_Register*             materials_Register;
	Bool                            recalculateWeights;
	PICelleratorContext*		context;

	/* This will also call _Swarm_Init */
	_Swarm_AssignFromXML( self, cf, data );

	mesh           = Stg_ComponentFactory_ConstructByKey( cf, self->name, "FeMesh", FeMesh, True, data );
	timeIntegrator = Stg_ComponentFactory_ConstructByKey( cf, self->name, "TimeIntegrator", TimeIntegrator, True, data );
	weights        = Stg_ComponentFactory_ConstructByKey( cf, self->name, "WeightsCalculator", WeightsCalculator, False, data );
	mapper         = Stg_ComponentFactory_ConstructByKey( cf, self->name, "IntegrationPointMapper", IntegrationPointMapper, True, data );
	recalculateWeights = Stg_ComponentFactory_GetBool( cf, self->name, "recalculateWeights", True );

	Journal_Firewall (
			weights != NULL ||
			(weights == NULL && (Stg_Class_IsInstance( mapper, GaussMapper_Type ) ||
					     Stg_Class_IsInstance( mapper, GaussCoincidentMapper_Type) ||
					     !strcmp( mapper->type, "PCDVCGaussMapper"))),
			Journal_MyStream( Error_Type, self ),
			"In func %s, %s which is a %s must either have a %s or use %s\n",
			__func__,
			self->name,
			self->type,
			WeightsCalculator_Type,
			GaussMapper_Type );
	
	context = (PICelleratorContext*)self->context;
	assert( Stg_CheckType( context, PICelleratorContext ) );
	materials_Register = context->materials_Register;
	assert( materials_Register );

	_IntegrationPointsSwarm_Init( self, mesh, timeIntegrator, weights, mapper, materials_Register, recalculateWeights );
}


void _IntegrationPointsSwarm_Init( 
		void*                                 swarm,
		FeMesh*                   mesh, 
		TimeIntegrator*                       timeIntegrator,
		WeightsCalculator*                    weights,
		IntegrationPointMapper*               mapper,
		Materials_Register*                   materials_Register,
		Bool                                  recalculateWeights )
{
	IntegrationPointsSwarm* self = (IntegrationPointsSwarm*)swarm;
	LocalParticle localParticle;
	IntegrationPoint particle;

	self->mesh               = mesh;
	self->timeIntegrator     = timeIntegrator;
	self->weights            = weights;
	self->mapper             = mapper;
	self->materials_Register = materials_Register;

	self->recalculateWeights = recalculateWeights;

	/* Disable checkpointing and reloading of IP swarms - currently they can't be reloaded if the particles
	don't have a global coord. We assume there is no history info on them which means we're happy to re-create
	them from scratch given the position the material points were in when the checkpoint was made as input
	-- PatrickSunter 12 June 2006 */
	self->isSwarmTypeToCheckPointAndReload = False;

	self->weightVariable = Swarm_NewScalarVariable( 
			self,
			"Weight",
			GetOffsetOfMember( particle , weight ), 
			Variable_DataType_Double );

	self->localCoordVariable = Swarm_NewVectorVariable(
		self,
		"LocalElCoord",
		GetOffsetOfMember( localParticle , xi ),
		Variable_DataType_Double,
		self->dim,
		"Xi",
		"Eta",
		"Zeta" );

	if ( timeIntegrator ) {
		/* Assuming this is called from _IntegrationPointsSwarm_AssignFromXML, it would have always called construct
		 * on the mapper which in turn would have constructed any MaterialPointsSwarms with it.
		 * The MaterialPointsSwarms would have already appended their update routines to the EP, and hence this
		 * ensures that the _IntegrationPointsSwarm_UpdateHook will always be called last */
		TimeIntegrator_InsertAfterFinishEP(
				timeIntegrator,
				"MaterialPointsSwarm_Update", /* Needs to be after a the material update */
				"IntegrationPointsSwarm_Update",
				_IntegrationPointsSwarm_UpdateHook,
				self->name,
				self );
	}
	
	/* _Construct calls _Swarm_Init */

	/* Lock down the extension manager.
	 * It doesn't make sense for the IntegrationPointsSwarm to allow IntegrationPoints to be extended
	 * This means attempts to extend integration points are firewalled to pickup errors.
	 * -- Alan 20060506
	 */
	ExtensionManager_SetLockDown( self->particleExtensionMgr, True );
}

void _IntegrationPointsSwarm_Delete( void* integrationPoints ) {
	IntegrationPointsSwarm* self = (IntegrationPointsSwarm*)integrationPoints;

	_Swarm_Delete( self );
}


void _IntegrationPointsSwarm_Print( void* integrationPoints, Stream* stream ) {
	IntegrationPointsSwarm* self = (IntegrationPointsSwarm*)integrationPoints;
	
	_Swarm_Print( self, stream );
}

void* _IntegrationPointsSwarm_Copy( void* integrationPoints, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	IntegrationPointsSwarm*	self = (IntegrationPointsSwarm*)integrationPoints;
	IntegrationPointsSwarm*	newIntegrationPointsSwarm;
	
	newIntegrationPointsSwarm = (IntegrationPointsSwarm*)_Swarm_Copy( self, dest, deep, nameExt, ptrMap );

	newIntegrationPointsSwarm->mesh = (FeMesh*)Stg_Class_Copy(
				self->mesh,
				NULL,
				deep,
				nameExt,
				ptrMap );
	if ( self->weights != NULL ) {
		newIntegrationPointsSwarm->weights = (WeightsCalculator*)Stg_Class_Copy(
					self->weights,
					NULL,
					deep,
					nameExt,
					ptrMap );
	}
	newIntegrationPointsSwarm->mapper = (IntegrationPointMapper*)Stg_Class_Copy(
				self->mapper,
				NULL,
				deep,
				nameExt,
				ptrMap );
	newIntegrationPointsSwarm->materials_Register = self->materials_Register;
	newIntegrationPointsSwarm->weightVariable = (SwarmVariable*)Stg_Class_Copy(
				self->weightVariable,
				NULL,
				deep,
				nameExt,
				ptrMap );
	
	return (void*)newIntegrationPointsSwarm;
}


void _IntegrationPointsSwarm_Build( void* integrationPoints, void* data ) {
	IntegrationPointsSwarm*	self = (IntegrationPointsSwarm*) integrationPoints;

	_Swarm_Build( self, data );

	Stg_Component_Build( self->localCoordVariable, data, False );
	Stg_Component_Build( self->weightVariable, data, False );
	Stg_Component_Build( self->mapper, data, False );
   Stg_Component_Build( self->mesh, data, False );
   Stg_Component_Build( self->timeIntegrator, data, False );
   if ( self->weights != NULL ) {
		Stg_Component_Build( self->weights, data, False );
	}
   
}
void _IntegrationPointsSwarm_Initialise( void* integrationPoints, void* data ) {
	IntegrationPointsSwarm*	self = (IntegrationPointsSwarm*) integrationPoints;

	Journal_DPrintf( self->debug, "In %s(): for swarm \"%s\":\n",
		__func__, self->name );
	Stream_IndentBranch( Swarm_Debug );

	_Swarm_Initialise( self, data );

	Stg_Component_Initialise( self->localCoordVariable, data, False );
	Stg_Component_Initialise( self->weightVariable, data, False );
	Stg_Component_Initialise( self->mapper, data, False );
   Stg_Component_Initialise( self->mesh, data, False );
   Stg_Component_Initialise( self->timeIntegrator, data, False );

	if ( self->weights != NULL ) {
		Stg_Component_Initialise( self->weights, data, False );
	}

	/* We call this function to actually set up the integration point positions and 
	weights, based on the now set up material point swarm */
	IntegrationPointsSwarm_RemapIntegrationPointsAndRecalculateWeights( self );

	Stream_UnIndentBranch( Swarm_Debug );
	Journal_DPrintf( self->debug, "...done in %s() for swarm \"%s\".\n",
		__func__, self->name );
}
void _IntegrationPointsSwarm_Execute( void* integrationPoints, void* data ) {
	IntegrationPointsSwarm*	self = (IntegrationPointsSwarm*)integrationPoints;
	
	_Swarm_Execute( self, data );
}
void _IntegrationPointsSwarm_Destroy( void* integrationPoints, void* data ) {
	IntegrationPointsSwarm*	self = (IntegrationPointsSwarm*)integrationPoints;

	Stg_Component_Destroy( self->localCoordVariable, data, False );
	Stg_Component_Destroy( self->weightVariable, data, False );
	Stg_Component_Destroy( self->mapper, data, False );
   Stg_Component_Destroy( self->mesh, data, False );
   Stg_Component_Destroy( self->timeIntegrator, data, False );
   if ( self->weights != NULL ) {
		Stg_Component_Destroy( self->weights, data, False );
	}
	
	_Swarm_Destroy( self, data );
}

void _IntegrationPointsSwarm_UpdateHook( void* timeIntegrator, void* swarm ) {
	IntegrationPointsSwarm* self = (IntegrationPointsSwarm*)swarm;

	Journal_DPrintf( self->debug, "In %s(): for swarm \"%s\":\n",
		__func__, self->name );
	Stream_IndentBranch( Swarm_Debug );

	IntegrationPointsSwarm_RemapIntegrationPointsAndRecalculateWeights( self );

	Stream_UnIndentBranch( Swarm_Debug );
	Journal_DPrintf( self->debug, "...done in %s() for swarm \"%s\".\n",
		__func__, self->name );
}


void IntegrationPointsSwarm_RemapIntegrationPointsAndRecalculateWeights( void* swarm ) {	
	IntegrationPointsSwarm* self = (IntegrationPointsSwarm*)swarm;
	double                  mapStartTime, mapTime;
	double                  weightsUpdateStartTime, weightsUpdateTime;

	Journal_DPrintf( self->debug, "In %s(): for swarm \"%s\":\n",
		__func__, self->name );
	Stream_IndentBranch( Swarm_Debug );

	Journal_DPrintf( self->debug, "Calling IntegrationPointsMapper \"%s\" (of type %s) to set up mappings\n"
		"\tfrom I.P.s to M.P.s, and calculate local coords:\n", self->mapper->name, self->mapper->type );
	mapStartTime = MPI_Wtime();	
	IntegrationPointMapper_Map( self->mapper );
	mapTime = MPI_Wtime() - mapStartTime;
	Journal_DPrintf( self->debug, "...done - took %g secs.\n", mapTime );

	if ( self->weights != NULL ) {
		Journal_DPrintf( self->debug, "Calling WeightsCalculator \"%s\" (of type %s)\n"
			"\tto calculate and set integration weights:\n",
			self->weights->name, self->weights->type );
		weightsUpdateStartTime = MPI_Wtime();
		WeightsCalculator_CalculateAll(self->weights, self );
		weightsUpdateTime = MPI_Wtime() - weightsUpdateStartTime;
		Journal_DPrintf( self->debug, "...weights updating finished - took %g secs.\n", weightsUpdateTime );
	}	
	else {
		Stream* errorStream = Journal_Register( Error_Type, self->type );
		Journal_Firewall( Stg_Class_IsInstance( self->mapper, GaussMapper_Type ) ||
				  Stg_Class_IsInstance( self->mapper, GaussCoincidentMapper_Type ) ||
				  !strcmp(self->mapper->type, "PCDVCGaussMapper"), errorStream,
			"Error - in %s(): for IntegrationPointSwarm \"%s\", no weights calculator provided "
			"and mapper is not a %s.\n", GaussMapper_Type );

		Journal_DPrintf( self->debug, "not recalculating weights since we are using a %s mapper and "
			"assume the points are not being advected.\n", GaussMapper_Type );
	}

	Stream_UnIndentBranch( Swarm_Debug );
	Journal_DPrintf( self->debug, "...done in %s() for swarm \"%s\"\n",
		__func__, self->name );
}


Material_Index IntegrationPointsSwarm_GetMaterialIndexOn( IntegrationPointsSwarm* swarm, IntegrationPoint* point ) {
	return IntegrationPointMapper_GetMaterialIndexOn( swarm->mapper, point );
}


Material* IntegrationPointsSwarm_GetMaterialOn( IntegrationPointsSwarm* swarm, IntegrationPoint* point ) {
	return Materials_Register_GetByIndex(
		swarm->materials_Register, 
		IntegrationPointsSwarm_GetMaterialIndexOn( swarm, point ) );
}


