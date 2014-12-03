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
** $Id: CalculateParticleDisplacement.c 518 2007-10-11 08:07:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "CalculateParticleDisplacement.h"
#include <math.h>
#include <string.h>
#include <assert.h>

const Type PICellerator_CalculateParticleDisplacement_Type = "PICellerator_CalculateParticleDisplacement";

void _PICellerator_CalculateParticleDisplacement_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
	CalculateParticleDisplacementPlugin*  self = (CalculateParticleDisplacementPlugin*)component;
	DomainContext*                context;
	MaterialPointsSwarm*                  materialPointsSwarm;
	StandardParticle                      particle;
	ParticleDisplacementInfo*             particleExt;
	
	context = (DomainContext*)Stg_ComponentFactory_ConstructByName( cf, (Name)"context", DomainContext, True, data  );

	ContextEP_Append( context, AbstractContext_EP_Initialise, _PICellerator_CalculateParticleDisplacement_StoreOriginalPos );
	ContextEP_Append( context, AbstractContext_EP_Step, _PICellerator_CalculateParticleDisplacement_UpdateDisplacement );

	/* Extend particle with original pos */
	materialPointsSwarm = (MaterialPointsSwarm*) LiveComponentRegister_Get( context->CF->LCRegister, (Name)"materialSwarm" );
	assert( materialPointsSwarm  );
	self->particleDisplacementInfo_Handle = ExtensionManager_Add( materialPointsSwarm->particleExtensionMgr, (Name)CURR_MODULE_NAME, sizeof( ParticleDisplacementInfo )  );

	/* now register these guys as swarm variables */
	particleExt = (ParticleDisplacementInfo*)ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, &particle,
		self->particleDisplacementInfo_Handle );

	self->particleOriginalCoordSwarmVariable = Swarm_NewVectorVariable( materialPointsSwarm, (Name)"OriginalCoord", (ArithPointer) &particleExt->originalCoord - (ArithPointer) &particle,
		Variable_DataType_Double,
		materialPointsSwarm->dim,
		"OriginalCoordX",
		"OriginalCoordY",
		"OriginalCoordZ" );
		
	self->particleDisplacementSwarmVariable = Swarm_NewVectorVariable( materialPointsSwarm, (Name)"Displacement", (ArithPointer) &particleExt->displacement - (ArithPointer) &particle,
		Variable_DataType_Double,
		materialPointsSwarm->dim,
		"DisplacementX",
		"DisplacementY",
		"DisplacementZ" );

	/* Create this here, rather than in XML, since if we tried to create it in XML it couldn't
	find the Displacement swarm variable at configure-time */
	/* OperatorSwarmVariables don't automatically prefix the swarm's name to the name
	you give, so we have to add it manually here */
	self->particleDisplacementMagSwarmVariable = OperatorSwarmVariable_NewUnary(
		"materialSwarm-DisplacementMagnitude",
		(AbstractContext*)context, 
		self->particleDisplacementSwarmVariable,
		"Magnitude" );
	
	/* Need to make sure this guy is build/initialised */
}


void _PICellerator_CalculateParticleDisplacement_StoreOriginalPos( PICelleratorContext* context ) {
	CalculateParticleDisplacementPlugin*  self;
	MaterialPointsSwarm*                  materialPointsSwarm;
	GlobalParticle*                       particle;
	double*                               originalCoord;
	ParticleDisplacementInfo*             particleDisplacementInfo; 
	Particle_Index                        lParticle_I;

	materialPointsSwarm = (MaterialPointsSwarm*) LiveComponentRegister_Get( context->CF->LCRegister, (Name)"materialSwarm" );
	assert( materialPointsSwarm );
	self = (CalculateParticleDisplacementPlugin* ) LiveComponentRegister_Get( context->CF->LCRegister, (Name)"PICellerator_CalculateParticleDisplacement" );
	assert( self );

	for ( lParticle_I = 0 ; lParticle_I < materialPointsSwarm->particleLocalCount ; lParticle_I++ ) {
		particle = (GlobalParticle* )Swarm_ParticleAt( materialPointsSwarm, lParticle_I );
		particleDisplacementInfo = (ParticleDisplacementInfo*)ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr,
			particle, self->particleDisplacementInfo_Handle );
		originalCoord = particleDisplacementInfo->originalCoord;

		memcpy( originalCoord, particle->coord, sizeof(Coord) );
		memset( particleDisplacementInfo->displacement, 0, sizeof(XYZ) );
	}
}


void _PICellerator_CalculateParticleDisplacement_UpdateDisplacement( PICelleratorContext* context ) {
	CalculateParticleDisplacementPlugin*  self;
	MaterialPointsSwarm*                  materialPointsSwarm;
	GlobalParticle*                       particle;
	double*                               originalCoord;
	double*                               coord;
	Particle_Index                        lParticle_I;
	ParticleDisplacementInfo*             particleDisplacementInfo = NULL;
	Dimension_Index                       dim_I;

	/* Add original pos to particle */
	materialPointsSwarm = (MaterialPointsSwarm*) LiveComponentRegister_Get( context->CF->LCRegister, (Name)"materialSwarm" );
	assert( materialPointsSwarm );
	self = (CalculateParticleDisplacementPlugin* ) LiveComponentRegister_Get( context->CF->LCRegister, (Name)"PICellerator_CalculateParticleDisplacement" );
	assert( self );

	for ( lParticle_I = 0 ; lParticle_I < materialPointsSwarm->particleLocalCount ; lParticle_I++ ) {
		particle      = (GlobalParticle* )Swarm_ParticleAt( materialPointsSwarm, lParticle_I );
		coord         = particle->coord;
		particleDisplacementInfo = (ParticleDisplacementInfo*)ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, particle,
			self->particleDisplacementInfo_Handle );
		originalCoord = particleDisplacementInfo->originalCoord;	

		for ( dim_I = 0; dim_I < materialPointsSwarm->dim; dim_I++ ) {
			particleDisplacementInfo->displacement[dim_I] = coord[dim_I] - originalCoord[dim_I];
		}
	}
}


void* _PICellerator_CalculateParticleDisplacement_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof( CalculateParticleDisplacementPlugin );
	Type                                                      type = PICellerator_CalculateParticleDisplacement_Type;
	Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
	Stg_Class_PrintFunction*                                _print = _Codelet_Print;
	Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _PICellerator_CalculateParticleDisplacement_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _PICellerator_CalculateParticleDisplacement_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _Codelet_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _Codelet_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _Codelet_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _Codelet_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return _Codelet_New(  CODELET_PASSARGS  );
}


Index PICellerator_CalculateParticleDisplacement_Register( PluginsManager* pluginsManager ) {
	Index result;

	result = PluginsManager_Submit( pluginsManager, PICellerator_CalculateParticleDisplacement_Type, "0", _PICellerator_CalculateParticleDisplacement_DefaultNew  );

	return result;
}


