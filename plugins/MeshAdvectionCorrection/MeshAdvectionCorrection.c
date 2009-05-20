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
*+
** $Id: MeshAdvectionCorrection.c 610 2007-10-11 08:09:29Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include <string.h>

/* Each component type needs a unique identifier (as a string for us to read and as an integer for realtime comparisions) */
const Type Underworld_MeshAdvectionCorrection_Type = "Underworld_MeshAdvectionCorrection_Type";

ExtensionInfo_Index Underworld_MeshAdvectionCorrection_ContextExtHandle;

typedef struct {
	/* Save the "old" execute function pointer (which we replace) */ 
	Stg_Component_ExecuteFunction* energySolverExecute;
	FeVariable* artDisplacement;
	double dt;
} Underworld_MeshAdvectionCorrection_ContextExt;


typedef enum {
	ADD,
	MINUS
} CorrectionFlag;

void MeshAdvectionCorrection_AddCorrection( FeVariable* velocityField, double* oldVelocity, double* artVelocity, CorrectionFlag flag ) {
	/* This function offsets the VelocityField depending on the the flag:
	 * if flag == ADD
	 * 		velocity = oldVelocity
	 * else {
	 * 		if( artificalVelocity ) {	velocity = artVelocity }
	 * 		else { velocity = (0,0) } 
	 */
	FeVariable*  self          = (FeVariable*)  velocityField;
	FeMesh*      mesh          = self->feMesh;
	int          lNodeCount    = FeMesh_GetNodeLocalSize( mesh );
	int          dof           = self->fieldComponentCount;
	int          dim           = self->dim;
	double       oldV[3], artV[3], zero[3];
	int lNode_I;

	zero[0]=zero[1]=zero[2] = 0;

	for( lNode_I = 0 ; lNode_I < lNodeCount ; lNode_I++ ) {

		if ( flag == ADD ) {
			/* velocity = stored values from oldVelocity */
			memcpy( oldV, &oldVelocity[lNode_I*dof], dof*sizeof(double) );
			FeVariable_SetValueAtNode( velocityField, lNode_I, oldV );
		}
		else {
			if( artVelocity != 0 ) { /* use artificial velocity */
				memcpy( artV, &artVelocity[lNode_I*dof], dof*sizeof(double) );
				FeVariable_SetValueAtNode( velocityField, lNode_I, artV );
			} else { /* just zero the velocity */
				FeVariable_SetValueAtNode( velocityField, lNode_I, zero );
			}
		}

	}
}

MeshAdvectionCorrection_StoreCurrentVelocity( FeVariable *velocityField, double *oldVelocity ) {
	/* save the current values of the velocity field in the oldVelocity array */
	FeMesh *mesh = velocityField->feMesh;
	unsigned dof = velocityField->fieldComponentCount;
	unsigned numLocalNodes = FeMesh_GetNodeLocalSize( mesh );
	unsigned lNode_I;
	double vel[3];

	for( lNode_I = 0 ; lNode_I < numLocalNodes ; lNode_I++ ) {
		FeVariable_GetValueAtNode( velocityField, lNode_I, vel );
		memcpy( &oldVelocity[lNode_I*dof], vel, dof*sizeof(double) );
	}
}

void MeshAdvectionCorrection_EulerDeformCorrection( FeVariable *artDField, double *artVelocity, double dt ) {
	/* save the purely artificial bit of the remeshing in the artVelocity */
	FeMesh*      mesh = artDField->feMesh;
	Dof_Index    dof  = artDField->fieldComponentCount;
	double       artV[3], artD[3];
	int numLocalNodes = FeMesh_GetNodeLocalSize( mesh );
	Dof_Index    dof_I, lNode_I;

		/* INITIAL CONDITION: artV = 0 */
	if( dt == 0 ) {
		for( lNode_I = 0 ; lNode_I < numLocalNodes ; lNode_I++ ) {
			for( dof_I = 0 ; dof_I < dof ; dof_I++ ) 
				artV[dof_I] = 0;

			memcpy( &artVelocity[lNode_I*dof] , artV, dof*sizeof(double) );
		}
		return;
	}
	
	 /* GENERAL algorithm artV = -1*artD / dt. It's -ve because we 
		* want to reverse the effects of the artificial displacement */
	for( lNode_I = 0 ; lNode_I < numLocalNodes ; lNode_I++ ) {
		FeVariable_GetValueAtNode( artDField, lNode_I, artD );
		/* artV = artD / dt */
		for( dof_I = 0 ; dof_I < dof ; dof_I++ ) 
			artV[dof_I] = -1*artD[dof_I] / dt;

		memcpy( &artVelocity[lNode_I*dof] , artV, dof*sizeof(double) );
	}

}

void MeshAdvectionCorrection( void* sle, void* data ) {
	UnderworldContext*                                      context                 = (UnderworldContext*) data;
	Underworld_MeshAdvectionCorrection_ContextExt*          plugin;
	FeVariable*                                             velocityField           = context->velocityField;
	double dt = context->dt;
	double *artVelocity, *oldVelocity;
	int lNodeCount;
	lNodeCount = FeMesh_GetNodeLocalSize( velocityField->feMesh );

	artVelocity = oldVelocity = NULL;

	/* store the current velocity in oldVelocity */
	oldVelocity = Memory_Alloc_Array( double, 
			lNodeCount * velocityField->fieldComponentCount, 
			"artificial nodal velocities" );

	/* get the plugin from the context */
	plugin = ExtensionManager_Get( 
		context->extensionMgr, 
		context, 
		Underworld_MeshAdvectionCorrection_ContextExtHandle );

	MeshAdvectionCorrection_StoreCurrentVelocity( velocityField, oldVelocity );

	if( plugin->artDisplacement ) {
		artVelocity = Memory_Alloc_Array( double, lNodeCount * velocityField->fieldComponentCount, "artificial nodal velocities" );
		MeshAdvectionCorrection_EulerDeformCorrection( plugin->artDisplacement, artVelocity, dt );
	}

	/* Correct velocity and re-sync shadow space */
	MeshAdvectionCorrection_AddCorrection( velocityField, oldVelocity, artVelocity, MINUS );
	FeVariable_SyncShadowValues( velocityField );

	/* Solve Energy equation */
	plugin->energySolverExecute( sle, context );

	/* Reverse correction and re-sync */
	MeshAdvectionCorrection_AddCorrection( velocityField, oldVelocity, artVelocity, ADD );
	FeVariable_SyncShadowValues( velocityField );

	if( plugin->artDisplacement )
		Memory_Free( artVelocity );
	Memory_Free( oldVelocity );
}

void _Underworld_MeshAdvectionCorrection_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	UnderworldContext*                                      context = 
	Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data ); 
	Underworld_MeshAdvectionCorrection_ContextExt*       plugin;
	
	Journal_DFirewall( 
		(Bool)context, 
		Journal_Register( Error_Type, Underworld_MeshAdvectionCorrection_Type ), 
		"No context found\n" );
	Journal_DFirewall( 
		(Bool)context->energySLE, 
		Journal_Register( Error_Type, Underworld_MeshAdvectionCorrection_Type ), 
		"The required energy SLE component has not been created or placed on the context.\n");	
	
	/* Add the extension to the context */
	Underworld_MeshAdvectionCorrection_ContextExtHandle = ExtensionManager_Add(
		context->extensionMgr, 
		Underworld_MeshAdvectionCorrection_Type, 
		sizeof( Underworld_MeshAdvectionCorrection_ContextExt ) );
	plugin = ExtensionManager_Get( 
		context->extensionMgr, 
		context, 
		Underworld_MeshAdvectionCorrection_ContextExtHandle );

	if( Stg_ComponentFactory_GetRootDictBool( cf, "MeshAdvectionCorrection_UseArtDisplacementField", False) ) {
		/* get the artificial displacement field */
		plugin->artDisplacement = Stg_ComponentFactory_ConstructByName( cf, "ArtDisplacementField", FeVariable, True, data );
	}

	/* Replace the energy SLE's execute with this one. Save the old value for use later. */
	plugin->energySolverExecute = context->energySLE->_execute;
	context->energySLE->_execute = MeshAdvectionCorrection;
}

/* This function will provide StGermain the abilty to instantiate (create) this codelet on demand. */
void* _Underworld_MeshAdvectionCorrection_DefaultNew( Name name ) {
	return Codelet_New(
			Underworld_MeshAdvectionCorrection_Type,
			_Underworld_MeshAdvectionCorrection_DefaultNew,
			_Underworld_MeshAdvectionCorrection_Construct, /* SQ NOTE: Used to be a construct extensions. */
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
}

/* This function is automatically run by StGermain when this plugin is loaded. The name must be "<plugin-name>_Register". */
Index Underworld_MeshAdvectionCorrection_Register( PluginsManager* pluginsManager ) {
	/* A plugin is only properly registered once it returns the handle provided when submitting a codelet to StGermain. */
	return PluginsManager_Submit( 
		pluginsManager, 
		Underworld_MeshAdvectionCorrection_Type, 
		"0", 
		_Underworld_MeshAdvectionCorrection_DefaultNew );
}
