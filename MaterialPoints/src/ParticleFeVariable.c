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
** $Id: ParticleFeVariable.c 595 2008-09-11 05:26:49Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PopulationControl/PopulationControl.h>
#include <PICellerator/Weights/Weights.h>

#include "MaterialPoints.h"

#include <assert.h>

const Type ParticleFeVariable_Type = "ParticleFeVariable";

char* ParticleFeVariable_names[10];
int ParticleFeVariable_nNames = 0;
int ParticleFeVariable_curName = 0;

ParticleFeVariable* _ParticleFeVariable_New(  PARTICLEFEVARIABLE_DEFARGS  ) {
	ParticleFeVariable* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(ParticleFeVariable) );
	/* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
	/* This means that any values of these parameters that are passed into this function are not passed onto the parent function
	   and so should be set to ZERO in any children of this class. */
	nameAllocationType = NON_GLOBAL;
	_syncShadowValues  = _FeVariable_SyncShadowValues;

	self = (ParticleFeVariable*) _FeVariable_New(  FEVARIABLE_PASSARGS  );

	self->_valueAtParticle = _valueAtParticle;
	
	return self;
}

void _ParticleFeVariable_Init( ParticleFeVariable* self, IntegrationPointsSwarm* swarm ) {
	/* Create Vector */
	self->assemblyVectorName = Stg_Object_AppendSuffix( self, "assemblyVector" );
	self->assemblyVector = ForceVector_New( 
		self->assemblyVectorName,
		(FiniteElementContext*) self->context,
		(FeVariable*) self, 
		self->dim, 
		self->context->entryPoint_Register, 
		self->communicator );
	self->assemblyTerm = ForceTerm_New( "assemblyTerm", (FiniteElementContext*) self->context, self->assemblyVector, (Swarm*)swarm, (Stg_Component*) self );

	self->massMatrixName = Stg_Object_AppendSuffix( self, "massMatrix" );
	self->massMatrix = ForceVector_New( 
		self->massMatrixName,
		(FiniteElementContext*) self->context,
		(FeVariable*) self, 
		self->dim, 
		self->context->entryPoint_Register, 
		self->communicator );
	self->massMatrixForceTerm = ForceTerm_New( "massMatrixForceTerm", (FiniteElementContext*) self->context, self->massMatrix, (Swarm*)swarm, (Stg_Component*) self );
	ForceTerm_SetAssembleElementFunction( self->massMatrixForceTerm, ParticleFeVariable_AssembleElementShapeFunc );
   
   /* Changed by PatrickSunter, 8/7/2009 - used to append onto stokesEqn-execute. However
    * this component should work for other FE systems other than the Stokes solver */
	EP_InsertClassHookBefore( Context_GetEntryPoint( self->context, AbstractContext_EP_UpdateClass ), "TimeIntegrator_UpdateClass", _ParticleFeVariable_Execute, self );
	ParticleFeVariable_names[ParticleFeVariable_nNames++] = self->name;

	self->useDeriv = False;
	self->GNx = NULL;
}

/* --- Virtual Function Implementations --- */
void _ParticleFeVariable_Delete( void* materialFeVariable ) {
	ParticleFeVariable* self = (ParticleFeVariable*) materialFeVariable;

	_FeVariable_Delete( self );
}

void _ParticleFeVariable_Print( void* materialFeVariable, Stream* stream ) {
	ParticleFeVariable* self = (ParticleFeVariable*) materialFeVariable;
	
	/* General info */
	Journal_Printf( stream, "ParticleFeVariable (ptr): %p\n", self );
	
	/* Print parent */
	_FeVariable_Print( self, stream );
	
	/* ParticleFeVariable info */
}

void* _ParticleFeVariable_Copy( void* feVariable, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	abort();
	
	return NULL;
}

void _ParticleFeVariable_AssignFromXML( void* materialFeVariable, Stg_ComponentFactory* cf, void* data ){
	ParticleFeVariable*		self = (ParticleFeVariable*) materialFeVariable;
	IntegrationPointsSwarm*	swarm;
	FiniteElementContext*	context;
	FeMesh*						mesh;

	_FieldVariable_AssignFromXML( self, cf, data );

	swarm = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Swarm", IntegrationPointsSwarm, True, data );
	context = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Context", FiniteElementContext, False, data );
	mesh = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Mesh", FeMesh, True, data );

	_FeVariable_Init( (FeVariable*)self, mesh, NULL, NULL, NULL, NULL, NULL, NULL, False, False );
	_ParticleFeVariable_Init( self, swarm );
}

void _ParticleFeVariable_Build( void* materialFeVariable, void* data ) {
	ParticleFeVariable* self = (ParticleFeVariable*) materialFeVariable;
	int dataSize;

	ForceTerm_SetAssembleElementFunction( self->assemblyTerm, ParticleFeVariable_AssembleElement );
#if 0
	if( self->useDeriv )
	   ForceTerm_SetAssembleElementFunction( self->assemblyTerm, ParticleFeVariable_AssembleElement_Deriv );
#endif
	
	Stg_Component_Build( self->feMesh, data, False );
	dataSize = FeMesh_GetNodeDomainSize( self->feMesh ) * self->fieldComponentCount;
	self->data = Memory_Alloc_Array( double, dataSize, "data" );

	/* Do a Variable_Update() first as well as last, since if we are loading from checkpoint we need
	to make sure the variable exists to put ICs onto - and we only just allocated it */
	Stg_Component_Build( self->dataVariable, data, False );
	Variable_Update( self->dataVariable );

	_FeVariable_Build( self, data );

	Stg_Component_Build( self->assemblyVector, data, False );
	Stg_Component_Build( self->massMatrix, data, False );
	Stg_Component_Build( self->assemblyTerm, data, False );

	Variable_Update( self->dataVariable );
}

void _ParticleFeVariable_Initialise( void* materialFeVariable, void* data ) {
	ParticleFeVariable*      self = (ParticleFeVariable*) materialFeVariable;
	DomainContext*   context = (DomainContext*)data;

	/* Do a Variable_Update() first as well as last, since if we are loading from checkpoint we need
	to make sure the variable exists to put ICs onto */

	Stg_Component_Initialise( self->assemblyVector, data, False );
	Stg_Component_Initialise( self->massMatrix, data, False );
	Stg_Component_Initialise( self->assemblyTerm, data, False );
	Stg_Component_Initialise( self->dataVariable, data, False );
	Variable_Update( self->dataVariable );

	_FeVariable_Initialise( self, data );

	Variable_Update( self->dataVariable );
	/* If loading from CP, _don't_ recalculate the field as we've already just loaded it!
		-- PatrickSunter 22 Nov 2006 */
	if ( !(context && (True == context->loadFromCheckPoint)  ) ) {
		ParticleFeVariable_Update( self );
	}
}

void _ParticleFeVariable_Execute( void* materialFeVariable, void* _ctx ) {
   AbstractContext* ctx = (AbstractContext*)_ctx;

   if( ParticleFeVariable_curName >= ParticleFeVariable_nNames )
      ParticleFeVariable_curName = 0;
   ParticleFeVariable_Update( LiveComponentRegister_Get( ctx->CF->LCRegister, ParticleFeVariable_names[ParticleFeVariable_curName++] ) );
}

void _ParticleFeVariable_Destroy( void* materialFeVariable, void* data ) {
	ParticleFeVariable* self = (ParticleFeVariable*) materialFeVariable;

	Memory_Free( self->assemblyVectorName );
	Memory_Free( self->massMatrixName );
	Memory_Free( self->data );

	Stg_Component_Destroy( self->assemblyVector, data, False );
	Stg_Component_Destroy( self->assemblyTerm, data, False );

	Stg_Component_Destroy( self->massMatrix, data, False );
	Stg_Component_Destroy( self->massMatrixForceTerm, data, False );

	_FeVariable_Destroy( self, data );
}

void ParticleFeVariable_Update( void* materialFeVariable ) {
	ParticleFeVariable* self = (ParticleFeVariable*) materialFeVariable;

	/* Initialise Vectors */
	VecSet( self->assemblyVector->vector, 0.0 );
	VecSet( self->massMatrix->vector, 0.0 );

	ForceVector_Assemble( self->assemblyVector );
	ForceVector_Assemble( self->massMatrix );

	VecPointwiseDivide( self->assemblyVector->vector, self->assemblyVector->vector, self->massMatrix->vector );

	SolutionVector_UpdateSolutionOntoNodes( self->assemblyVector );
}

void ParticleFeVariable_AssembleElement( void* _forceTerm, ForceVector* forceVector, Element_LocalIndex lElement_I, double* elForceVector ) {
	ForceTerm*                 forceTerm         = (ForceTerm*) _forceTerm;
	ParticleFeVariable*        self              = Stg_CheckType( forceVector->feVariable, ParticleFeVariable );
	IntegrationPointsSwarm*    swarm             = (IntegrationPointsSwarm*)forceTerm->integrationSwarm;
	FeMesh*							mesh              = self->feMesh;
	Element_NodeIndex          elementNodeCount  = FeMesh_GetElementNodeSize( mesh, lElement_I );
	ElementType*               elementType       = FeMesh_GetElementType( mesh, lElement_I );
	Cell_Index                 cell_I            = CellLayout_MapElementIdToCellId( swarm->cellLayout, lElement_I );
	Particle_InCellIndex       cellParticleCount;
	Particle_InCellIndex       cParticle_I;
	IntegrationPoint*          particle;
	Node_Index                 node_I;
	Dof_Index                  dofCount          = self->fieldComponentCount;
	Dof_Index                  dof_I;
	Dof_Index                  dim               = self->dim;
	double                     shapeFunc[8], detJac;
	double                     particleValue[9];

	cellParticleCount = swarm->cellParticleCountTbl[ cell_I ];
	
	for( cParticle_I = 0 ; cParticle_I < cellParticleCount; cParticle_I++ ) {
		/* Find this particle in the element */
		particle = (IntegrationPoint*) Swarm_ParticleInCellAt( swarm, cell_I, cParticle_I );

		self->currentParticleIndex = swarm->cellParticleTbl[cell_I][cParticle_I];
		ParticleFeVariable_ValueAtParticle( self, swarm, lElement_I, particle, particleValue );

		/* get shape function and detJac */
		ElementType_EvaluateShapeFunctionsAt( elementType, particle->xi, shapeFunc );
		detJac = ElementType_JacobianDeterminant( elementType, mesh, lElement_I, particle->xi, dim );

		for ( dof_I = 0 ; dof_I < dofCount ; dof_I++ ) {
			for ( node_I = 0 ; node_I < elementNodeCount ; node_I++ ) {
				elForceVector[ node_I * dofCount + dof_I ] += particle->weight * detJac * shapeFunc[ node_I ] * particleValue[ dof_I ]; 
			}
		}
	}
}

void ParticleFeVariable_AssembleElementShapeFunc( void* _forceTerm, ForceVector* forceVector, Element_LocalIndex lElement_I, double* elForceVector ) {
	ForceTerm*                 forceTerm         = (ForceTerm*) _forceTerm;
	ParticleFeVariable*        self              = Stg_CheckType( forceVector->feVariable, ParticleFeVariable );
	Swarm*                     swarm             = forceTerm->integrationSwarm;
	FeMesh*							mesh              = self->feMesh;
	Element_NodeIndex          elementNodeCount  = FeMesh_GetElementNodeSize( mesh, lElement_I );
	ElementType*               elementType       = FeMesh_GetElementType( mesh, lElement_I );
	Cell_Index                 cell_I            = CellLayout_MapElementIdToCellId( swarm->cellLayout, lElement_I );
	int                        dim = Mesh_GetDimSize( mesh );
	Particle_InCellIndex       cellParticleCount;
	Particle_InCellIndex       cParticle_I;
	IntegrationPoint*          particle;
	double                     detJac;
	double                     shapeFunc[8];
	Node_Index                 node_I;
	Dof_Index                  dofCount          = self->fieldComponentCount;
	Dof_Index                  dof_I;

	cellParticleCount = swarm->cellParticleCountTbl[ cell_I ];
	
	for( cParticle_I = 0 ; cParticle_I < cellParticleCount; cParticle_I++ ) {
		/* Find this particle in the element */
		particle = (IntegrationPoint*) Swarm_ParticleInCellAt( swarm, cell_I, cParticle_I );

		ElementType_EvaluateShapeFunctionsAt( elementType, particle->xi, shapeFunc );
		detJac = ElementType_JacobianDeterminant( elementType, mesh, lElement_I, particle->xi, dim );

		for ( dof_I = 0 ; dof_I < dofCount ; dof_I++ ) {
			for ( node_I = 0 ; node_I < elementNodeCount ; node_I++ ) {
			   elForceVector[ node_I * dofCount + dof_I ] += particle->weight * detJac * shapeFunc[node_I];
			}
		}
	}
}


