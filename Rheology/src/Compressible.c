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
** $Id: Compressible.c 610 2007-10-11 08:09:29Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "Compressible.h"
#include "RheologyMaterial.h"

#include <assert.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type Compressible_Type = "Compressible";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
Compressible* _Compressible_New( 
		SizeT                                              sizeOfSelf,
		Type                                               type,
		Stg_Class_DeleteFunction*                          _delete,
		Stg_Class_PrintFunction*                           _print,
		Stg_Class_CopyFunction*                            _copy, 
		Stg_Component_DefaultConstructorFunction*          _defaultConstructor,
		Stg_Component_ConstructFunction*                   _construct,
		Stg_Component_BuildFunction*                       _build,
		Stg_Component_InitialiseFunction*                  _initialise,
		Stg_Component_ExecuteFunction*                     _execute,
		Stg_Component_DestroyFunction*                     _destroy,
		StiffnessMatrixTerm_AssembleElementFunction*       _assembleElement,
		Name                                               name ) 
{
	Compressible*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(Compressible) );
	self = (Compressible*) _StiffnessMatrixTerm_New( 
			sizeOfSelf,
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
			_assembleElement,
			name );
	
	return self;
}

void _Compressible_Init(
		Compressible*        self, 
		FeMesh*		     geometryMesh, 
		Materials_Register*  materials_Register,
		double               oneOnLambda )
{
	self->isConstructed = True;
	
	self->oneOnLambda            = oneOnLambda; 
	self->geometryMesh           = geometryMesh;
	self->materials_Register     = materials_Register;
}


void _Compressible_Delete( void* compressible ) {
	Compressible*    self = (Compressible*)compressible;

	_StiffnessMatrixTerm_Delete( self );
}
void _Compressible_Print( void* compressible, Stream* stream ) {
	Compressible*    self = (Compressible*)compressible;

	_StiffnessMatrixTerm_Print( self, stream );

	Journal_PrintValue( stream, self->oneOnLambda );
}

void* _Compressible_DefaultNew( Name name ) {
	return (void*) _Compressible_New(
		sizeof(Compressible),
		Compressible_Type,
		_Compressible_Delete,
		_Compressible_Print,
		NULL,
		_Compressible_DefaultNew,
		_Compressible_Construct,
		_Compressible_Build,
		_Compressible_Initialise,
		_Compressible_Execute,
		_Compressible_Destroy,
		_Compressible_AssembleElement,
		name );
}


void _Compressible_Construct( void* compressible, Stg_ComponentFactory* cf, void* data ){
	Compressible*    self = (Compressible*)compressible;
	FeMesh*		 geometryMesh;
	Materials_Register* materials_Register;

	_StiffnessMatrixTerm_Construct( self, cf, data );

	geometryMesh = Stg_ComponentFactory_ConstructByKey( cf, self->name, "GeometryMesh", FeMesh, True, data );

	materials_Register = Stg_ObjectList_Get( cf->registerRegister, "Materials_Register" );
	assert( materials_Register );

	_Compressible_Init( 
			self, 
			geometryMesh,
			materials_Register,
			Stg_ComponentFactory_GetDouble( cf, self->name, "oneOnLambda", 10.0 ) );

	/* Make sure that we are using the correct type of swarm 
	 * SHOULDN'T THIS BEEN DONE AT THE STIFFNESSMATRIXTERM LEVEL? */
	Stg_CheckType( self->integrationSwarm, IntegrationPointsSwarm );
}

void _Compressible_Build( void* compressible, void* data ){
	_StiffnessMatrixTerm_Build( compressible, data );
}
void _Compressible_Initialise( void* compressible, void* data ){
	_StiffnessMatrixTerm_Initialise( compressible, data );
}
void _Compressible_Execute( void* compressible, void* data ){
	_StiffnessMatrixTerm_Execute( compressible, data );
}
void _Compressible_Destroy( void* compressible, void* data ){
	_StiffnessMatrixTerm_Destroy( compressible, data );
}

void _Compressible_AssembleElement(
		void*                                              compressible,
		StiffnessMatrix*                                   stiffnessMatrix,
	 	Element_LocalIndex                                 lElement_I,
		SystemLinearEquations*                             sle,
		FiniteElementContext*                              context,
		double**                                           elStiffMat )
{
	Compressible*             self                = (Compressible*)compressible;
	IntegrationPointsSwarm*   swarm               = (IntegrationPointsSwarm*) self->integrationSwarm;
	RheologyMaterial*         material;
	FeVariable*               variable1           = stiffnessMatrix->rowVariable;
	Dimension_Index           dim                 = stiffnessMatrix->dim;
	IntegrationPoint*         particle;
	Particle_InCellIndex      cParticle_I;
	Particle_InCellIndex      cellParticleCount;
	Element_NodeIndex         elementNodeCount;
	Index                     row_I;
	Index                     col_I;
	double                    detJac;
	Cell_Index                cell_I;
	ElementType*              elementType;
	Dof_Index                 dofCount;
	FeMesh*       		  mesh                = variable1->feMesh;
	double                    Ni[27];
	double*                   xi;
	double                    factor;
	FeMesh*  		  geometryMesh        = self->geometryMesh;
	ElementType*              geometryElementType;
	Particle_Index            lParticle_I;
	double oneOnLambda = 0.0;
	Bool oneToMany;

	/* Set the element type */
	elementType         = FeMesh_GetElementType( mesh, lElement_I );
	geometryElementType = FeMesh_GetElementType( geometryMesh, lElement_I );
	elementNodeCount    = elementType->nodeCount;
	dofCount            = elementNodeCount;

	/* Get number of particles per element */
	cell_I            = CellLayout_MapElementIdToCellId( swarm->cellLayout, lElement_I );
	cellParticleCount = swarm->cellParticleCountTbl[ cell_I ];

	/*
	 * Keep a flag indicating whether we are usinga one-to-one swarm mapper or not.
	 */

	oneToMany = Stg_Class_IsInstance(((IntegrationPointsSwarm*)self->integrationSwarm)->mapper, OneToManyMapper_Type);
	
	/* Loop over points to build Stiffness Matrix */
	for ( cParticle_I = 0 ; cParticle_I < cellParticleCount ; cParticle_I++ ) {
		particle = (IntegrationPoint*) Swarm_ParticleInCellAt( swarm, cell_I, cParticle_I );
		lParticle_I = swarm->cellParticleTbl[cell_I][cParticle_I];

		if(oneToMany) {
		    /*
		     * We're dealing with a one-to-many mapper. We will assemble each material point's
		     * constitutive matrix and combine them using their weights.
		     */

		    OneToManyRef *ref;
		    MaterialPointsSwarm *matSwarm;
		    int isComp = 0;
		    int ii;

		    matSwarm = ((OneToManyMapper*)((IntegrationPointsSwarm*)swarm)->mapper)->materialSwarm;
		    ref = OneToManyMapper_GetMaterialRef(((IntegrationPointsSwarm*)swarm)->mapper, particle);
		    for(ii = 0; ii < ref->numParticles; ii++) {
			material = (RheologyMaterial*)MaterialPointsSwarm_GetMaterialAt(matSwarm, ref->particleInds[ii]);
			if(!material->compressible)
			    continue;
			isComp++;
			oneOnLambda += material->compressible->oneOnLambda;
		    }

		    if(((float)isComp)/((float)ref->numParticles) < 0.5)
			continue;
		    oneOnLambda /= ((double)ref->numParticles);
		}
		else {

		    material = (RheologyMaterial*) IntegrationPointsSwarm_GetMaterialAt( swarm, lParticle_I );

		    /* Only make contribution to the compressibility matrix if this material is compressible */
		    if ( !material->compressible ) 
			continue;

		    oneOnLambda = material->compressible->oneOnLambda;
		}

		/* Calculate Determinant of Jacobian and Shape Functions */
		xi = particle->xi;
		detJac = ElementType_JacobianDeterminant( geometryElementType, geometryMesh, lElement_I, xi, dim );
		ElementType_EvaluateShapeFunctionsAt( elementType, xi, Ni );
		factor = detJac * particle->weight * oneOnLambda;

		for( row_I = 0 ; row_I < dofCount ; row_I++ ) {
			for( col_I = 0 ; col_I < dofCount ; col_I++ ) {
				elStiffMat[ row_I ][ col_I ] -= factor * Ni[ row_I ] * Ni[ col_I ];
			}
		}
	}
}
