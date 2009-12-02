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
** $Id: ConstitutiveMatrixCartesian.c 803 2008-09-11 05:22:20Z LukeHodkinson $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/



#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "Rheology_Register.h"
#include "ConstitutiveMatrix.h"
#include "ConstitutiveMatrixCartesian.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>


/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type ConstitutiveMatrixCartesian_Type = "ConstitutiveMatrixCartesian";

ConstitutiveMatrixCartesian* ConstitutiveMatrixCartesian_New( 
		Name                                                name,
		StiffnessMatrix*                                    stiffnessMatrix,
		Swarm*                                              swarm,
		Dimension_Index                                     dim,
		FiniteElementContext*                               context,
		Materials_Register*                                 materials_Register )
{
	ConstitutiveMatrixCartesian* self = (ConstitutiveMatrixCartesian*) _ConstitutiveMatrixCartesian_DefaultNew( name );

	ConstitutiveMatrixCartesian_InitAll( self, stiffnessMatrix, swarm, dim, context, materials_Register );

	return self;
}

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
ConstitutiveMatrixCartesian* _ConstitutiveMatrixCartesian_New( 
		SizeT                                               sizeOfSelf,  
		Type                                                type,
		Stg_Class_DeleteFunction*                           _delete,
		Stg_Class_PrintFunction*                            _print,
		Stg_Class_CopyFunction*                             _copy, 
		Stg_Component_DefaultConstructorFunction*           _defaultConstructor,
		Stg_Component_ConstructFunction*                    _construct,
		Stg_Component_BuildFunction*                        _build,
		Stg_Component_InitialiseFunction*                   _initialise,
		Stg_Component_ExecuteFunction*                      _execute,
		Stg_Component_DestroyFunction*                      _destroy,
		StiffnessMatrixTerm_AssembleElementFunction*        _assembleElement,		
		ConstitutiveMatrix_SetValueFunc*                    _setValue,
		ConstitutiveMatrix_GetValueFunc*                    _getViscosity,
		ConstitutiveMatrix_SetValueFunc*                    _isotropicCorrection,
		ConstitutiveMatrix_SetSecondViscosityFunc*          _setSecondViscosity,
		ConstitutiveMatrix_Assemble_D_B_Func*               _assemble_D_B,
		ConstitutiveMatrix_CalculateStressFunc*             _calculateStress,
		Name                                                name )
{
	ConstitutiveMatrixCartesian* self;
	
	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(ConstitutiveMatrixCartesian) );
	self = (ConstitutiveMatrixCartesian*) _ConstitutiveMatrix_New( 
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
		_setValue,
		_getViscosity,
		_isotropicCorrection,
		_setSecondViscosity,
		_assemble_D_B,
		_calculateStress,
		name );
	
	/* Function pointers for this class that are not on the parent class should be set here */
	
	return self;
}

void _ConstitutiveMatrixCartesian_Init( 
		ConstitutiveMatrixCartesian*                 self )
{
	self->rowSize = self->columnSize = StGermain_nSymmetricTensorVectorComponents( self->dim );
	self->Dtilda_B = Memory_Alloc_2DArray( double, self->rowSize, self->dim, "D~ times B matrix" );

	if( self->dim == 2 ) {
		self->_setValue = _ConstitutiveMatrixCartesian2D_SetValueInAllEntries;
		self->_setSecondViscosity = _ConstitutiveMatrixCartesian2D_SetSecondViscosity;
		self->_getViscosity = _ConstitutiveMatrixCartesian2D_GetIsotropicViscosity;
		self->_isotropicCorrection = _ConstitutiveMatrixCartesian2D_IsotropicCorrection;
		self->_assemble_D_B = _ConstitutiveMatrixCartesian2D_Assemble_D_B;
		self->_calculateStress = _ConstitutiveMatrixCartesian2D_CalculateStress;
	} else {
		self->_setValue = _ConstitutiveMatrixCartesian3D_SetValueInAllEntries;
		self->_setSecondViscosity = _ConstitutiveMatrixCartesian3D_SetSecondViscosity;
		self->_getViscosity = _ConstitutiveMatrixCartesian3D_GetIsotropicViscosity;
		self->_isotropicCorrection = _ConstitutiveMatrixCartesian3D_IsotropicCorrection;
		self->_assemble_D_B = _ConstitutiveMatrixCartesian3D_Assemble_D_B;
		self->_calculateStress = _ConstitutiveMatrixCartesian3D_CalculateStress;
	}

  /* store each particle's constitutiveMatrix */
  if( self->storeConstitutiveMatrix ) 
		ConstitutiveMatrixCartesian_SetupParticleStorage( self );

}

void ConstitutiveMatrixCartesian_InitAll( 
		void*                                        constitutiveMatrix,
		StiffnessMatrix*                             stiffnessMatrix,
		Swarm*                                       swarm,
		Dimension_Index                              dim,
		FiniteElementContext*                        context,
		Materials_Register*                          materials_Register )
{
	ConstitutiveMatrixCartesian* self = (ConstitutiveMatrixCartesian*) constitutiveMatrix;

	ConstitutiveMatrix_InitAll( self, stiffnessMatrix, swarm, dim, context, materials_Register );
	_ConstitutiveMatrixCartesian_Init( self );
}

void _ConstitutiveMatrixCartesian_Delete( void* constitutiveMatrix ) {
	ConstitutiveMatrixCartesian* self = (ConstitutiveMatrixCartesian*)constitutiveMatrix;

	_ConstitutiveMatrix_Delete( self );
}

void _ConstitutiveMatrixCartesian_Print( void* constitutiveMatrix, Stream* stream ) {
	ConstitutiveMatrixCartesian* self = (ConstitutiveMatrixCartesian*)constitutiveMatrix;
	
	_ConstitutiveMatrix_Print( self, stream );

	/* General info */
}

void* _ConstitutiveMatrixCartesian_DefaultNew( Name name ) {
	return (void*)_ConstitutiveMatrixCartesian_New( 
		sizeof(ConstitutiveMatrixCartesian), 
		ConstitutiveMatrixCartesian_Type,
		_ConstitutiveMatrixCartesian_Delete,
		_ConstitutiveMatrixCartesian_Print,
		NULL,
		_ConstitutiveMatrixCartesian_DefaultNew,
		_ConstitutiveMatrixCartesian_Construct,
		_ConstitutiveMatrixCartesian_Build,
		_ConstitutiveMatrixCartesian_Initialise,
		_ConstitutiveMatrixCartesian_Execute,
		_ConstitutiveMatrixCartesian_Destroy,
		_ConstitutiveMatrixCartesian_AssembleElement,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		name );
}

void _ConstitutiveMatrixCartesian_Construct( void* constitutiveMatrix, Stg_ComponentFactory* cf, void* data ) {
	ConstitutiveMatrixCartesian*            self             = (ConstitutiveMatrixCartesian*)constitutiveMatrix;

	/* Construct Parent */
	_ConstitutiveMatrix_Construct( self, cf, data );

	_ConstitutiveMatrixCartesian_Init( self );
}

void _ConstitutiveMatrixCartesian_Build( void* constitutiveMatrix, void* data ) {
	ConstitutiveMatrixCartesian*             self             = (ConstitutiveMatrixCartesian*)constitutiveMatrix;

	_ConstitutiveMatrix_Build( self, data );
}

void _ConstitutiveMatrixCartesian_Initialise( void* constitutiveMatrix, void* data ) {
	ConstitutiveMatrixCartesian*             self             = (ConstitutiveMatrixCartesian*)constitutiveMatrix;

	_ConstitutiveMatrix_Initialise( self, data );
}

void _ConstitutiveMatrixCartesian_Execute( void* constitutiveMatrix, void* data ) {
	_ConstitutiveMatrix_Execute( constitutiveMatrix, data );
}

void _ConstitutiveMatrixCartesian_Destroy( void* constitutiveMatrix, void* data ) {
	ConstitutiveMatrixCartesian* self = (ConstitutiveMatrixCartesian*)constitutiveMatrix;

	_ConstitutiveMatrix_Destroy( constitutiveMatrix, data );

	Memory_Free( self->Dtilda_B );
	Memory_Free( self->Ni );
}

void _ConstitutiveMatrixCartesian_AssembleElement( 
		void*                                              constitutiveMatrix,
		StiffnessMatrix*                                   stiffnessMatrix, 
		Element_LocalIndex                                 lElement_I, 
		SystemLinearEquations*                             sle,
		FiniteElementContext*                              context,
		double**                                           elStiffMat ) 
{
	ConstitutiveMatrixCartesian*     self       = (ConstitutiveMatrixCartesian*) constitutiveMatrix;
	Swarm*                  swarm               = self->integrationSwarm;
	FeVariable*             variable1           = stiffnessMatrix->rowVariable;
	Dimension_Index         dim                 = stiffnessMatrix->dim;
	IntegrationPoint*       particle;
	Particle_InCellIndex	cParticle_I;
	Particle_InCellIndex	cellParticleCount;
	Element_NodeIndex       elementNodeCount;
	Node_ElementLocalIndex  rowNode_I;
	Node_ElementLocalIndex  colNode_I;
	double**                GNx;
	double                  detJac;
	Cell_Index              cell_I;
	ElementType*            elementType;
        double                  Bj_x, Bj_y;
	double                  Bi_x; 
	double                  Bi_y;
	double                  Bi_z;
	Dof_Index               rowNodeDof_I;
	Dof_Index               colNodeDof_I;
	Dof_Index               nodeDofCount;
	double**                Dtilda_B;
	double                  vel[3], velDerivs[9], *Ni, eta;
	Bool oneToMany;

	self->sle = sle;

	/* Set the element type */
	elementType       = FeMesh_GetElementType( variable1->feMesh, lElement_I );
	elementNodeCount  = elementType->nodeCount;
	nodeDofCount      = dim;

	/* allocate */
	if( elementNodeCount > self->max_nElNodes ) {
		 self->max_nElNodes = elementNodeCount;
		 self->GNx = ReallocArray2D( self->GNx, double, dim, elementNodeCount );
		 self->Ni = ReallocArray( self->Ni, double, elementNodeCount );
	}
	GNx = self->GNx;
	Ni = self->Ni;
	Dtilda_B = self->Dtilda_B;

	/* Get number of particles per element */
	cell_I            = CellLayout_MapElementIdToCellId( swarm->cellLayout, lElement_I );
	cellParticleCount = swarm->cellParticleCountTbl[ cell_I ];

	/* Determine whether this is the first solve for not */
	Journal_Firewall( sle != NULL, Journal_Register( Error_Type, ConstitutiveMatrix_Type ), 
			"In func %s: SLE is NULL.\n", __func__ );

	/* Note: we may have deliberately set the previousSolutionExists flag to true in the
		parent ConstitutiveMatrix constructor if in restart mode, even if the SLE hasn't executed yet
		in this run - so only update to the sle's value when SLE is confirming it has
		executed */
	if ( True == sle->hasExecuted ) {
		self->previousSolutionExists = sle->hasExecuted;
	}
	self->sleNonLinearIteration_I = sle->nonLinearIteration_I;


	/*
	 * Keep a flag indicating whether we are usinga one-to-one swarm mapper or not.
	 */

	oneToMany = Stg_Class_IsInstance(((IntegrationPointsSwarm*)self->integrationSwarm)->mapper, OneToManyMapper_Type);

	
	/* Loop over points to build Stiffness Matrix */
	for ( cParticle_I = 0 ; cParticle_I < cellParticleCount ; cParticle_I++ ) {
		particle = (void*) Swarm_ParticleInCellAt( swarm, cell_I, cParticle_I );

		/* Calculate Determinant of Jacobian and Shape Function Global Derivatives */
		ElementType_ShapeFunctionsGlobalDerivs( 
			elementType,
			variable1->feMesh, lElement_I,
			particle->xi, dim, &detJac, GNx );

                /* Evalulate velocity and velocity derivatives at this particle. */
                FeVariable_InterpolateWithinElement(
                   variable1, lElement_I, particle->xi, vel );
                FeVariable_InterpolateDerivatives_WithGNx(
                   variable1, lElement_I, GNx, velDerivs );

		/*
		 * Assemble constitutive matrix. Note that we have to handle one-to-many swarms
		 * differently.
		 */

		if(!oneToMany) {
		    // TODO : pass in the context here?
		    ConstitutiveMatrix_Assemble( constitutiveMatrix, lElement_I,
						 swarm->cellParticleTbl[cell_I][cParticle_I], particle );
		}
		else {
		    /*
		     * We're dealing with a one-to-many mapper. We will assemble each material point's
		     * constitutive matrix and combine them using their weights.
		     */

		    OneToManyRef *ref;
		    double **matrixData;
		    int ii, jj, kk;

		    matrixData = Memory_Alloc_2DArray(double, self->columnSize, self->rowSize, self->name);
		    memset(matrixData[0], 0, self->columnSize*self->rowSize*sizeof(double));
		    ref = OneToManyMapper_GetMaterialRef(((IntegrationPointsSwarm*)swarm)->mapper, particle);
		    for(ii = 0; ii < ref->numParticles; ii++) {
			/* Assemble this material point. */
			ConstitutiveMatrix_AssembleMaterialPoint(
			    constitutiveMatrix, lElement_I,
			    ((OneToManyMapper*)((IntegrationPointsSwarm*)swarm)->mapper)->materialSwarm,
			    ref->particleInds[ii]);
			/* Add to cumulative matrix. */
			for(jj = 0; jj < self->rowSize; jj++) {
			    for(kk = 0; kk < self->columnSize; kk++)
				matrixData[jj][kk] += ref->weights[ii]*self->matrixData[jj][kk];
			}
		    }
		    /* Copy matrix data and free temporary array. */
		    memcpy(self->matrixData[0], matrixData[0], self->columnSize*self->rowSize*sizeof(double));
		    Memory_Free(matrixData);
		}

		eta = self->matrixData[2][2];

		/* Turn D Matrix into D~ Matrix by multiplying in the weight and the detJac (this is a shortcut for speed) */
		ConstitutiveMatrix_MultiplyByValue( constitutiveMatrix, detJac * particle->weight );

		for( rowNode_I = 0 ; rowNode_I < elementNodeCount ; rowNode_I++ ) {
			rowNodeDof_I = rowNode_I*nodeDofCount;
                        Bj_x = GNx[0][rowNode_I];
                        Bj_y = GNx[1][rowNode_I];
			
			/* Build D~ * B */
			ConstitutiveMatrix_Assemble_D_B( constitutiveMatrix, GNx, rowNode_I, Dtilda_B );
				
			for( colNode_I = 0 ; colNode_I < elementNodeCount ; colNode_I++ ) {
				colNodeDof_I = colNode_I*nodeDofCount;
				Bi_x = GNx[ I_AXIS ][colNode_I];
				Bi_y = GNx[ J_AXIS ][colNode_I];
				
				/* Build BTrans * ( D~ * B ) */
				if ( dim == 2 ) {
          if( !sle->nlFormJacobian ) {

             elStiffMat[ colNodeDof_I     ][ rowNodeDof_I     ] += Bi_x * Dtilda_B[0][0] + Bi_y * Dtilda_B[2][0];
             elStiffMat[ colNodeDof_I     ][ rowNodeDof_I + 1 ] += Bi_x * Dtilda_B[0][1] + Bi_y * Dtilda_B[2][1];
             elStiffMat[ colNodeDof_I + 1 ][ rowNodeDof_I     ] += Bi_y * Dtilda_B[1][0] + Bi_x * Dtilda_B[2][0];
             elStiffMat[ colNodeDof_I + 1 ][ rowNodeDof_I + 1 ] += Bi_y * Dtilda_B[1][1] + Bi_x * Dtilda_B[2][1];

          }
          else {
             double DuDx, DuDy, DvDx, DvDy;
             double DetaDu, DetaDv;
             double intFac, fac;

             DuDx = velDerivs[0]; DuDy = velDerivs[1];
             DvDx = velDerivs[2]; DvDy = velDerivs[3];
             DetaDu = self->derivs[0] * Bj_x + self->derivs[1] * Bj_y + self->derivs[2] * Ni[rowNode_I];
             DetaDv = self->derivs[3] * Bj_x + self->derivs[4] * Bj_y + self->derivs[5] * Ni[rowNode_I];
             intFac = particle->weight * detJac;

             fac = eta * Bj_y + DuDy * DetaDu + DvDx * DetaDu;
             elStiffMat[colNodeDof_I][rowNodeDof_I] +=
                intFac * (2.0 * Bi_x * (eta * Bj_x + DuDx * DetaDu) + Bi_y * fac);
             elStiffMat[colNodeDof_I + 1][rowNodeDof_I] +=
                intFac * (2.0 * Bi_y * DvDy * DetaDu + Bi_x * fac);

             fac = eta * Bj_x + DvDx * DetaDv + DuDy * DetaDv;
             elStiffMat[colNodeDof_I][rowNodeDof_I + 1] +=
                intFac * (2.0 * Bi_x * DuDx * DetaDv + Bi_y * fac);
             elStiffMat[colNodeDof_I + 1][rowNodeDof_I + 1] +=
                intFac * (2.0 * Bi_y * (eta * Bj_y + DvDy * DetaDv) + Bi_x * fac);

          }
				}
				else {
					Bi_z = GNx[ K_AXIS ][colNode_I];

					elStiffMat[ colNodeDof_I     ][ rowNodeDof_I     ] += 
						Bi_x * Dtilda_B[0][0] + Bi_y * Dtilda_B[3][0] + Bi_z * Dtilda_B[4][0];
					elStiffMat[ colNodeDof_I     ][ rowNodeDof_I + 1 ] += 
						Bi_x * Dtilda_B[0][1] + Bi_y * Dtilda_B[3][1] + Bi_z * Dtilda_B[4][1];
					elStiffMat[ colNodeDof_I     ][ rowNodeDof_I + 2 ] += 
						Bi_x * Dtilda_B[0][2] + Bi_y * Dtilda_B[3][2] + Bi_z * Dtilda_B[4][2];
						
					elStiffMat[ colNodeDof_I + 1 ][ rowNodeDof_I     ] += 
						Bi_y * Dtilda_B[1][0] + Bi_x * Dtilda_B[3][0] + Bi_z * Dtilda_B[5][0];
					elStiffMat[ colNodeDof_I + 1 ][ rowNodeDof_I + 1 ] += 
						Bi_y * Dtilda_B[1][1] + Bi_x * Dtilda_B[3][1] + Bi_z * Dtilda_B[5][1];
					elStiffMat[ colNodeDof_I + 1 ][ rowNodeDof_I + 2 ] += 
						Bi_y * Dtilda_B[1][2] + Bi_x * Dtilda_B[3][2] + Bi_z * Dtilda_B[5][2];
						
					elStiffMat[ colNodeDof_I + 2 ][ rowNodeDof_I     ] += 
						Bi_z * Dtilda_B[2][0] + Bi_x * Dtilda_B[4][0] + Bi_y * Dtilda_B[5][0];
					elStiffMat[ colNodeDof_I + 2 ][ rowNodeDof_I + 1 ] += 
						Bi_z * Dtilda_B[2][1] + Bi_x * Dtilda_B[4][1] + Bi_y * Dtilda_B[5][1];
					elStiffMat[ colNodeDof_I + 2 ][ rowNodeDof_I + 2 ] += 
						Bi_z * Dtilda_B[2][2] + Bi_x * Dtilda_B[4][2] + Bi_y * Dtilda_B[5][2];
				}
			}
		}
	}
}

void _ConstitutiveMatrixCartesian2D_SetValueInAllEntries( void* constitutiveMatrix, double value ) {
	ConstitutiveMatrix* self   = (ConstitutiveMatrix*) constitutiveMatrix;

	if ( fabs( value ) < 1.0e-20 ) 
		ConstitutiveMatrix_ZeroMatrix( self );
	else {
		double**            D      = self->matrixData;

		D[0][0] = D[0][1] = D[0][2] = value;
		D[1][0] = D[1][1] = D[1][2] = value;
		D[2][0] = D[2][1] = D[2][2] = value;
	
		self->isDiagonal = False;
	}
}

void _ConstitutiveMatrixCartesian3D_SetValueInAllEntries( void* _constitutiveMatrix, double value ) {
	ConstitutiveMatrix* self   = (ConstitutiveMatrix*)_constitutiveMatrix;

	if ( fabs( value ) < 1.0e-20 ) 
		ConstitutiveMatrix_ZeroMatrix( self );
	else {
		double**            D      = self->matrixData;

		D[0][0] = D[0][1] = D[0][2] = D[0][3] = D[0][4] = D[0][5] = value;
		D[1][0] = D[1][1] = D[1][2] = D[1][3] = D[1][4] = D[1][5] = value;
		D[2][0] = D[2][1] = D[2][2] = D[2][3] = D[2][4] = D[2][5] = value;
		D[3][0] = D[3][1] = D[3][2] = D[3][3] = D[3][4] = D[3][5] = value;
		D[4][0] = D[4][1] = D[4][2] = D[4][3] = D[4][4] = D[4][5] = value;
		D[5][0] = D[5][1] = D[5][2] = D[5][3] = D[5][4] = D[5][5] = value;
	
		self->isDiagonal = False;
	}
}

double _ConstitutiveMatrixCartesian2D_GetIsotropicViscosity( void* constitutiveMatrix ) {
	ConstitutiveMatrix* self = (ConstitutiveMatrix*) constitutiveMatrix;

	return self->matrixData[2][2];
}

double _ConstitutiveMatrixCartesian3D_GetIsotropicViscosity( void* constitutiveMatrix ) {
	ConstitutiveMatrix* self = (ConstitutiveMatrix*) constitutiveMatrix;

	return self->matrixData[3][3];
}

void _ConstitutiveMatrixCartesian2D_IsotropicCorrection( void* constitutiveMatrix, double isotropicCorrection ) {
	ConstitutiveMatrix* self   = (ConstitutiveMatrix*) constitutiveMatrix;
	double**            D      = self->matrixData;
		
	D[0][0] += 2.0 * isotropicCorrection;
	D[1][1] += 2.0 * isotropicCorrection;
	D[2][2] += isotropicCorrection;
}

void _ConstitutiveMatrixCartesian3D_IsotropicCorrection( void* constitutiveMatrix, double isotropicCorrection ) {
	ConstitutiveMatrix* self   = (ConstitutiveMatrix*) constitutiveMatrix;
	double**            D      = self->matrixData;

	D[0][0] += 2.0 * isotropicCorrection;
	D[1][1] += 2.0 * isotropicCorrection;
	D[2][2] += 2.0 * isotropicCorrection;
	
	D[3][3] += isotropicCorrection;
	D[4][4] += isotropicCorrection;
	D[5][5] += isotropicCorrection;
}

void _ConstitutiveMatrixCartesian2D_SetSecondViscosity( void* constitutiveMatrix, double deltaViscosity, XYZ director ) {
	ConstitutiveMatrix* self      = (ConstitutiveMatrix*) constitutiveMatrix;
	double**            D         = self->matrixData;
	double              n1        = director[ I_AXIS ];
	double              n2        = director[ J_AXIS ];
	double              a0;
	double              a1;

	a0 = 4.0 * deltaViscosity * n1 * n1 * n2 * n2;
	a1 = 2.0 * deltaViscosity * n1 * n2 * (n2*n2 - n1*n1);

	D[0][0] += -a0 ;	D[0][1] +=  a0 ;	D[0][2] += -a1 ;
	D[1][0] +=  a0 ;	D[1][1] += -a0 ;	D[1][2] +=  a1 ;
	D[2][0] += -a1 ;	D[2][1] +=  a1 ;	D[2][2] +=  a0 - deltaViscosity ;

	self->isDiagonal = False;
}

void _ConstitutiveMatrixCartesian3D_SetSecondViscosity( void* constitutiveMatrix, double deltaViscosity, XYZ director ) {
	ConstitutiveMatrix* self      = (ConstitutiveMatrix*) constitutiveMatrix;
	double**            D         = self->matrixData;
	double              n1        = director[ I_AXIS ];
	double              n2        = director[ J_AXIS ];
	double              n3        = director[ K_AXIS ];
	double              a00,a01,a02,a03,a04,a05;
	double                  a11,a12,a13,a14,a15;
	double                      a22,a23,a24,a25;
	double                          a33,a34,a35;
	double                              a44,a45;
	double                                  a55;	

	a00 = -4 * n1*n1 * ( 1 - n1*n1 ) * deltaViscosity; 
	a01 =  4 * n1*n1 * n2*n2 * deltaViscosity; 
	a02 =  4 * n1*n1 * n3*n3 * deltaViscosity; 
	a03 =  2 * n1*n2 * (2*n1*n1-1) * deltaViscosity;
	a04 =  2 * n1*n3 * (2*n1*n1-1) * deltaViscosity; 
	a05 =  4 * n1*n1 * n2*n3 * deltaViscosity;
		 
	a11= 4 * n2*n2 * (n2*n2-1) * deltaViscosity; 
	a12= 4 * n2*n2 * n3*n3 * deltaViscosity; 
	a13= 2 * n1*n2 * (2*n2*n2-1) * deltaViscosity; 
	a14= 4 * n1*n2 * n2*n3 * deltaViscosity;  
	a15= 2 * n2*n3 * (2*n2*n2-1) * deltaViscosity;
		
	a22 = 4 * n3*n3 * (n3*n3-1) * deltaViscosity; 
	a23 = 4 * n1*n2 * n3*n3 * deltaViscosity; 
	a24 = 2 * n1*n3 * (2*n3*n3-1) * deltaViscosity;
	a25 = 2 * n2*n3 * (2*n3*n3-1) * deltaViscosity;
	 
	a33 = (4 * n1*n1 * n2*n2 - n1*n1 - n2*n2) * deltaViscosity; 
	a34 = (4 * n1*n1 * n2*n3 - n2*n3) * deltaViscosity; 
	a35 = (4 * n1*n2 * n2*n3 - n1*n3) * deltaViscosity;
	
	a44 = (4 * n1*n1 * n3*n3 - n1*n1 -n3*n3) * deltaViscosity; 
	a45 = (4 * n1*n2 * n3*n3 - n1*n2) * deltaViscosity;

	a55 = (4 * n3*n3 * n2*n2 - n3*n3 - n2*n2) * deltaViscosity;
		
	/* D_{anisotropic} to D */
	D[0][0] += a00 ; D[0][1] += a01 ; D[0][2] += a02 ; D[0][3] += a03 ; D[0][4] += a04 ; D[0][5] += a05 ;
	D[1][0] += a01 ; D[1][1] += a11 ; D[1][2] += a12 ; D[1][3] += a13 ; D[1][4] += a14 ; D[1][5] += a15 ;
	D[2][0] += a02 ; D[2][1] += a12 ; D[2][2] += a22 ; D[2][3] += a23 ; D[2][4] += a24 ; D[2][5] += a25 ;
	D[3][0] += a03 ; D[3][1] += a13 ; D[3][2] += a23 ; D[3][3] += a33 ; D[3][4] += a34 ; D[3][5] += a35 ;
	D[4][0] += a04 ; D[4][1] += a14 ; D[4][2] += a24 ; D[4][3] += a34 ; D[4][4] += a44 ; D[4][5] += a45 ;
	D[5][0] += a05 ; D[5][1] += a15 ; D[5][2] += a25 ; D[5][3] += a35 ; D[5][4] += a45 ; D[5][5] += a55 ;

	self->isDiagonal = False;
}

/*
[B] = [ d/dx,     0  ]
      [    0,  d/dy  ]
      [ d/dy,  d/dx  ]  */
void _ConstitutiveMatrixCartesian2D_Assemble_D_B( void* constitutiveMatrix, double** GNx, Node_Index node_I, double** D_B ){
	ConstitutiveMatrix* self = (ConstitutiveMatrix*) constitutiveMatrix;
	double**            D    = self->matrixData;
	double              d_dx = GNx[ I_AXIS ][ node_I ];
	double              d_dy = GNx[ J_AXIS ][ node_I ];

	if (self->isDiagonal) {
		D_B[0][0] = D[0][0] * d_dx;
		D_B[0][1] = 0.0;
				
		D_B[1][0] = 0.0;
		D_B[1][1] = D[1][1] * d_dy;
				
		D_B[2][0] = D[2][2] * d_dy;
		D_B[2][1] = D[2][2] * d_dx;		
	}
	else {
		D_B[0][0] = D[0][0] * d_dx + D[0][2] * d_dy;
		D_B[0][1] = D[0][1] * d_dy + D[0][2] * d_dx;
				
		D_B[1][0] = D[1][0] * d_dx + D[1][2] * d_dy;
		D_B[1][1] = D[1][1] * d_dy + D[1][2] * d_dx;
				
		D_B[2][0] = D[2][0] * d_dx + D[2][2] * d_dy;
		D_B[2][1] = D[2][1] * d_dy + D[2][2] * d_dx;
	}
}


/*
[B] = [ d/dx,     0,      0  ]
      [    0,  d/dy,      0  ]
      [    0,     0,   d/dx  ]
      [ d/dy,  d/dx,      0  ] 
      [ d/dz,     0,   d/dx  ]
      [    0,  d/dz,   d/dy  ] */
void _ConstitutiveMatrixCartesian3D_Assemble_D_B( void* constitutiveMatrix, double** GNx, Node_Index node_I, double** D_B ){
	ConstitutiveMatrix* self = (ConstitutiveMatrix*) constitutiveMatrix;
	double**            D    = self->matrixData;
	double              d_dx = GNx[ I_AXIS ][ node_I ];
	double              d_dy = GNx[ J_AXIS ][ node_I ];
	double              d_dz = GNx[ K_AXIS ][ node_I ];

	if (self->isDiagonal) {
		D_B[0][0] = D[0][0] * d_dx;
		D_B[0][1] = 0.0;
		D_B[0][2] = 0.0;
		
		D_B[1][0] = 0.0;
		D_B[1][1] = D[1][1] * d_dy;
		D_B[1][2] = 0.0;

		D_B[2][0] = 0.0;
		D_B[2][1] = 0.0;
		D_B[2][2] = D[2][2] * d_dz;

		D_B[3][0] = D[3][3] * d_dy;
		D_B[3][1] = D[3][3] * d_dx;
		D_B[3][2] = 0.0;

		D_B[4][0] = D[4][4] * d_dz;
		D_B[4][1] = 0.0;
		D_B[4][2] = D[4][4] * d_dx;
			
		D_B[5][0] = 0.0;
		D_B[5][1] = D[5][5] * d_dz;
		D_B[5][2] = D[5][5] * d_dy;
	}
	else {
		D_B[0][0] = D[0][0] * d_dx + D[0][3] * d_dy + D[0][4] * d_dz;
		D_B[0][1] = D[0][1] * d_dy + D[0][3] * d_dx + D[0][5] * d_dz;
		D_B[0][2] = D[0][2] * d_dz + D[0][4] * d_dx + D[0][5] * d_dy;
		
		D_B[1][0] = D[1][0] * d_dx + D[1][3] * d_dy + D[1][4] * d_dz;
		D_B[1][1] = D[1][1] * d_dy + D[1][3] * d_dx + D[1][5] * d_dz;
		D_B[1][2] = D[1][2] * d_dz + D[1][4] * d_dx + D[1][5] * d_dy;

		D_B[2][0] = D[2][0] * d_dx + D[2][3] * d_dy + D[2][4] * d_dz;
		D_B[2][1] = D[2][1] * d_dy + D[2][3] * d_dx + D[2][5] * d_dz;
		D_B[2][2] = D[2][2] * d_dz + D[2][4] * d_dx + D[2][5] * d_dy;

		D_B[3][0] = D[3][0] * d_dx + D[3][3] * d_dy + D[3][4] * d_dz;
		D_B[3][1] = D[3][1] * d_dy + D[3][3] * d_dx + D[3][5] * d_dz;
		D_B[3][2] = D[3][2] * d_dz + D[3][4] * d_dx + D[3][5] * d_dy;

		D_B[4][0] = D[4][0] * d_dx + D[4][3] * d_dy + D[4][4] * d_dz;
		D_B[4][1] = D[4][1] * d_dy + D[4][3] * d_dx + D[4][5] * d_dz;
		D_B[4][2] = D[4][2] * d_dz + D[4][4] * d_dx + D[4][5] * d_dy;
			
		D_B[5][0] = D[5][0] * d_dx + D[5][3] * d_dy + D[5][4] * d_dz;
		D_B[5][1] = D[5][1] * d_dy + D[5][3] * d_dx + D[5][5] * d_dz;
		D_B[5][2] = D[5][2] * d_dz + D[5][4] * d_dx + D[5][5] * d_dy;
	}
}

void _ConstitutiveMatrixCartesian2D_CalculateStress( void* constitutiveMatrix, SymmetricTensor strainRate, SymmetricTensor stress ) {
	ConstitutiveMatrix* self = (ConstitutiveMatrix*) constitutiveMatrix;
	double**            D    = self->matrixData;

	if (self->isDiagonal) {
		stress[0] = D[0][0] * strainRate[0];
		stress[1] = D[1][1] * strainRate[1];
		stress[2] = D[2][2] * 2.0 * strainRate[2];
	}
	else {
		stress[0] = D[0][0] * strainRate[0] + D[0][1] * strainRate[1] + D[0][2] * 2.0 * strainRate[2];
		stress[1] = D[1][0] * strainRate[0] + D[1][1] * strainRate[1] + D[1][2] * 2.0 * strainRate[2];
		stress[2] = D[2][0] * strainRate[0] + D[2][1] * strainRate[1] + D[2][2] * 2.0 * strainRate[2];
	}
}



void _ConstitutiveMatrixCartesian3D_CalculateStress( void* constitutiveMatrix, SymmetricTensor strainRate, SymmetricTensor stress ) {
	ConstitutiveMatrix* self = (ConstitutiveMatrix*) constitutiveMatrix;
	double**            D    = self->matrixData;
	          	
	if (self->isDiagonal) {
		stress[0] = D[0][0] * strainRate[0];
		stress[1] = D[1][1] * strainRate[1];
		stress[2] = D[2][2] * strainRate[2];
		stress[3] = D[3][3] * 2.0 * strainRate[3];
		stress[4] = D[4][4] * 2.0 * strainRate[4];
		stress[5] = D[5][5] * 2.0 * strainRate[5];
	}
	else {
		stress[0] = D[0][0] * strainRate[0] + D[0][1] * strainRate[1] + D[0][2] * strainRate[2] 
			+ 2.0 * (D[0][3] * strainRate[3] + D[0][4] * strainRate[4] + D[0][5] * strainRate[5]);

		stress[1] = D[1][0] * strainRate[0] + D[1][1] * strainRate[1] + D[1][2] * strainRate[2] 
			+ 2.0 * (D[1][3] * strainRate[3] + D[1][4] * strainRate[4] + D[1][5] * strainRate[5]);

		stress[2] = D[2][0] * strainRate[0] + D[2][1] * strainRate[1] + D[2][2] * strainRate[2] 
			+ 2.0 * (D[2][3] * strainRate[3] + D[2][4] * strainRate[4] + D[2][5] * strainRate[5]);

		stress[3] = D[3][0] * strainRate[0] + D[3][1] * strainRate[1] + D[3][2] * strainRate[2] 
			+ 2.0 * (D[3][3] * strainRate[3] + D[3][4] * strainRate[4] + D[3][5] * strainRate[5]);

		stress[4] = D[4][0] * strainRate[0] + D[4][1] * strainRate[1] + D[4][2] * strainRate[2] 
			+ 2.0 * (D[4][3] * strainRate[3] + D[4][4] * strainRate[4] + D[4][5] * strainRate[5]);

		stress[5] = D[5][0] * strainRate[0] + D[5][1] * strainRate[1] + D[5][2] * strainRate[2] 
			+ 2.0 * (D[5][3] * strainRate[3] + D[5][4] * strainRate[4] + D[5][5] * strainRate[5]);
	}
}

void ConstitutiveMatrixCartesian_SetupParticleStorage( ConstitutiveMatrixCartesian* self ) {
	/* a function which defines the storage of each particle's constitutive information on the particle, 
	 * should be called before the "Build" phase */

	static int beenHere = 0; /* don't want to do this routine twice */

	IntegrationPointsSwarm* swarm = (IntegrationPointsSwarm*)self->integrationSwarm;
	MaterialPointsSwarm **materialSwarms, *materialSwarm;
	MaterialPoint particle;
	int materialSwarmCount;
	double *cMatrix = NULL;

	if( beenHere ) return; 

	beenHere = 1;

	/* get material swram mapped to the integration points,
	*      * currently only one material point is mapped 26 FEB 09 */
	materialSwarms = IntegrationPointMapper_GetMaterialPointsSwarms( swarm->mapper, &materialSwarmCount );
	assert( materialSwarmCount < 2 );
	materialSwarm = materialSwarms[0];

	/* add extension to material swarm */
	self->storedConstHandle = ExtensionManager_Add( 
	materialSwarm->particleExtensionMgr, 
	self->type, 
	self->rowSize * self->columnSize * sizeof(double) );

	cMatrix = ExtensionManager_Get( materialSwarm->particleExtensionMgr, &particle, self->storedConstHandle );

	if( self->dim == 2 ) {
		/* TODO: clean up this vector logic. The only reson there's an if is because
		*        * of the list of names the must be given as the final arguments to this function.  */ 
		self->storedConstSwarmVar = Swarm_NewVectorVariable( materialSwarm, "ConstitutiveMatrix",
		(ArithPointer)cMatrix - (ArithPointer)&particle,
		Variable_DataType_Double, self->rowSize * self->columnSize,
		"c00", "c01", "c02", "c10", "c11", "c12", "c20", "c21", "c22" );
	} else {
		self->storedConstSwarmVar = Swarm_NewVectorVariable( materialSwarm, "ConstitutiveMatrix",
		(ArithPointer)cMatrix - (ArithPointer)&particle,
		Variable_DataType_Double, self->rowSize * self->columnSize,
		"c00", "c01", "c02", "c03", "c04", "c05",
		"c10", "c11", "c12", "c13", "c14", "c15",
		"c20", "c21", "c22", "c23", "c24", "c25",
		"c30", "c31", "c32", "c33", "c34", "c35",
		"c40", "c41", "c42", "c43", "c44", "c45",
		"c50", "c51", "c52", "c53", "c54", "c55" );
	}

	/* set the storedConstitutive matrix NOT to be checkpointed */
	self->storedConstSwarmVar->isCheckpointedAndReloaded = False;
}
