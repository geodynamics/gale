/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
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
** $Id: RecoveredFeVariable.c 518 2007-10-11 08:07:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include "Underworld/Rheology/Rheology.h"

#include "types.h"
#include "BaseRecoveryFeVar.h"
#include "RecoveredFeVariable.h"
#include "REP_Algorithm.h"

#include <math.h>
#include <assert.h>
#include <string.h>

const Type RecoveredFeVariable_Type = "RecoveredFeVariable";

void* _RecoveredFeVariable_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(RecoveredFeVariable);
	Type                                                      type = RecoveredFeVariable_Type;
	Stg_Class_DeleteFunction*                              _delete = _RecoveredFeVariable_Delete;
	Stg_Class_PrintFunction*                                _print = _RecoveredFeVariable_Print;
	Stg_Class_CopyFunction*                                  _copy = NULL;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _RecoveredFeVariable_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _RecoveredFeVariable_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _RecoveredFeVariable_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _RecoveredFeVariable_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _RecoveredFeVariable_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _RecoveredFeVariable_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType                                         nameAllocationType = (AllocationType)ZERO;
	FieldVariable_InterpolateValueAtFunction*             _interpolateValueAt = ZERO;
	FieldVariable_GetValueFunction*               _getMinGlobalFieldMagnitude = ZERO;
	FieldVariable_GetValueFunction*               _getMaxGlobalFieldMagnitude = ZERO;
	FieldVariable_GetCoordFunction*                  _getMinAndMaxLocalCoords = ZERO;
	FieldVariable_GetCoordFunction*                 _getMinAndMaxGlobalCoords = ZERO;
	FeVariable_InterpolateWithinElementFunction*    _interpolateWithinElement = ZERO;
	FeVariable_GetValueAtNodeFunction*                        _getValueAtNode = ZERO;
	FeVariable_SyncShadowValuesFunc*                        _syncShadowValues = ZERO;

	return (void*)_RecoveredFeVariable_New(  RECOVEREDFEVARIABLE_PASSARGS  );
}

RecoveredFeVariable* _RecoveredFeVariable_New(  RECOVEREDFEVARIABLE_DEFARGS  )
{
	RecoveredFeVariable*		self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(RecoveredFeVariable) );
	self = (RecoveredFeVariable*)
		_BaseRecoveryFeVar_New(  BASERECOVERYFEVAR_PASSARGS  );

	return self;
}

/* --- Virtual Function Implementations --- */
void _RecoveredFeVariable_Delete( void* recFeVariable ) {
	RecoveredFeVariable* self = (RecoveredFeVariable*) recFeVariable;
	_BaseRecoveryFeVar_Delete( self );
}

void _RecoveredFeVariable_Print( void* recFeVariable, Stream* stream ) {
	RecoveredFeVariable* self = (RecoveredFeVariable*) recFeVariable;
	
	/* General info */
	Journal_Printf( stream, "RecoveredFeVariable (ptr): %p\n", self );
	
	/* Print parent */
	_FeVariable_Print( self, stream );
	
	/* RecoveredFeVariable info */
}

void* _RecoveredFeVariable_Copy( const void* recFeVariable, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	abort();
}

void _RecoveredFeVariable_Init( RecoveredFeVariable* self, StiffnessMatrix* stiffnessMatrix, FeVariable* rawPressureField, FeVariable* velGradField, Bool refreshMeshConnectivity, Bool recoverStrain ) {
	/* Setup basic pointers and functionPtrs that are specific to problem spec  */

	self->rawPressureField = rawPressureField;
	self->velGradField = velGradField;
	self->refreshMeshConnectivity = refreshMeshConnectivity;
	self->recoverStrain = recoverStrain;
	self->fieldComponentCount = StGermain_nSymmetricTensorVectorComponents(self->dim);
   self->nonLinearProblem = stiffnessMatrix->isNonLinear;

	/* assign functionPtrs for 2D and 3D here */
	if( self->dim == 2 ) {    
		self->_calcHi = _RecoveredFeVariable_CalcHi2D;
		self->_calcFi = _RecoveredFeVariable_CalcFi2D;
		self->nodesInPatch = 9; /*TODO is static, should be dynamic later */
	} else {
		self->_calcHi =	_RecoveredFeVariable_CalcHi3D;
		self->_calcFi =	_RecoveredFeVariable_CalcFi3D;
		self->nodesInPatch = 27; /*TODO is static, should be dynamic later */
	}

	self->coeffCount = self->orderOfInterpolation * self->fieldComponentCount;
}

void _RecoveredFeVariable_AssignFromXML( void* recFeVariable, Stg_ComponentFactory* cf, void* data ) {
   RecoveredFeVariable*  self        = (RecoveredFeVariable*) recFeVariable;
   StiffnessMatrix*      stiffnessMatrix = NULL;
   FeVariable*           rawPressureField = NULL;
   FeVariable*           velGradField = NULL;
   Bool                  refreshMeshConnectivity, recoverStrain, recoverPressure;

   _BaseRecoveryFeVar_AssignFromXML( self, cf, data );

   /* assume there this always a stiffness matrix called k_matrix */
   stiffnessMatrix = Stg_ComponentFactory_ConstructByName( cf, (Name)"k_matrix", StiffnessMatrix, True, data  ); 

   recoverPressure = Stg_ComponentFactory_GetBool( cf, self->name, (Dictionary_Entry_Key)"recoverPressure", True ); 
   if( recoverPressure == True  ) {
      rawPressureField = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"RawPressureField", FeVariable, True, data );
      /* record the rawField that is required on the repRequired_RawField register */
      if( (unsigned )-1 == Stg_ObjectList_GetIndex( repRequiredRawFields_Reg, rawPressureField->name ) )
         Stg_ObjectList_Append( repRequiredRawFields_Reg, rawPressureField );
   }

   velGradField = Stg_ComponentFactory_ConstructByName( cf, (Name)"VelocityGradientsField", FeVariable, False, data  );

   refreshMeshConnectivity = Stg_ComponentFactory_GetBool( cf, self->name , (Dictionary_Entry_Key)"refreshMeshConnectivity", False  );
   recoverStrain = Stg_ComponentFactory_GetBool( cf, self->name , (Dictionary_Entry_Key)"recoverStrain", False  );

   _RecoveredFeVariable_Init( self, stiffnessMatrix, rawPressureField, velGradField, refreshMeshConnectivity, recoverStrain );
}

void _RecoveredFeVariable_Build( void* recFeVariable, void* data ) {
	RecoveredFeVariable* self = (RecoveredFeVariable*) recFeVariable;

	/* Build the parent class */
	_BaseRecoveryFeVar_Build( self, data );

	/* check and build for pressure field */
	if( self->rawPressureField )
		Stg_Component_Build( self->rawPressureField, data, False );

   Stg_Component_Build( self->velGradField, data, False );
	/* Make static, for current release of code.
	 * When AMR comes this will need to be addressed */
	self->nodesPerEl = (self->rawField->dim == 2) ? 4 : 8;
	if( self->refreshMeshConnectivity == False ) {
		self->elementRep_H = Memory_Alloc_Array( double, 
						FeMesh_GetNodeDomainSize( self->feMesh )*
						self->fieldComponentCount*
						self->dim * self->nodesPerEl * self->orderOfInterpolation, 
						"domain Element REP RHS");
		self->elementRep_F = Memory_Alloc_Array( double, 
						FeMesh_GetNodeDomainSize( self->feMesh )*
						self->fieldComponentCount*
						self->nodesPerEl * self->dim,
						"domain Element REP LHS");
	}
}

void _RecoveredFeVariable_Initialise( void* recFeVariable, void* data ) {
   RecoveredFeVariable*      self = (RecoveredFeVariable*) recFeVariable;

   /* Initialise parent */
   _BaseRecoveryFeVar_Initialise( self, data );

   /* Initialise class specific stuff */
   if( self->rawPressureField ) 
      Stg_Component_Initialise( self->rawPressureField, data, False );
   Stg_Component_Initialise( self->dofLayout, data, False );
   Stg_Component_Initialise( self->velGradField, data, False );

   self->pMatrix = Memory_Alloc_2DArray( double, self->fieldComponentCount, self->orderOfInterpolation , (Name)"P_Matrix" );
   self->CPmat = Memory_Alloc_2DArray( double, self->fieldComponentCount, self->orderOfInterpolation, "CP matrix, for strainRate only");
   self->tmpC  = Memory_Alloc_2DArray( double, self->fieldComponentCount, self->fieldComponentCount, "tmp C matrix, for strainRate only");

   /* this functionPtr could be useful later, JG, 12May09
   * self->_assembleOnParticle = _RecoveredFeVariable_AssembleAtParticle;
   */
   self->_putElIntoProc = _RecoveredFeVariable_PutElIntoProc;
}

void _RecoveredFeVariable_Execute( void* recFeVariable, void* data ) {
	RecoveredFeVariable* self = (RecoveredFeVariable*) recFeVariable; 

	_BaseRecoveryFeVar_Execute( self, data );
}

void _RecoveredFeVariable_Destroy( void* recFeVariable, void* data ) {
	RecoveredFeVariable* self = (RecoveredFeVariable*) recFeVariable;

	Memory_Free( self->elementRep_H );
	Memory_Free( self->elementRep_F );
	Memory_Free( self->pMatrix );
	Memory_Free( self->CPmat );
	Memory_Free( self->tmpC );

	Stg_Component_Destroy( self->velGradField, data, False );

	if( self->rawPressureField )
		Stg_Component_Destroy( self->rawPressureField, data, False );
	
	_BaseRecoveryFeVar_Destroy( self, data );
}

void _RecoveredFeVariable_AssembleAtParticle(
   RecoveredFeVariable* self,
   ConstitutiveMatrix*  constitutiveMatrix,
   IntegrationPoint*    particle,
   int                  lElement_I, 
   double*              globalCoord,
   double**             GNx,
   double               detJac,
   double***            Hi_Mat,
   double**             Fi_Mat )
{
	SymmetricTensor p_StrainRate, p_Stress, theCopy;
	int dim  = self->dim;

        IntegrationPointsSwarm* NNswarm((IntegrationPointsSwarm*)constitutiveMatrix->integrationSwarm);
        IntegrationPoint* NNparticle(particle);
        NearestNeighbor_Replace(&NNswarm,&NNparticle,lElement_I,dim);
        RheologyMaterial* material = (RheologyMaterial*) IntegrationPointsSwarm_GetMaterialOn(NNswarm, NNparticle);
	int dofThatExist = self->fieldComponentCount;
	int order = self->orderOfInterpolation;
	int dof_I;
	double pressure, factor;
	double* pVec   = self->pVec;
	double** pMatrix = self->pMatrix; 
	double** tmpC = self->tmpC; 

	/* integration factor */
	factor = particle->weight * detJac;

	self->_makePoly( globalCoord, pVec );

	/* get raw field values on particles */
	FeVariable_InterpolateWithinElement( self->rawField, lElement_I, particle->xi, p_StrainRate );

   if( constitutiveMatrix ) {
      if( self->nonLinearProblem ) {
         /*use the stored constMatrix to match current velocity and pressure solution */
         ConstitutiveMatrix_GetStoredMatrixOnParticle( constitutiveMatrix, particle, tmpC );
      } else {
         /* build total constitutive matrix - could diagree with current velcoity and pressure solution */
         ConstitutiveMatrix_Assemble( constitutiveMatrix, lElement_I, particle ); 
      }

      if( !self->recoverStrain ) {
         /* if not recovering strain rate we build p_Stress now */
         if( self->nonLinearProblem ) {
            /* perform stress calculation; note, this use the cumulative call because p_Stress */
            if (dim == 2 ) _TmpConstMat2D_CalculateStress( tmpC, p_StrainRate, p_Stress, constitutiveMatrix->isDiagonal );
            else _TmpConstMat3D_CalculateStress( tmpC, p_StrainRate, p_Stress, constitutiveMatrix->isDiagonal );

         } else { ConstitutiveMatrix_CalculateStress( constitutiveMatrix, p_StrainRate, p_Stress ); }
      }
   } else {
      /* if not consitutiveMatrix was given assume viscosity = 1 */
      for( dof_I = 0 ; dof_I < dofThatExist ; dof_I++ ) 
         p_Stress[ dof_I ] = 2 * p_StrainRate[ dof_I ]; 
   }

   if( material->compressible ) { 
   /* Check if material is compressible, if so adjust direct stress components, here it's assumed the
   compressible material has an isotropic viscosity */
      double du_dx, dv_dy, dw_dz;
      double compressibleBit;
      TensorArray p_velGradField;
      FeVariable_InterpolateWithinElement( self->velGradField, lElement_I, particle->xi, p_velGradField );
      du_dx = p_velGradField[0];

      if( dim == 2 ) {
         dv_dy = p_velGradField[2];
         compressibleBit = (-2/3)*constitutiveMatrix->matrixData[0][0]*(du_dx+dv_dy);
         p_Stress[0] = p_Stress[0] + compressibleBit;
         p_Stress[1] = p_Stress[1] + compressibleBit;
      } else {
         dv_dy = p_velGradField[4];
         dw_dz = p_velGradField[8];
         compressibleBit = (-2/3)*constitutiveMatrix->matrixData[0][0]*(du_dx+dv_dy+dw_dz);
         p_Stress[0] = p_Stress[0] + compressibleBit;
         p_Stress[1] = p_Stress[1] + compressibleBit;
         p_Stress[2] = p_Stress[2] + compressibleBit;
      }
   }

	/* If pressure exists, then recovered total pressure by taking the raw pressure field as input */
	if ( self->rawPressureField != NULL ) {	
		FeVariable_InterpolateWithinElement( self->rawPressureField, lElement_I, particle->xi, &pressure );
		p_Stress[0] = p_Stress[0] - pressure; // xx 
		p_Stress[1] = p_Stress[1] - pressure; // yy 
		if( dim == 3 )
			p_Stress[2] = p_Stress[2] - pressure; // zz
	}

#if DEBUG_REP
	p_Stress[0] = 1;
	p_Stress[1] = 1;
	p_Stress[2] = 1;
	p_Stress[3] = 1;
	p_Stress[4] = 1;
	p_Stress[5] = 1;
#endif

	if( self->recoverStrain )
		memcpy( theCopy, p_StrainRate, dofThatExist*sizeof(double) );
	else
		memcpy( theCopy, p_Stress, dofThatExist*sizeof(double) );

	for( dof_I = 0 ; dof_I < dofThatExist ; dof_I++ ) {
    /* now assemble things per dof, see eq. 20 in 
     * B.Boroomand & O.C.Zienkiewicz, 
     * "An Improved REP Recovery and the Effectivity Robustness Test",
     * Int. J. for Numerical Methods in Engineering, vol. 40, pages 3247-3277, 1997. */

		if( self->recoverStrain ) {
			memset( p_StrainRate, 0, dofThatExist*sizeof(double) );
			p_StrainRate[dof_I] = theCopy[dof_I];
		} else {
			/* Setup raw stress tensor for RHS */
			memset( p_Stress, 0, dofThatExist*sizeof(double) );
			p_Stress[dof_I] = theCopy[dof_I];
		}

		/* Setup polynomial matrix for creating LHS */
		ZeroMatrix( pMatrix, dofThatExist, order );
		memcpy( pMatrix[dof_I], pVec, order*sizeof(double) );

		if( self->recoverStrain ) {
			if( !constitutiveMatrix ) {
        /* if no consitutive matrix make it up - assume viscosity = 1 */
				for(int ii = 0 ; ii < dofThatExist; ii++ ) { memset( tmpC[ii], 0, dofThatExist*sizeof(double)); }
				if (dim == 2) { tmpC[0][0] = tmpC[1][1] = tmpC[2][2] = 1; } 
        else { tmpC[0][0] = tmpC[1][1] = tmpC[2][2] = tmpC[3][3] = tmpC[4][4] = tmpC[5][5] = 1; }

			} else {
				/* copy constitutive matrix */
				for(int ii = 0 ; ii < dofThatExist ; ii++)
					memcpy( tmpC[ii], constitutiveMatrix->matrixData[ii], dofThatExist*sizeof(double) ); 
			}
      if( dim == 2 ){
        /* modify stored constitutive matrix to account for symmetry (see CalculateStress routines)  */
        tmpC[0][2] *= 2.0;
        tmpC[1][2] *= 2.0;
        tmpC[2][2] *= 2.0;
      } else {
        tmpC[0][3] *= 2.0;
        tmpC[1][3] *= 2.0;
        tmpC[2][3] *= 2.0;
        tmpC[3][3] *= 2.0;
        tmpC[4][3] *= 2.0;
        tmpC[5][3] *= 2.0;
        tmpC[0][4] *= 2.0;
        tmpC[1][4] *= 2.0;
        tmpC[2][4] *= 2.0;
        tmpC[3][4] *= 2.0;
        tmpC[4][4] *= 2.0;
        tmpC[5][4] *= 2.0;
        tmpC[0][5] *= 2.0;
        tmpC[1][5] *= 2.0;
        tmpC[2][5] *= 2.0;
        tmpC[3][5] *= 2.0;
        tmpC[4][5] *= 2.0;
        tmpC[5][5] *= 2.0;
      }
			/* Multiply in tmpC so that we solve for strainrates */ 
			NonSquareMatrix_MultiplicationByNonSquareMatrix(
					tmpC, dofThatExist, dofThatExist, 
					pMatrix, dofThatExist, order, 
					self->CPmat);
			ConstitutiveMatrix_CalculateStress( constitutiveMatrix, p_StrainRate, p_Stress );
			self->_calcHi( self, GNx, self->CPmat, factor, Hi_Mat[dof_I] );
		} else {
			/* cal Hi, build lefthand Operator */
			self->_calcHi( self, GNx, pMatrix, factor, Hi_Mat[dof_I] );
		}
		/* cal Fi, build righthand Operator */
		self->_calcFi( self, GNx, p_Stress, factor, Fi_Mat[dof_I] );
	}
}

void _RecoveredFeVariable_ZeroHiFi( RecoveredFeVariable* self, double*** elHi_Mat, double** elFi_Mat ) {
	/* Function description:
   * Zeros object Hi and Fi, probably a bit expensive */
	int dim, order, nodesPerEl, dofThatExist, rowsInH, dof_I;
	
	dim = self->dim;
	dofThatExist = self->fieldComponentCount;
	nodesPerEl = self->nodesPerEl;
	order = self->orderOfInterpolation;
	rowsInH = nodesPerEl * dim;

  ZeroMatrix( elFi_Mat, dofThatExist, rowsInH );
	for( dof_I = 0 ; dof_I < dofThatExist ; dof_I++ ) 
		ZeroMatrix( elHi_Mat[dof_I], rowsInH, order );
}


void _RecoveredFeVariable_CalcHi2D( RecoveredFeVariable* self, double** GNx, double** pMatrix, double factor, double** Hi ) {
  /* See eq 24 in Boroomand and Zienk, 1997 */
	int order, dim, row, order_I, node_I, nodesPerEl;
	double d_dx, d_dy;

	nodesPerEl = self->nodesPerEl;
	dim = self->dim;
	order = self->orderOfInterpolation;

	for( node_I = 0 ; node_I < nodesPerEl ; node_I++ ) {
		d_dx = GNx[0][node_I] * factor;
		d_dy = GNx[1][node_I] * factor;
		row  = dim*node_I;
		for( order_I = 0 ; order_I < order ; order_I++ ) {
			Hi[row][order_I]   += d_dx * pMatrix[0][order_I] + d_dy * pMatrix[2][order_I];
			Hi[row+1][order_I] += d_dy * pMatrix[1][order_I] + d_dx * pMatrix[2][order_I];
		}
	}
}

void _RecoveredFeVariable_CalcHi3D( RecoveredFeVariable* self, double** GNx, double** pMatrix, double factor, double** Hi ) {
  /* See eq 24 in Boroomand and Zienk, 1997 */
	int order, dim, row, order_I, node_I, nodesPerEl;
	double d_dx, d_dy, d_dz;

	nodesPerEl = self->nodesPerEl;
	dim = self->dim;
	order = self->orderOfInterpolation;

	for( order_I = 0 ; order_I < order ; order_I++ ) {
		for( node_I = 0 ; node_I < nodesPerEl ; node_I++ ) {
			d_dx = GNx[0][node_I] * factor;
			d_dy = GNx[1][node_I] * factor;
			d_dz = GNx[2][node_I] * factor;
			row  = dim*node_I;
			Hi[row][order_I]   += d_dx * pMatrix[0][order_I] + d_dy * pMatrix[3][order_I] + d_dz*pMatrix[4][order_I];
			Hi[row+1][order_I] += d_dy * pMatrix[1][order_I] + d_dx * pMatrix[3][order_I] + d_dz*pMatrix[5][order_I];
			Hi[row+2][order_I] += d_dz * pMatrix[2][order_I] + d_dx * pMatrix[4][order_I] + d_dy*pMatrix[5][order_I];
		}
	}
}

void _RecoveredFeVariable_CalcFi2D( RecoveredFeVariable* self, double** GNx, double* p_Stress, double factor, double* Fi ) {
  /* See eq 25 in Boroomand and Zienk, 1997 */
	int node_I, nodesPerEl, dim, row;
	double d_dx, d_dy;

	dim = self->dim;
	nodesPerEl = self->nodesPerEl;

	for( node_I = 0 ; node_I < nodesPerEl ; node_I ++ ) {
		d_dx = GNx[0][node_I] * factor;
		d_dy = GNx[1][node_I] * factor;
		row  = dim*node_I;
		/* TODO: Only in 2-D at present */
		Fi[row]   += d_dx * p_Stress[0] + d_dy * p_Stress[2];
		Fi[row+1] += d_dy * p_Stress[1] + d_dx * p_Stress[2];
	}
}

void _RecoveredFeVariable_CalcFi3D( RecoveredFeVariable* self, double** GNx, double* p_Stress, double factor, double* Fi ) {
  /* See eq 25 in Boroomand and Zienk, 1997 */
	int node_I, nodesPerEl, dim, row;
	double d_dx, d_dy, d_dz;

	dim = self->dim;
	nodesPerEl = self->nodesPerEl;

	for( node_I = 0 ; node_I < nodesPerEl ; node_I ++ ) {
		d_dx = GNx[0][node_I] * factor;
		d_dy = GNx[1][node_I] * factor;
		d_dz = GNx[2][node_I] * factor;
		row  = dim*node_I;
		Fi[row]   += d_dx * p_Stress[0] + d_dy * p_Stress[3] + d_dz*p_Stress[4];
		Fi[row+1] += d_dy * p_Stress[1] + d_dx * p_Stress[3] + d_dz*p_Stress[5];
		Fi[row+2] += d_dz * p_Stress[2] + d_dx * p_Stress[4] + d_dy*p_Stress[5];
	}
}

void _RecoveredFeVariable_PutElIntoProc( RecoveredFeVariable* self, int lElement_I, double*** elHi_Mat, double** elFi_Mat ) {
	int keyF, dof_I, row_I, order_I;
	int dim, dofThatExist, order, nodesPerEl, rowsInH;
	double* ptrH = self->elementRep_H;
	double* ptrF = self->elementRep_F;

	dim = self->dim;
	nodesPerEl = self->nodesPerEl;
	dofThatExist = self->fieldComponentCount;
	order = self->orderOfInterpolation;
	rowsInH = nodesPerEl * dim;
	
	for( dof_I = 0 ; dof_I < dofThatExist ; dof_I++ ) {  
		keyF = (lElement_I*dofThatExist*rowsInH) + (dof_I*rowsInH);
		for( row_I = 0 ; row_I < rowsInH ; row_I++ ) {
			ptrF[ keyF + row_I ] =+ elFi_Mat[dof_I][row_I];

			for( order_I = 0 ; order_I < order ; order_I++ )
				/* lElement_I*dofThatExist*rowsInH*order + row_I*order + order_I */
				ptrH[ (keyF + row_I)*order + order_I ] =+ elHi_Mat[dof_I][row_I][order_I];
		}
	}
}

void RecoveredFeVariable_CommunicateHF( RecoveredFeVariable* self ) {
	/* Update all other procs. */
	Sync* sync;
	int chunkSize_H, chunkSize_F, nLocalEl;

	if( self->dim ==2 ) sync = Mesh_GetSync( self->feMesh, MT_FACE );
	else sync = Mesh_GetSync( self->feMesh, MT_VOLUME );

	nLocalEl = FeMesh_GetElementLocalSize( self->feMesh );
	chunkSize_H = self->fieldComponentCount * self->nodesPerEl * self->dim * self->orderOfInterpolation;
	chunkSize_F = self->fieldComponentCount * self->nodesPerEl * self->dim;
	/* Argument Describtion below:
	 * 0) the sync object
	 * 1) start of local address
	 * 2) size of chunk in array
	 * 3) location + size of local information in array
	 * 4) size of chunk to send???
	 * 5) size of chunk to receive??
	 */
	/* Sync of the elementRep_H */
	Sync_SyncArray( sync, self->elementRep_H, chunkSize_H * sizeof(double), 
			self->elementRep_H + nLocalEl * chunkSize_H, chunkSize_H * sizeof(double), 
			chunkSize_H * sizeof(double) );

	/* Sync of the elementRep_H */
	Sync_SyncArray( sync, self->elementRep_F, chunkSize_F * sizeof(double), 
			self->elementRep_F + nLocalEl * chunkSize_F, chunkSize_F * sizeof(double), 
			chunkSize_F * sizeof(double) );
}

void RecoveredFeVariable_SetupWorkSpace( RecoveredFeVariable* self ) {
	int dofThatExist, order, nodesInPatch, rowsInH;

	dofThatExist = self->fieldComponentCount;
	order = self->orderOfInterpolation;
	nodesInPatch = self->nodesInPatch;
	rowsInH = nodesInPatch * self->dim;

	/* will be useful if AMR begins, currently assumes a static mesh */
	self->patch_H    = Memory_Alloc_3DArray( double, dofThatExist, rowsInH, order , (Name)"H_Matrix" );
	self->patch_F    = Memory_Alloc_2DArray( double, dofThatExist, rowsInH, (Name)"F Vector" );
	self->H_t        = Memory_Alloc_2DArray( double, order, rowsInH/*rowsInH, dofThatExist*/, " H transpose");
	self->AMat       = Memory_Alloc_3DArray( double, dofThatExist, order, order , (Name)"AMatrix" );
	self->bVec       = Memory_Alloc_2DArray( double, dofThatExist, order , (Name)"BVector" );
	self->tmpEl_H    = Memory_Alloc_2DArray( double, rowsInH, order , (Name)"tmp Element H" );
	self->tmpEl_F    = Memory_Alloc_Array( double, rowsInH, "tmp Element H");
}

void RecoveredFeVariable_RemoveWorkSpace( RecoveredFeVariable* self ) {
	/* will be useful if AMR begins, currently assumes a static mesh */
	Memory_Free( self->patch_H );
	Memory_Free( self->patch_F );
	Memory_Free( self->H_t );
	Memory_Free( self->AMat );
	Memory_Free( self->bVec );
	Memory_Free( self->tmpEl_H );
	Memory_Free( self->tmpEl_F );
}

void _RecoveredFeVariable_UnpackElement( double* ptrH, double* ptrF, int rowsInH, int order, double **el_H, double *el_F ) {
	int row_I, order_I;
		/* Unpack */
		for( row_I = 0 ; row_I < rowsInH ; row_I++ ) {
			for(order_I = 0 ; order_I < order ; order_I++ ) {
				el_H[row_I][order_I] = ptrH[row_I*order + order_I];
			}
			el_F[row_I]= ptrF[row_I];
		}
}
void RecoveredFeVariable_SolvePatch( RecoveredFeVariable* self, int pNodeID, int *elList, int nEl, LmStruct* lmStruct ) {
	double **patch_H, *patch_F;
	/* Used for unpacking */
	double **el_H = self->tmpEl_H;
	double *el_F = self->tmpEl_F;
	double *ptrH, *ptrF;
        double coeff[50];
	IArray *inc;
	int dim, dofThatExist, order, nodesPerEl, rowsInElH, rowsInPatchH, row, elID;
	int dof_I, el_I, node_I, nodeID, indexOfEntry, order_I;
	int *nodeList;
	int chunkSize_H;
	int chunkSize_F;

	dim = self->dim;
	nodesPerEl = self->nodesPerEl;
	dofThatExist = self->fieldComponentCount;
	order = self->orderOfInterpolation;
	rowsInElH = nodesPerEl * self->dim;
	rowsInPatchH = self->nodesInPatch * self->dim;

	chunkSize_H =  self->nodesPerEl * self->dim * self->orderOfInterpolation;
	chunkSize_F =  self->nodesPerEl * self->dim;

	/* use the inc on the FeVariable */
	inc = self->inc;

	for( dof_I = 0 ; dof_I < dofThatExist ; dof_I++ ) {
			ZeroMatrix( self->patch_H[dof_I], rowsInPatchH, order );
	}
	ZeroMatrix(self->patch_F, dofThatExist, rowsInPatchH );
		
	for( dof_I = 0 ; dof_I < dofThatExist ; dof_I++ ) {
		patch_H = self->patch_H[dof_I];
		patch_F = self->patch_F[dof_I];
		for( el_I = 0 ; el_I < nEl ; el_I++ ) {
			elID = elList[el_I];
			ptrH = &(self->elementRep_H[ chunkSize_H * (elID*dofThatExist + dof_I) ]);
			ptrF = &(self->elementRep_F[ chunkSize_F * (elID*dofThatExist + dof_I) ]);
			_RecoveredFeVariable_UnpackElement( ptrH, ptrF, rowsInElH, order, el_H, el_F );
			/* get nodes in element */
			FeMesh_GetElementNodes( self->feMesh, elID, inc );
			nodesPerEl = IArray_GetSize( inc );
			nodeList = IArray_GetPtr( inc );
			for( node_I = 0 ; node_I < nodesPerEl ; node_I++ ) {
				nodeID = nodeList[node_I];	
				indexOfEntry = _REP_Algorithm_locateInPatchListStruct( lmStruct, nodeID );	
				assert( indexOfEntry != -1 );
				row = indexOfEntry * dim;

				for( order_I = 0 ; order_I < order ; order_I++ ) {
					patch_H[row][order_I]   += el_H[node_I*dim][order_I]; 
					patch_H[row+1][order_I] += el_H[node_I*dim+1][order_I]; 
					if( dim == 3 )
						patch_H[row+2][order_I] += el_H[node_I*dim+2][order_I]; 
				}

				patch_F[row] += el_F[node_I*dim];
				patch_F[row+1] += el_F[node_I*dim+1];
				if( dim == 3 ) 
					patch_F[row+2] += el_F[node_I*dim+2];
			}
		}
		/* Now the actual maths */
		NonSquareMatrix_Transpose( patch_H, rowsInPatchH, order, self->H_t );
	
		/* 5) Construct A_Matrix = Ht_Matrix . H_Matrix */
		NonSquareMatrix_MultiplicationByNonSquareMatrix( self->H_t, order, rowsInPatchH, 
																patch_H, rowsInPatchH,  order, 
																self->AMat[dof_I] );
	
		/* 6) Construct bMatrix = Ht_Matrix . F_Vector */
		NonSquareMatrix_MatrixVectorMultiplication( self->H_t, order, rowsInPatchH,
																patch_F, rowsInPatchH,
																self->bVec[dof_I] );

		_REP_Algorithm_Solver(self->AMat[dof_I], self->bVec[dof_I], order);
		memcpy( &coeff[order*dof_I], self->bVec[dof_I], sizeof(double)*order );
	}
	FeVariable_SetValueAtNode( self, pNodeID , coeff );
}

void _TmpConstMat2D_CalculateStress( double** D, SymmetricTensor strainRate, SymmetricTensor stress, Bool isDiagonal ) {
	if (isDiagonal) {
		stress[0] = D[0][0] * strainRate[0];
		stress[1] = D[1][1] * strainRate[1];
		stress[2] = D[2][2] * 2.0 * strainRate[2];
	} else {
		stress[0] = D[0][0] * strainRate[0] + D[0][1] * strainRate[1] + D[0][2] * 2.0 * strainRate[2];
		stress[1] = D[1][0] * strainRate[0] + D[1][1] * strainRate[1] + D[1][2] * 2.0 * strainRate[2];
		stress[2] = D[2][0] * strainRate[0] + D[2][1] * strainRate[1] + D[2][2] * 2.0 * strainRate[2];
	}
}

void _TmpConstMat3D_CalculateStress( double** D, SymmetricTensor strainRate, SymmetricTensor stress, Bool isDiagonal ) {
	          	
	if (isDiagonal) {
		stress[0] = D[0][0] * strainRate[0];
		stress[1] = D[1][1] * strainRate[1];
		stress[2] = D[2][2] * strainRate[2];
		stress[3] = D[3][3] * 2.0 * strainRate[3];
		stress[4] = D[4][4] * 2.0 * strainRate[4];
		stress[5] = D[5][5] * 2.0 * strainRate[5];
	} else {
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


