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
** $Id:  $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <StgFEM/FrequentOutput/FrequentOutput.h>

#include "Plateness.h"

const Type Underworld_Plateness_Type = "Underworld_Plateness";

void _Underworld_Plateness_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	UnderworldContext* context;

	context = Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data ); 

	/* Add functions to entry points */
	Underworld_Plateness_Setup( context );
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput, Underworld_Plateness_Output );
}

void _Underworld_Plateness_Build( void* component, void* data ) {
	Underworld_Plateness*	self = (Underworld_Plateness*)component;

	Stg_Component_Build( self->reducedStrainRateFieldInvariantRoot, data, False );

}

void* _Underworld_Plateness_DefaultNew( Name name ) {
	return _Codelet_New(
		sizeof(Underworld_Plateness),
		Underworld_Plateness_Type,
		_Codelet_Delete,
		_Codelet_Print,
		_Codelet_Copy,
		_Underworld_Plateness_DefaultNew,
		_Underworld_Plateness_Construct,
		_Underworld_Plateness_Build,
		_Codelet_Initialise,
		_Codelet_Execute,
		_Codelet_Destroy,
		name );
}

Index Underworld_Plateness_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Underworld_Plateness_Type, "0", _Underworld_Plateness_DefaultNew );
}

void Underworld_Plateness_Setup( UnderworldContext* context ) {
	FieldVariable_Register*              fV_Register               = context->fieldVariable_Register;
	FieldVariable*                       strainRateField;
	Func_Ptr                             _carryOut;
	Dof_Index                            resultDofs;
	Dof_Index                            operandDofs;
	Index                                numberOfOperands;
	Operator*                            ownOperator;
	Dimension_Index                      dim;
	
	Underworld_Plateness* self;

	self = (Underworld_Plateness*)LiveComponentRegister_Get(
					context->CF->LCRegister,
					Underworld_Plateness_Type );
	
	StgFEM_FrequentOutput_PrintString( context, "Plateness" );

	Journal_Firewall( 
			context->gaussSwarm != NULL, 
			Underworld_Error,
			"Cannot find gauss swarm. Cannot use %s.\n", CURR_MODULE_NAME );

	/* Create Some FeVariables to determine plateness */
	strainRateField  = FieldVariable_Register_GetByName( fV_Register, "StrainRateField" );

	dim = strainRateField->dim;
	/* setup required parameters to create new operate */
	resultDofs       = 1 ; //( dim == 2 ? 1 : 3 );
	numberOfOperands = 1 ;
	operandDofs      = ( dim == 2 ? 3 : 6 );
	_carryOut        = ( dim == 2 ? Underworld_Plateness_SymmetricTensor_LowerDimension_InvariantRoot_2d : Underworld_Plateness_SymmetricTensor_LowerDimension_InvariantRoot_3d );
	
	ownOperator = Operator_New( "Plateness_SymmetricTensor_LowerDimension_InvariantRoot", _carryOut, numberOfOperands, operandDofs, resultDofs, dim );

	self->reducedStrainRateFieldInvariantRoot = OperatorFeVariable_NewUnary_OwnOperator(
			"ReducedStrainRateFieldInvariantRoot",
			strainRateField, 
			ownOperator );

	/* Add the variables to register so we can checkpoint & examine if necessary */
	FieldVariable_Register_Add( fV_Register, self->reducedStrainRateFieldInvariantRoot );
}

void Underworld_Plateness_Output( UnderworldContext* context ) {
	Underworld_Plateness* self;
	double   minValue;
	double   maxValue;
	double   totalIntegral;
	double   integralSoFar;
	double   weightSoFar;
	double   weightSoFar2;
	double   *min, *max;
	double   **value_weight_Matrix;
	double   **value_weight_Matrix_Global;
	int ii;
	int jj;
	
	self = (Underworld_Plateness*)LiveComponentRegister_Get(
				context->CF->LCRegister,
				Underworld_Plateness_Type );

	min = Memory_Alloc_Array_Unnamed( double, Mesh_GetDimSize( self->reducedStrainRateFieldInvariantRoot->feMesh ) );
	max = Memory_Alloc_Array_Unnamed( double, Mesh_GetDimSize( self->reducedStrainRateFieldInvariantRoot->feMesh ) );

	Mesh_GetGlobalCoordRange( self->reducedStrainRateFieldInvariantRoot->feMesh, min, max );
	
	//minValue = self->reducedStrainRateFieldInvariantRoot->_getMinGlobalFieldMagnitude( self->reducedStrainRateFieldInvariantRoot );
	//maxValue = self->reducedStrainRateFieldInvariantRoot->_getMaxGlobalFieldMagnitude( self->reducedStrainRateFieldInvariantRoot );

	value_weight_Matrix        = Memory_Alloc_2DArray( double, 2, 100 , "value_weight_Matrix");
	value_weight_Matrix_Global = Memory_Alloc_2DArray( double, 2, 100 , "value_weight_Matrix_Global");
	
	for(ii=0 ; ii<100; ii++) for(jj=0 ; jj<2; jj++)	value_weight_Matrix[jj][ii] = 0;
	for(ii=0 ; ii<100; ii++) for(jj=0 ; jj<2; jj++)	value_weight_Matrix_Global[jj][ii] = 0;
	
	totalIntegral = Plateness_IntegratePlane( self->reducedStrainRateFieldInvariantRoot, J_AXIS, max[ J_AXIS ], value_weight_Matrix); //, minValue, maxValue );
	
	MPI_Allreduce( value_weight_Matrix[0], value_weight_Matrix_Global[0], 100, MPI_DOUBLE, MPI_SUM, context->communicator );
	MPI_Allreduce( value_weight_Matrix[1], value_weight_Matrix_Global[1], 100, MPI_DOUBLE, MPI_SUM, context->communicator );
	
	integralSoFar = 0.;
	weightSoFar   = 0.;
	for(ii = 99; ii >= 0 ; ii--){
		integralSoFar += value_weight_Matrix_Global[0][ii];
		weightSoFar   += value_weight_Matrix_Global[1][ii];
		if (integralSoFar > 0.8*totalIntegral) break;
	}

	integralSoFar = 0.;
	weightSoFar2   = 0.;
	for(ii = 99; ii >= 0 ; ii--){
		integralSoFar += value_weight_Matrix_Global[0][ii];
		weightSoFar2   += value_weight_Matrix_Global[1][ii];
	}
	
	self->plateness = weightSoFar/weightSoFar2;
	//printf("integral count, totalintegral, totalweight, plateness = %f, %f, %f, %f \n", integralSoFar, totalIntegral, weightSoFar2, self->plateness); 
	StgFEM_FrequentOutput_PrintValue( context, self->plateness );
	
	FreeArray( min );
	FreeArray( max );
	FreeArray( value_weight_Matrix );
	FreeArray( value_weight_Matrix_Global );

}

void Underworld_Plateness_SymmetricTensor_LowerDimension_InvariantRoot_2d( void* operator, double* operand0, double* result ) {
	Operator* self = (Operator*) operator;
	double    temp;
	
	Operator_FirewallUnary( self );
	Operator_FirewallResultDofs( self, self->operandDofs );

	temp      = operand0[0];
	result[0] = fabs( temp );
}

void Underworld_Plateness_SymmetricTensor_LowerDimension_InvariantRoot_3d( void* operator, double* operand0, double* result ) {
	Operator* self = (Operator*) operator;
	SymmetricTensor temp;
	double temp2;
	
	Operator_FirewallUnary( self );
	Operator_FirewallResultDofs( self, self->operandDofs );

	temp[0] = operand0[0];
	temp[1] = operand0[1];
	temp[2] = operand0[6];
	
	temp2 = SymmetricTensor_2ndInvariant( temp, 2 );
	
	*result = sqrt( temp2 );
}


double Plateness_IntegratePlane( void* feVariable, Axis planeAxis, double planeHeight, double** value_weight_Matrix ){//, double minValue, double maxValue ) {
	FeVariable*                self               = (FeVariable*)         feVariable;
	IJK                        planeIJK;
	Element_LocalIndex         lElement_I;
	Element_GlobalIndex        gElement_I;
	Element_LocalIndex         elementLocalCount  = FeMesh_GetElementLocalSize( self->feMesh );
	Axis                       aAxis              = ( planeAxis == I_AXIS ? J_AXIS : I_AXIS );
	Axis                       bAxis              = ( planeAxis == K_AXIS ? J_AXIS : K_AXIS );
	double                     integral;
	/* Swarm Stuff */
	Swarm*                     tmpSwarm;
	Bool                       dimExists[]        = { False, False, False };
	ExtensionManager_Register* extensionMgr_Register;
	SingleCellLayout*          singleCellLayout;
	GaussParticleLayout*       gaussParticleLayout;
	Particle_Index             lParticle_I;
	IntegrationPoint*          particle;
	/* Plane location stuff */
	double                     storedXi_J_AXIS;
	Coord                      planeCoord;
	double                     planeXi           = -1;
	double                     planeXiGlobal;
	Index                      planeLayer        = 0;
	Index                      planeLayerGlobal;
	Particle_InCellIndex       particlesPerDim[] = {2,2,2};

	/* Find Elements which plane cuts through */
	memcpy( planeCoord, Mesh_GetVertex( self->feMesh, 0 ), sizeof( Coord ) );
	planeCoord[ planeAxis ] = planeHeight;

	if( Mesh_Algorithms_SearchElements( self->feMesh->algorithms, planeCoord, &lElement_I ) && 
	    lElement_I < elementLocalCount )
	{
		Coord		planeXiCoord;

		gElement_I = FeMesh_ElementDomainToGlobal( self->feMesh, lElement_I );
		RegularMeshUtils_Element_1DTo3D( self->feMesh, gElement_I, planeIJK );
		planeLayer = planeIJK[ planeAxis ];
		
		/* Find Local Coordinate of plane */
		FeMesh_CoordGlobalToLocal( self->feMesh, lElement_I, planeCoord, planeXiCoord );
		planeXi = planeXiCoord[ planeAxis ];
	}
	
	/* Should be broadcast */
	MPI_Allreduce( &planeXi,    &planeXiGlobal, 1, MPI_DOUBLE, MPI_MAX, self->communicator );
	MPI_Allreduce( &planeLayer, &planeLayerGlobal, 1, MPI_UNSIGNED, MPI_MAX, self->communicator );

	/* Create Swarm in plane */
	extensionMgr_Register = ExtensionManager_Register_New();
	dimExists[ aAxis ] = True;
	if (self->dim == 3)
		dimExists[ bAxis ] = True;
	
	singleCellLayout = SingleCellLayout_New( "cellLayout", dimExists, NULL, NULL );
	particlesPerDim[ planeAxis ] = 1;
	gaussParticleLayout = GaussParticleLayout_New( "particleLayout", self->dim - 1, particlesPerDim );
	tmpSwarm = Swarm_New( 
			"tmpgaussSwarm",
			singleCellLayout, 
			gaussParticleLayout,
			self->dim,
			sizeof(IntegrationPoint), 
			extensionMgr_Register, 
			NULL,
			self->communicator,
		        NULL	);
	Stg_Component_Build( tmpSwarm, NULL, False );

	/* Change Positions of the particles */
	Stg_Component_Initialise( tmpSwarm, NULL, False );
	for ( lParticle_I = 0 ; lParticle_I < tmpSwarm->particleLocalCount ; lParticle_I++ ) {
		particle = (IntegrationPoint*) Swarm_ParticleAt( tmpSwarm, lParticle_I );

		storedXi_J_AXIS = particle->xi[ J_AXIS ];
		particle->xi[ aAxis ]     = particle->xi[ I_AXIS ];
		particle->xi[ bAxis ]     = storedXi_J_AXIS;
		particle->xi[ planeAxis ] = planeXiGlobal;
	}
	
	integral = Plateness_IntegrateLayer_AxisIndependent( self, tmpSwarm, planeAxis, planeLayerGlobal, 
			self->dim - 1, aAxis, bAxis, planeAxis, value_weight_Matrix); //, minValue, maxValue );

	/* Delete */
	Stg_Class_Delete( tmpSwarm );
	Stg_Class_Delete( gaussParticleLayout );
	Stg_Class_Delete( singleCellLayout );
	Stg_Class_Delete( extensionMgr_Register );
	
	return integral;
}

double Plateness_IntegrateLayer_AxisIndependent( 
		void* feVariable, void* _swarm,
		Axis layerAxis, Index layerIndex, Dimension_Index dim, 
		Axis axis0, Axis axis1, Axis axis2, double** value_weight_Matrix) //, double minValue, double maxValue ) 
{ 
	FeVariable*                self               = (FeVariable*)         feVariable;
	Swarm*                     swarm              = (Swarm*)              _swarm;
	Element_LocalIndex         lElement_I;
	Element_GlobalIndex        gElement_I;
	IJK                        elementIJK;
	double                     elementIntegral;
	double                     integral;
	double                     integralGlobal;
	double                     minValue;
	double                     maxValue;
	IntegrationPoint*          particle;
	Cell_LocalIndex            cell_I;
	Particle_InCellIndex       cParticle_I;
	Particle_InCellIndex       cellParticleCount;
	double                     value;
	
	Journal_DPrintf( self->debug, "In %s() for FeVariable \"%s\":\n", __func__, self->name );

	/* Initialise Sumation of Integral */
	integral = 0.0;

	Stream_Indent( self->debug );
	
	/* Loop over all particles in element to get minValue and maxValue*/
	minValue = 1.e300;
	maxValue = 0.;
	for ( gElement_I = 0 ; gElement_I < FeMesh_GetElementGlobalSize( self->feMesh ); gElement_I++ ) {
		RegularMeshUtils_Element_1DTo3D( self->feMesh, gElement_I, elementIJK );

		/* Check if element is in layer plane */
		if ( elementIJK[ layerAxis ] != layerIndex )
			continue;

		/* Check if element is local */
		if( !FeMesh_ElementGlobalToDomain( self->feMesh, gElement_I, &lElement_I ) || 
		    lElement_I >= FeMesh_GetElementLocalSize( self->feMesh ) )
		{
			continue;
		}

		/* Determine number of particles in element */
		cell_I = CellLayout_MapElementIdToCellId( swarm->cellLayout, lElement_I );
		cellParticleCount = swarm->cellParticleCountTbl[ cell_I ];

		for( cParticle_I = 0 ; cParticle_I < cellParticleCount ; cParticle_I++ ) {
			/* Get Pointer to particle */
			particle = (IntegrationPoint*) Swarm_ParticleInCellAt( swarm, cell_I, cParticle_I );
	
			/* Interpolate Value of Field at Particle */
			FeVariable_InterpolateWithinElement( feVariable, lElement_I, particle->xi, &value );
			
			if(value < minValue) minValue = value;
			if(value > maxValue) maxValue = value;
	
		}
		
	}
	
	for ( gElement_I = 0 ; gElement_I < FeMesh_GetElementGlobalSize( self->feMesh ); gElement_I++ ) {
		RegularMeshUtils_Element_1DTo3D( self->feMesh, gElement_I, elementIJK );

		/* Check if element is in layer plane */
		if ( elementIJK[ layerAxis ] != layerIndex )
			continue;

		/* Check if element is local */
		if( !FeMesh_ElementGlobalToDomain( self->feMesh, gElement_I, &lElement_I ) || 
		    lElement_I >= FeMesh_GetElementLocalSize( self->feMesh ) )
		{
			continue;
		}

		elementIntegral = Plateness_IntegrateElement_AxisIndependent( self, swarm, lElement_I, dim, axis0, axis1, axis2,  value_weight_Matrix, minValue, maxValue );
		Journal_DPrintfL( self->debug, 2, "Integral of element %d was %f\n", lElement_I, elementIntegral );
		integral += elementIntegral;
	}
	Stream_UnIndent( self->debug );


	/* Gather and sum integrals from other processors */
	MPI_Allreduce( &integral, &integralGlobal, 1, MPI_DOUBLE, MPI_SUM, self->communicator );

	Journal_DPrintf( self->debug, "Calculated global integral of layer %d in Axis %d was %f\n", layerIndex, layerAxis, integralGlobal );
	return integralGlobal;
}

double Plateness_IntegrateElement_AxisIndependent( 
		void* feVariable, void* _swarm, 
		Element_DomainIndex dElement_I, Dimension_Index dim, 
		Axis axis0, Axis axis1, Axis axis2, double** value_weight_Matrix, double minValue, double maxValue ) 
{
	FeVariable*          self               = (FeVariable*)         feVariable;
	Swarm*               swarm              = (Swarm*)              _swarm;
	FeMesh*			feMesh             = self->feMesh;
	FeMesh*			mesh;
	ElementType*         elementType;
	Cell_LocalIndex      cell_I;
	Particle_InCellIndex cParticle_I;
	Particle_InCellIndex cellParticleCount;
	IntegrationPoint*    particle;
	double               detJac;
	double               integral;
	double               value;
	double               whereBlock;
	long                 whichBlock;
	
	/* Initialise Summation of Integral */
	integral = 0.0;

	/* Use feVariable's mesh as geometry mesh if one isn't passed in */
	if( Stg_Class_IsInstance( feMesh->algorithms, Mesh_CentroidAlgorithms_Type ) )
		mesh = (FeMesh*)((Mesh_CentroidAlgorithms*)feMesh->algorithms)->elMesh;
	else
		mesh = feMesh;
	elementType = FeMesh_GetElementType( mesh, dElement_I );

	/* Determine number of particles in element */
	cell_I = CellLayout_MapElementIdToCellId( swarm->cellLayout, dElement_I );
	cellParticleCount = swarm->cellParticleCountTbl[ cell_I ];

	/* Loop over all particles in element to get integral*/
	for( cParticle_I = 0 ; cParticle_I < cellParticleCount ; cParticle_I++ ) {
		/* Get Pointer to particle */
		particle = (IntegrationPoint*) Swarm_ParticleInCellAt( swarm, cell_I, cParticle_I );

		/* Interpolate Value of Field at Particle */
		FeVariable_InterpolateWithinElement( feVariable, dElement_I, particle->xi, &value );

		Journal_DPrintfL( self->debug, 3, "%s: Integrating element %d - particle %d - Value = %g\n", self->name, dElement_I, cParticle_I, value );

		/* Calculate Determinant of Jacobian */
		detJac = ElementType_JacobianDeterminant_AxisIndependent( 
				elementType, mesh, dElement_I, particle->xi, dim, axis0, axis1, axis2 );
		
		whereBlock = -minValue*100/(maxValue-minValue) + value*100/(maxValue-minValue);
		whichBlock = (long)whereBlock;
		if(whichBlock < 0 ) whichBlock = 0;
		if(whichBlock > 99) whichBlock = 99;
		
		value_weight_Matrix[0][whichBlock] += detJac * particle->weight * value;
		value_weight_Matrix[1][whichBlock] += detJac * particle->weight;
		
		/* Sum Integral */
		integral += detJac * particle->weight * value;
	}
	
	return integral;
}
