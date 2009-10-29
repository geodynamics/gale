/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**	Computational Infrastructure for Geodynamics - http://www.geodynamics.org
**
** Contributors:
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
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
** $Id: GradientStiffnessMatrixTerm.c 964 2007-10-11 08:03:06Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/Discretisation/Discretisation.h>
#include <StgFEM/SLE/SLE.h>

#include "types.h"
#include "GradientStiffnessMatrixTerm.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class */
const Type GradientStiffnessMatrixTerm_Type = "GradientStiffnessMatrixTerm";

GradientStiffnessMatrixTerm* GradientStiffnessMatrixTerm_New( 
		Name                                                name,
		StiffnessMatrix*                                    stiffnessMatrix,
		Swarm*                                              integrationSwarm )
{
	GradientStiffnessMatrixTerm* self = (GradientStiffnessMatrixTerm*) _GradientStiffnessMatrixTerm_DefaultNew( name );

	self->isConstructed = True;
	_StiffnessMatrixTerm_Init( self, stiffnessMatrix, integrationSwarm, NULL );
	_GradientStiffnessMatrixTerm_Init( self );

	return self;
}

/* Creation implementation / Virtual constructor */
GradientStiffnessMatrixTerm* _GradientStiffnessMatrixTerm_New( 
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
		Name                                                name )
{
	GradientStiffnessMatrixTerm* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(GradientStiffnessMatrixTerm) );
	self = (GradientStiffnessMatrixTerm*) _StiffnessMatrixTerm_New( 
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
	
	/* Virtual info */
	
	return self;
}

void _GradientStiffnessMatrixTerm_Init( 
		GradientStiffnessMatrixTerm*                                    self )
{
	self->max_nElNodes_col = 0;
	self->Ni_col = NULL;
}

void _GradientStiffnessMatrixTerm_Delete( void* matrixTerm ) {
	GradientStiffnessMatrixTerm* self = (GradientStiffnessMatrixTerm*)matrixTerm;

	_StiffnessMatrixTerm_Delete( self );
}

void _GradientStiffnessMatrixTerm_Print( void* matrixTerm, Stream* stream ) {
	GradientStiffnessMatrixTerm* self = (GradientStiffnessMatrixTerm*)matrixTerm;
	
	_StiffnessMatrixTerm_Print( self, stream );

	/* General info */
}

void* _GradientStiffnessMatrixTerm_DefaultNew( Name name ) {
	return (void*)_GradientStiffnessMatrixTerm_New( 
		sizeof(GradientStiffnessMatrixTerm), 
		GradientStiffnessMatrixTerm_Type,
		_GradientStiffnessMatrixTerm_Delete,
		_GradientStiffnessMatrixTerm_Print,
		NULL,
		_GradientStiffnessMatrixTerm_DefaultNew,
		_GradientStiffnessMatrixTerm_Construct,
		_GradientStiffnessMatrixTerm_Build,
		_GradientStiffnessMatrixTerm_Initialise,
		_GradientStiffnessMatrixTerm_Execute,
		_GradientStiffnessMatrixTerm_Destroy,
		_GradientStiffnessMatrixTerm_AssembleElement,
		name );
}

void _GradientStiffnessMatrixTerm_Construct( void* matrixTerm, Stg_ComponentFactory* cf, void* data ) {
	GradientStiffnessMatrixTerm* self = (GradientStiffnessMatrixTerm*)matrixTerm;

	/* Construct Parent */
	_StiffnessMatrixTerm_Construct( self, cf, data );

	_GradientStiffnessMatrixTerm_Init( self );
}

void _GradientStiffnessMatrixTerm_Build( void* matrixTerm, void* data ) {
	GradientStiffnessMatrixTerm* self = (GradientStiffnessMatrixTerm*)matrixTerm;

	_StiffnessMatrixTerm_Build( self, data );
}

void _GradientStiffnessMatrixTerm_Initialise( void* matrixTerm, void* data ) {
	GradientStiffnessMatrixTerm* self = (GradientStiffnessMatrixTerm*)matrixTerm;

	_StiffnessMatrixTerm_Initialise( self, data );
}

void _GradientStiffnessMatrixTerm_Execute( void* matrixTerm, void* data ) {
	_StiffnessMatrixTerm_Execute( matrixTerm, data );
}

void _GradientStiffnessMatrixTerm_Destroy( void* matrixTerm, void* data ) {
	GradientStiffnessMatrixTerm* self = (GradientStiffnessMatrixTerm*)matrixTerm;

	_StiffnessMatrixTerm_Destroy( matrixTerm, data );

	Memory_Free( self->Ni_col );
}


void _GradientStiffnessMatrixTerm_AssembleElement( 
		void*                                              matrixTerm,
		StiffnessMatrix*                                   stiffnessMatrix, 
		Element_LocalIndex                                 lElement_I, 
		SystemLinearEquations*                             sle,
		FiniteElementContext*                              context,
		double**                                           elStiffMat ) 
{
	GradientStiffnessMatrixTerm*        self         = Stg_CheckType( matrixTerm, GradientStiffnessMatrixTerm );
	Swarm*                              swarm        = self->integrationSwarm;
	FeVariable*                         variable_row = stiffnessMatrix->rowVariable;
	FeVariable*                         variable_col = stiffnessMatrix->columnVariable;
	Dimension_Index                     dim          = stiffnessMatrix->dim;
	double*                             xi;
	double                              weight;
	Particle_InCellIndex                cParticle_I, cellParticleCount;
	Node_ElementLocalIndex              nodesPerEl_row;
	Node_ElementLocalIndex              nodesPerEl_col;	
	Dof_Index                           totalDofsThisElement_row, totalDofsThisElement_col;
	
	Dof_Index                           dofPerNode_row, dofPerNode_col;
	Index                               row, col; /* Indices into the stiffness matrix */
	Node_ElementLocalIndex              rowNode_I;
	Node_ElementLocalIndex              colNode_I;
	Dof_Index                           rowDof_I, colDof_I;
	double**                            GNx_row;
	double*                             Ni_col;
	double                              detJac;
	IntegrationPoint*                   currIntegrationPoint;
	
	Cell_Index                          cell_I;
	ElementType*                        elementType_row;
	ElementType*                        elementType_col;
	
	/* Set the element type */
	elementType_row = FeMesh_GetElementType( variable_row->feMesh, lElement_I );
	nodesPerEl_row = elementType_row->nodeCount;
	
	elementType_col = FeMesh_GetElementType( variable_col->feMesh, lElement_I );
	nodesPerEl_col = elementType_col->nodeCount;
		
	dofPerNode_row = dim;	/* velocity */
	dofPerNode_col = 1;	/* pressure */
	
	totalDofsThisElement_row = nodesPerEl_row * dofPerNode_row;
	totalDofsThisElement_col = nodesPerEl_col * dofPerNode_col;
	
	if( nodesPerEl_row > self->max_nElNodes ) {
		/* reallocate */
		self->GNx = ReallocArray2D( self->GNx, double, dim, nodesPerEl_row );
		self->max_nElNodes = nodesPerEl_row;
	}
	GNx_row = self->GNx;

	if( nodesPerEl_col > self->max_nElNodes_col ) {
   /* allocate */
		self->Ni_col = ReallocArray( self->Ni_col, double, nodesPerEl_col );
		self->max_nElNodes_col = nodesPerEl_col;
	}
	Ni_col = self->Ni_col;
	
	cell_I = CellLayout_MapElementIdToCellId( swarm->cellLayout, lElement_I );
	cellParticleCount = swarm->cellParticleCountTbl[ cell_I ];
	
	for ( cParticle_I = 0 ; cParticle_I < cellParticleCount ; cParticle_I++ ) {
		currIntegrationPoint = (IntegrationPoint*)Swarm_ParticleInCellAt( swarm, cell_I, cParticle_I );
		xi = currIntegrationPoint->xi;
		weight = currIntegrationPoint->weight;
		
		/* get shape function derivs for the row (ie. velocity) */
		ElementType_ShapeFunctionsGlobalDerivs( 
				elementType_row,
				variable_row->feMesh, lElement_I,
				xi, dim, &detJac, GNx_row );
		
		/* get the shape functions for the col. (ie. pressure) */
		ElementType_EvaluateShapeFunctionsAt( elementType_col, xi, Ni_col );
		
		/* build stiffness matrix */
		for ( rowNode_I = 0; rowNode_I < nodesPerEl_row ; rowNode_I++) {  	
			for( rowDof_I=0; rowDof_I<dofPerNode_row; rowDof_I++) {	
				row = rowNode_I*dofPerNode_row + rowDof_I;

				for (colNode_I = 0; colNode_I < nodesPerEl_col; colNode_I++ ) { 
					for( colDof_I=0; colDof_I<dofPerNode_col; colDof_I++) {		
						col = colNode_I*dofPerNode_col + colDof_I;
						
						elStiffMat[row][col] +=
							+ weight * ( - detJac ) * ( GNx_row[rowDof_I][rowNode_I] * Ni_col[colNode_I] );
					}
				}		
			}
		}
	}
}
