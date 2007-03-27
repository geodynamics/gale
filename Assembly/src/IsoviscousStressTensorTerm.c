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
** $Id: IsoviscousStressTensorTerm.c 733 2007-02-07 00:55:26Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/Discretisation/Discretisation.h>
#include <StgFEM/SLE/SLE.h>

#include "types.h"
#include "IsoviscousStressTensorTerm.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class */
const Type IsoviscousStressTensorTerm_Type = "IsoviscousStressTensorTerm";

IsoviscousStressTensorTerm* IsoviscousStressTensorTerm_New( 
		Name                                                name,
		StiffnessMatrix*                                    stiffnessMatrix,
		Swarm*                                              integrationSwarm,
		double                                              viscosity )
{
	IsoviscousStressTensorTerm* self = (IsoviscousStressTensorTerm*) _IsoviscousStressTensorTerm_DefaultNew( name );

	IsoviscousStressTensorTerm_InitAll( 
			self,
			stiffnessMatrix,
			integrationSwarm,
			viscosity );

	return self;
}

/* Creation implementation / Virtual constructor */
IsoviscousStressTensorTerm* _IsoviscousStressTensorTerm_New( 
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
	IsoviscousStressTensorTerm* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(IsoviscousStressTensorTerm) );
	self = (IsoviscousStressTensorTerm*) _StiffnessMatrixTerm_New( 
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

void _IsoviscousStressTensorTerm_Init( 
		IsoviscousStressTensorTerm*                         self,
		double                                              viscosity )
{
	self->viscosity = viscosity;
}

void IsoviscousStressTensorTerm_InitAll( 
		void*                                               matrixTerm,
		StiffnessMatrix*                                    stiffnessMatrix,
		Swarm*                                              integrationSwarm,
		double                                              viscosity )
{
	IsoviscousStressTensorTerm* self = (IsoviscousStressTensorTerm*) matrixTerm;

	StiffnessMatrixTerm_InitAll( self, stiffnessMatrix, integrationSwarm, NULL );
	_IsoviscousStressTensorTerm_Init( self, viscosity );
}

void _IsoviscousStressTensorTerm_Delete( void* matrixTerm ) {
	IsoviscousStressTensorTerm* self = (IsoviscousStressTensorTerm*)matrixTerm;

	_StiffnessMatrixTerm_Delete( self );
}

void _IsoviscousStressTensorTerm_Print( void* matrixTerm, Stream* stream ) {
	IsoviscousStressTensorTerm* self = (IsoviscousStressTensorTerm*)matrixTerm;
	
	_StiffnessMatrixTerm_Print( self, stream );

	/* General info */
	Journal_PrintValue( stream, self->viscosity );
}

void* _IsoviscousStressTensorTerm_DefaultNew( Name name ) {
	return (void*)_IsoviscousStressTensorTerm_New( 
		sizeof(IsoviscousStressTensorTerm), 
		IsoviscousStressTensorTerm_Type,
		_IsoviscousStressTensorTerm_Delete,
		_IsoviscousStressTensorTerm_Print,
		NULL,
		_IsoviscousStressTensorTerm_DefaultNew,
		_IsoviscousStressTensorTerm_Construct,
		_IsoviscousStressTensorTerm_Build,
		_IsoviscousStressTensorTerm_Initialise,
		_IsoviscousStressTensorTerm_Execute,
		_IsoviscousStressTensorTerm_Destroy,
		_IsoviscousStressTensorTerm_AssembleElement,
		name );
}

void _IsoviscousStressTensorTerm_Construct( void* matrixTerm, Stg_ComponentFactory* cf, void* data ) {
	IsoviscousStressTensorTerm*            self             = (IsoviscousStressTensorTerm*)matrixTerm;

	/* Construct Parent */
	_StiffnessMatrixTerm_Construct( self, cf, data );

	_IsoviscousStressTensorTerm_Init( self, Stg_ComponentFactory_GetDouble( cf, self->name, "viscosity", 1.0 ) );
}

void _IsoviscousStressTensorTerm_Build( void* matrixTerm, void* data ) {
	IsoviscousStressTensorTerm*             self             = (IsoviscousStressTensorTerm*)matrixTerm;

	_StiffnessMatrixTerm_Build( self, data );
}

void _IsoviscousStressTensorTerm_Initialise( void* matrixTerm, void* data ) {
	IsoviscousStressTensorTerm*             self             = (IsoviscousStressTensorTerm*)matrixTerm;

	_StiffnessMatrixTerm_Initialise( self, data );
}

void _IsoviscousStressTensorTerm_Execute( void* matrixTerm, void* data ) {
	_StiffnessMatrixTerm_Execute( matrixTerm, data );
}

void _IsoviscousStressTensorTerm_Destroy( void* matrixTerm, void* data ) {
	_StiffnessMatrixTerm_Destroy( matrixTerm, data );
}


void _IsoviscousStressTensorTerm_AssembleElement( 
		void*                                              matrixTerm,
		StiffnessMatrix*                                   stiffnessMatrix, 
		Element_LocalIndex                                 lElement_I, 
		SystemLinearEquations*                             sle,
		FiniteElementContext*                              context,
		double**                                           elStiffMat ) 
{
	IsoviscousStressTensorTerm*         self         = Stg_CheckType( matrixTerm, IsoviscousStressTensorTerm );
	Swarm*                              swarm        = self->integrationSwarm;
	FeVariable*                         variable1    = stiffnessMatrix->rowVariable;
	Dimension_Index                     dim          = stiffnessMatrix->dim;
	
	IntegrationPoint*                   currIntegrationPoint;
	double*                             xi;
	double                              weight;
	Particle_InCellIndex                cParticle_I, cellParticleCount;
	Index                               nodesPerEl;
	Node_ElementLocalIndex              rowNode_I;
	Node_ElementLocalIndex              colNode_I;
	double**                            GNx;
	double                              detJac;
	double                              visc           = self->viscosity;
	
	Cell_Index                          cell_I;
	ElementType*                        elementType;
	int                                 nodalDofs;
	Dof_Index                           diagDof_I;
	Dof_Index                           rowDof_I;
	Dof_Index                           colDof_I;
	Dof_Index                           d;
	Dof_Index                           dofsPerNode;
	Index                               rowIndex;
	Index                               colIndex;
	
	/* Set the element type */
	elementType = FeMesh_ElementTypeAt( variable1->feMesh, lElement_I );
	nodesPerEl = elementType->nodeCount;

	/* assumes constant number of dofs per element */
	nodalDofs = nodesPerEl * dim;
	dofsPerNode = dim;
	
	/* allocate */
	GNx = Memory_Alloc_2DArray( double, dim, nodesPerEl, "GNx" );

	cell_I = CellLayout_MapElementIdToCellId( swarm->cellLayout, lElement_I );
	cellParticleCount = swarm->cellParticleCountTbl[ cell_I ];
	
	for( cParticle_I=0; cParticle_I < cellParticleCount; cParticle_I++ ) {
		
		currIntegrationPoint = (IntegrationPoint*)Swarm_ParticleInCellAt( swarm, cell_I, cParticle_I );
		xi = currIntegrationPoint->xi;
		weight = currIntegrationPoint->weight;
		
		ElementType_ShapeFunctionsGlobalDerivs( 
			elementType,
			variable1->feMesh, lElement_I,
			xi, dim, &detJac, GNx );
		
		/* Initially just build the diagonal (from a Dof perspective) entries */
		for( rowNode_I=0; rowNode_I < nodesPerEl; rowNode_I++ ) { 		
			for( colNode_I=0; colNode_I < nodesPerEl; colNode_I++ ) {		
				for( diagDof_I=0; diagDof_I < dofsPerNode; diagDof_I++ ) {	
					rowIndex =  rowNode_I*dofsPerNode + diagDof_I;
					colIndex = colNode_I*dofsPerNode + diagDof_I;
					elStiffMat[rowIndex][colIndex] +=
						detJac * visc * weight * ( GNx[diagDof_I][rowNode_I] * GNx[diagDof_I][colNode_I] );
					for( d=0; d<dofsPerNode; d++ ) {	
						elStiffMat[rowIndex][colIndex] +=
							detJac * visc * weight * ( GNx[d][rowNode_I] * GNx[d][colNode_I] );
					}
				}
			}
		}
		
		/* then build one set of the off diagonal blocks of K, and copy to the other off diags using symmetry */
		for( rowNode_I=0; rowNode_I < nodesPerEl; rowNode_I++ ) {
			for( rowDof_I=0; rowDof_I < dofsPerNode-1; rowDof_I++ ) {	// dont need to do anything in last row
				for( colNode_I=0; colNode_I < nodesPerEl; colNode_I++ ) {
					for( colDof_I=rowDof_I+1; colDof_I < dofsPerNode; colDof_I++ ) {
						rowIndex = rowNode_I*dofsPerNode + rowDof_I;
						colIndex = colNode_I*dofsPerNode + colDof_I;
						elStiffMat[rowIndex][colIndex] += detJac * visc * weight * 
								( GNx[colDof_I][rowNode_I] * GNx[rowDof_I][colNode_I] );
						/* symmetry	*/
						elStiffMat[colIndex][rowIndex] = elStiffMat[rowIndex][colIndex];
					}
				}
			}
		}
	}
	
	/* free */
	Memory_Free(GNx); 
	
	return;
}
