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
** $Id:  $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/Discretisation/Discretisation.h>
#include <StgFEM/SLE/LinearAlgebra/LinearAlgebra.h>

#include "types.h"
#include "StiffRemesher.h"


/* Textual name of this class */
const Type StiffRemesher_Type = "StiffRemesher";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

#define REMESHER_DEFARGS				\
	sizeof(StiffRemesher),				\
	StiffRemesher_Type,				\
	_StiffRemesher_Delete,				\
	_StiffRemesher_Print,				\
	NULL,						\
	(void*(*)(Name))_StiffRemesher_DefaultNew,	\
	_StiffRemesher_Construct,			\
	_StiffRemesher_Build,				\
	_StiffRemesher_Initialise,			\
	_StiffRemesher_Execute,				\
	_StiffRemesher_Destroy,				\
	name,						\
	False,						\
	_StiffRemesher_SetMesh


StiffRemesher* StiffRemesher_New( Name name ) {
	return _StiffRemesher_New( REMESHER_DEFARGS );
}


StiffRemesher* _StiffRemesher_New( CLASS_ARGS, 
				   COMPONENT_ARGS, 
				   REMESHER_ARGS )
{
	StiffRemesher*	self;

	/* Allocate memory. */
	self = (StiffRemesher*)_Remesher_New( _sizeOfSelf,
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
					    name, 
					      initFlag, 
					      setMeshFunc );

	/* StiffRemesher info */
	if( initFlag ) {
		_StiffRemesher_Init( self );
	}

	return self;
}


void StiffRemesher_Init( StiffRemesher* self ) {
	assert( 0 ); /* TODO */
#if 0
	/* General info */
	self->type = StiffRemesher_Type;
	self->_sizeOfSelf = sizeof(StiffRemesher);
	self->_deleteSelf = False;
	
	/* Virtual info */
	self->_delete = _StiffRemesher_Delete;
	self->_print = _StiffRemesher_Print;
	self->_copy = NULL;
	_Stg_Class_Init( (Stg_Class*)self );
	
	/* StiffRemesher info */
	_StiffRemesher_Init( self );
#endif
}


void _StiffRemesher_Init( StiffRemesher* self ) {
	/* StiffRemesher info */
	memset( &self->nDims, 
		0, 
		(size_t)&self->matSolver - (size_t)&self->rests + sizeof(MatrixSolver*) );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _StiffRemesher_Delete( void* stiffRemesher ) {
	StiffRemesher*	self = (StiffRemesher*)stiffRemesher;

	/* Delete the class itself */
	_StiffRemesher_Free( self );

	/* Delete parent */
	_Stg_Component_Delete( stiffRemesher );
}


void _StiffRemesher_Print( void* stiffRemesher, Stream* stream ) {
	StiffRemesher*	self = (StiffRemesher*)stiffRemesher;
	Stream*		myStream;
	
	/* Set the Journal for printing informations */
	myStream = Journal_Register( InfoStream_Type, "StiffRemesherStream" );

	/* Print parent */
	_Stg_Component_Print( self, stream );

	/* General info */
	Journal_Printf( myStream, "StiffRemesher (ptr): (%p)\n", self );

	/* Virtual info */

	/* StiffRemesher info */
}


StiffRemesher* _StiffRemesher_DefaultNew( Name name ) {
	return _StiffRemesher_New( REMESHER_DEFARGS );
}


void _StiffRemesher_Construct( void* stiffRemesher, Stg_ComponentFactory* cf, void* data ) {
	StiffRemesher*	self = (StiffRemesher*)stiffRemesher;
	Dictionary*	dict;
	char*		swarmName;

	assert( self );
	assert( cf );
	assert( cf->componentDict );

	_Remesher_Construct( self, cf, data );

	/* Get the number of dimensions. */
	dict = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( cf->componentDict, self->name ) );
	self->nDims = Dictionary_GetUnsignedInt( dict, "dims" );
	assert( self->nDims > 0 );

	/* Get the swarm. */
	swarmName = Dictionary_GetString( dict, "swarm" );
	assert( swarmName );
	self->swarm = Stg_ComponentFactory_ConstructByName( cf, swarmName, Swarm, True, data ); 
}


void _StiffRemesher_Build( void* stiffRemesher, void* data ) {
	StiffRemesher*	self = (StiffRemesher*)stiffRemesher;

	assert( self );

	if( !self->meshType ) {
		return;
	}

	assert( self->mesh );

	/* Build parent. */
	_Remesher_Build( self, data );

	/* Set the restricted nodes based on mesh type. */
	if( !strcmp( self->meshType, "regular" ) ) {
		GRM		grm;
		unsigned*	dimInds;
		unsigned	n_i;

		/* When regular meshing is required, boundary nodes are restricted. Construct a boundary set of indices. */
		RegMesh_Generalise( self->mesh, &grm );
		dimInds = Memory_Alloc_Array_Unnamed( unsigned, grm.nDims );
		for( n_i = 0; n_i < self->mesh->nodeLocalCount; n_i++ ) {
			unsigned	gNodeInd = Mesh_NodeMapLocalToGlobal( self->mesh, n_i );
			unsigned	d_i;

			GRM_Lift( &grm, gNodeInd, dimInds );
			for( d_i = 0; d_i < grm.nDims; d_i++ ) {
				if( dimInds[d_i] == 0 || dimInds[d_i] == grm.nNodes[d_i] - 1 ) {
					StiffRemesher_SetRestriction( self, n_i, d_i, True );

					/* If on the top surface... */
					if( d_i == 1 && dimInds[d_i] == grm.nNodes[d_i] - 1 ) {
						unsigned	d_j;

						/* Restrict in every dimension. */
						for( d_j = 0; d_j < grm.nDims; d_j++ ) {
							StiffRemesher_SetRestriction( self, n_i, d_j, True );
						}
					}
				}
			}
		}

		/* Free dimensional indices array. */
		FreeArray( dimInds );
	}
	else {
		assert( 0 );
	}
}


void _StiffRemesher_Initialise( void* stiffRemesher, void* data ) {
	StiffRemesher*	self = (StiffRemesher*)stiffRemesher;

	assert( self );

	/* Initialise parent. */
	_Remesher_Initialise( self, data );

	/* Initialise the system. */
	StiffRemesher_BuildSystem( self );
}


void _StiffRemesher_Execute( void* stiffRemesher, void* data ) {
	StiffRemesher*	self = (StiffRemesher*)stiffRemesher;
	unsigned	nLNodes;
	double*		solArray;
	unsigned	n_i, d_i;

	assert( self );
	assert( self->mesh );
	assert( self->mesh->layout );
	assert( self->mesh->layout->elementLayout );
	/* TODO: remaining asserts */

	nLNodes = self->mesh->nodeLocalCount;

	/* Evaluate the weights. */
	_StiffRemesher_UpdateWeights( self );

	/* Fill RHS. */
	Matrix_Zero( self->stiffMat );
	Vector_Zero( self->rhsVec );
	Vector_Zero( self->solVec );
	for( n_i = 0; n_i < nLNodes; n_i++ ) {
		for( d_i = 0; d_i < self->nDims; d_i++ ) {
			unsigned	eqNumI;
			unsigned	nNbrs;
			unsigned*	nbrs;
			unsigned	nbr_i;

			/* If restricted in this dimension, skip. */
			if( self->rests[n_i][d_i] ) {
				continue;
			}

			/* Collect neighbour info. */
			nNbrs = self->mesh->nodeNeighbourCountTbl[n_i];
			nbrs = self->mesh->nodeNeighbourTbl[n_i];

			/* Calc equation number and coef. */
			eqNumI = self->baseEqNum + self->eqNums[n_i][d_i];

			/* Enter matrix values. */
			Matrix_AddValue( self->stiffMat, eqNumI, eqNumI, 1.0 );
			for( nbr_i = 0; nbr_i < nNbrs; nbr_i++ ) {
				/* If this neighbour is not restricted, add to matrix. */
				if( !self->rests[nbrs[nbr_i]][d_i] ) {
					unsigned	eqNumJ;

					eqNumJ = self->baseEqNum + self->eqNums[nbrs[nbr_i]][d_i];
					Matrix_AddValue( self->stiffMat, eqNumI, eqNumJ, -self->weights[n_i][nbr_i][d_i] );
				}
			}

			/* Enter RHS value. */
			for( nbr_i = 0; nbr_i < nNbrs; nbr_i++ ) {
				/* If this neighbour is restricted, add to RHS. */
				if( self->rests[nbrs[nbr_i]][d_i] ) {
					Vector_AddEntry( self->rhsVec, eqNumI, 
							 self->weights[n_i][nbr_i][d_i] * self->mesh->nodeCoord[nbrs[nbr_i]][d_i] );
				}
			}
		}
	}

	/* Assemble the system. */
	Matrix_AssemblyBegin( self->stiffMat );
	Vector_AssemblyBegin( self->solVec );
	Vector_AssemblyBegin( self->rhsVec );
	Matrix_AssemblyEnd( self->stiffMat );
	Vector_AssemblyEnd( self->solVec );
	Vector_AssemblyEnd( self->rhsVec );

	/* Solve the system. */
	MatrixSolver_Solve( self->matSolver, self->solVec, self->rhsVec );

	/* Update mesh coordinates. */
	Vector_Get( self->solVec, &solArray );
	for( n_i = 0; n_i < nLNodes; n_i++ ) {
		for( d_i = 0; d_i < self->nDims; d_i++ ) {
			unsigned	eqNum;

			if( self->rests[n_i][d_i] ) {
				continue;
			}

			/* Copy the coordinate from the vector. */
			eqNum = self->eqNums[n_i][d_i] - self->baseEqNum;
			self->mesh->nodeCoord[n_i][d_i] = solArray[eqNum];
		}
	}
	Vector_Restore( self->solVec, &solArray );
}


void _StiffRemesher_Destroy( void* stiffRemesher, void* data ) {
	StiffRemesher*	self = (StiffRemesher*)stiffRemesher;

	assert( self );

	/* TODO: If delete deletes, what does destroy do? */
}


void _StiffRemesher_SetMesh( void* stiffRemesher, Mesh* mesh ) {
	StiffRemesher*	self = (StiffRemesher*)stiffRemesher;

	assert( self );
	assert( self->mesh->layout );
	assert( self->mesh->layout->decomp );

	/* Kill all internals. */
	_StiffRemesher_Free( self );

	/* Store the mesh and communicator. */
	self->mesh = mesh;
	self->comm = mesh->layout->decomp->communicator;

	/* Allocate for element volume approximations. */
	_StiffRemesher_CalcWeights( self );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void StiffRemesher_SetRestriction( void* stiffRemesher, unsigned lNodeInd, unsigned dim, Bool state ) {
	StiffRemesher*	self = (StiffRemesher*)stiffRemesher;

	assert( self );
	assert( self->mesh );
	assert( lNodeInd < self->mesh->nodeLocalCount );
	assert( self->mesh->layout );
	assert( self->mesh->layout->elementLayout );
	assert( dim < ((HexaEL*)self->mesh->layout->elementLayout)->dim );

	/* If we this is the first restriction, allocate the arrays. */
	if( !self->rests ) {
		unsigned	nLNodes;
		unsigned	n_i;

		nLNodes = self->mesh->nodeLocalCount;

		self->rests = Memory_Alloc_2DArray( unsigned, nLNodes, self->nDims, "StiffRemesher->rests" );
		for( n_i = 0; n_i < nLNodes; n_i++ ) {
			memset( self->rests[n_i], 0, sizeof(Bool) * self->nDims );
		}
	}

	/* If the new state is different from existing state, destroy the current system. */
	if( self->rests[lNodeInd][dim] != state ) {
		_StiffRemesher_FreeSystem( self );
	}

	/* Set the restriction. */
	self->rests[lNodeInd][dim] = state;
}


void StiffRemesher_BuildSystem( void* stiffRemesher ) {
	const unsigned	eqNumTag = 1010;
	StiffRemesher*	self = (StiffRemesher*)stiffRemesher;
	unsigned	nRows = 0;
	unsigned	nonZeros = 0;
	unsigned	nLNodes;
	unsigned	rank, nProcs;
	unsigned	n_i, d_i;

	assert( self );
	assert( self->mesh );
	assert( self->mesh->nodeNeighbourCountTbl );
	assert( self->mesh->nodeNeighbourTbl );
	assert( self->mesh->layout );
	assert( self->mesh->layout->elementLayout );
	/* TODO: remaining assertions */

	MPI_Comm_rank( self->comm, (int*)&rank );
	MPI_Comm_size( self->comm, (int*)&nProcs );

	nLNodes = self->mesh->nodeLocalCount;

	/* Determine the number of rows and build local equation number mappings. */
	self->eqNums = Memory_Alloc_2DArray( unsigned, nLNodes, self->nDims, "StiffRemesher->eqNums" );
	for( n_i = 0; n_i < nLNodes; n_i++ ) {
		for( d_i = 0; d_i < self->nDims; d_i++ ) {
			if( !self->rests[n_i][d_i] ) {
				self->eqNums[n_i][d_i] = nRows++;
			}
			else {
				self->eqNums[n_i][d_i] = (unsigned)-1;
			}
		}
	}

	/* Determine the average number of non-zeros per row. */
	for( n_i = 0; n_i < nLNodes; n_i++ ) {
		for( d_i = 0; d_i < self->nDims; d_i++ ) {
			/* Skip if restricted in this dimension. */
			if( !self->rests[n_i][d_i] ) {
				/* Set the non-zero count to be the number of neighbours. */
				nonZeros += self->mesh->nodeNeighbourCountTbl[n_i];
			}
		}
	}
	nonZeros /= nRows;

	/* Send and receive base equation numbers. */
	if( rank > 0 ) {
		MPI_Status	status;

		MPI_Recv( &self->baseEqNum, 1, MPI_UNSIGNED, rank - 1, eqNumTag, self->comm, &status );
	}
	else {
		self->baseEqNum = 0;
	}

	if( rank < nProcs - 1 ) {
		self->baseEqNum += nRows;
		MPI_Send( &self->baseEqNum, 1, MPI_UNSIGNED, rank + 1, eqNumTag, self->comm );
		self->baseEqNum -= nRows;
	}

	/* Create the system. */
	self->stiffMat = Matrix_New( self->comm, nRows, nRows, nonZeros );
	self->solVec = Vector_New_SpecifyLocalSize( self->comm, nRows );
	self->rhsVec = Vector_New_SpecifyLocalSize( self->comm, nRows );
	Matrix_Zero( self->stiffMat );
	Vector_Zero( self->solVec );
	Vector_Zero( self->rhsVec );

	/* Create the matrix solver. */
	self->matSolver = MatrixSolver_Build( self->comm, self->stiffMat );
	MatrixSolver_Setup( self->matSolver, self->stiffMat );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void _StiffRemesher_Free( StiffRemesher* self ) {
	assert( self );

	KillArray( self->weights ); /* FIX */
	KillArray( self->rests );
	KillArray( self->eqNums );
	_StiffRemesher_FreeSystem( self );
}


void _StiffRemesher_FreeSystem( StiffRemesher* self ) {
	assert( self );

	if(  self->stiffMat ) {
		Matrix_Destroy( self->stiffMat );
		self->stiffMat = NULL;
	}

	if( self->solVec ) {
		Vector_Destroy( self->solVec );
		self->solVec = NULL;
	}

	if( self->rhsVec ) {
		Vector_Destroy( self->rhsVec );
		self->rhsVec = NULL;
	}

	if( self->matSolver ) {
		MatrixSolver_Destroy( self->matSolver );
		self->matSolver = NULL;
	}
}


void _StiffRemesher_UpdateWeights( StiffRemesher* self ) {
	unsigned	nNodes;
	Coord*		nodeCrds;
	double*		dispSums;
	unsigned	n_i;

	assert( self );

	nodeCrds = self->mesh->nodeCoord;
	nNodes = self->mesh->nodeLocalCount;
	dispSums = Memory_Alloc_Array_Unnamed( double, self->nDims );
	for( n_i = 0; n_i < nNodes; n_i++ ) {
		unsigned	nNbrs;
		unsigned*	nbrs;
		double*		crd = nodeCrds[n_i];
		double**	disps;
		unsigned	nbr_i, d_i;

		/* Get neighbours. */
		nNbrs = self->mesh->nodeNeighbourCountTbl[n_i];
		nbrs = self->mesh->nodeNeighbourTbl[n_i];

		/* Allocate some stuff. */
		disps = Memory_Alloc_2DArray_Unnamed( double, nNbrs, self->nDims );

		/* Calculate unit vectors'n'such. */
		memset( dispSums, 0, self->nDims * sizeof(double) );
		for( nbr_i = 0; nbr_i < nNbrs; nbr_i++ ) {
			double	mag = 0.0;

			for( d_i = 0; d_i < self->nDims; d_i++ ) {
				/* Calculate dimensional influence. */
				disps[nbr_i][d_i] = nodeCrds[nbrs[nbr_i]][d_i] - crd[d_i];
				disps[nbr_i][d_i] *= (disps[nbr_i][d_i] < 0.0) ? -1.0 : 1.0;
				mag += disps[nbr_i][d_i] * disps[nbr_i][d_i];
			}
			mag = 1.0 / sqrt( mag );
			for( d_i = 0; d_i < self->nDims; d_i++ ) {
				disps[nbr_i][d_i] *= mag;
				dispSums[d_i] += disps[nbr_i][d_i];
			}
		}

		/* Build the weights. */
		for( d_i = 0; d_i < self->nDims; d_i++ ) {
			double	inv = 1.0 / dispSums[d_i];

			for( nbr_i = 0; nbr_i < nNbrs; nbr_i++ ) {
				self->weights[n_i][nbr_i][d_i] = disps[nbr_i][d_i] * inv;
			}
		}

		FreeArray( disps );
	}

	FreeArray( dispSums );
}


void _StiffRemesher_CalcWeights( StiffRemesher* self ) {
	unsigned	n_i;

	assert( self );

	self->weights = Memory_Alloc_Array_Unnamed( double**, self->mesh->nodeLocalCount );
	for( n_i = 0; n_i < self->mesh->nodeLocalCount; n_i++ ) {
		unsigned	nNbrs;

		nNbrs = self->mesh->nodeNeighbourCountTbl[n_i];
		self->weights[n_i] = Memory_Alloc_2DArray_Unnamed( double, nNbrs, self->nDims );
	}

	_StiffRemesher_UpdateWeights( self );
}


#if 0
void _StiffRemesher_CalcWeights( StiffRemesher* self ) {
	unsigned	nEls, nNodes;
	double*		volumes;
	double		avgVol = 0.0;
	unsigned	e_i, n_i;

	assert( self );
	assert( self->weights );
	assert( self->mesh->type == FiniteElement_Mesh_Type );

	nEls = self->mesh->elementLocalCount;
	nNodes = self->mesh->nodeLocalCount;

	/* Allocate for new volumes. */
	volumes = Memory_Alloc_Array_Unnamed( double, nEls );

	/* Calculate current volumes. */
	for( e_i = 0; e_i < nEls; e_i++ ) {
		ElementType*	elType;
		unsigned	nElNodes;
		unsigned	cellInd;
		unsigned	nParticles;
		double		vol = 0.0;
		double**	gDerivs;
		unsigned	p_i;

		/* Get the element type. */
		elType = FeMesh_ElementTypeAt( (FiniteElement_Mesh*)self->mesh, e_i );
		nElNodes = elType->nodeCount;

		/* Allocate for global derivatives. */
		gDerivs = Memory_Alloc_2DArray_Unnamed( double, self->nDims, nElNodes );

		/* Get particle info. */
		cellInd = CellLayout_MapElementIdToCellId( self->swarm->cellLayout, e_i );
		nParticles = self->swarm->cellParticleCountTbl[cellInd];

		/* Integrate domain. */
		for( p_i = 0; p_i < nParticles; p_i++ ) {
			IntegrationPoint*	intPoint;
			double			jacDet;

			/* Get the current integration point. */
			intPoint = (IntegrationPoint*)Swarm_ParticleInCellAt( self->swarm, cellInd, p_i );

			/* Evaluate the jacobian. */
			ElementType_ShapeFunctionsGlobalDerivs( elType, self->mesh, e_i, intPoint->xi, self->nDims, 
								&jacDet, gDerivs );

			/* Add to total. */
			vol += intPoint->weight * jacDet;
		}

		/* Free global derivs. */
		FreeArray( gDerivs );

		/* Store the volume. */
		volumes[e_i] = vol;
		printf( "volume[%d]: \t%g\n", e_i, vol );

		/* Calc average. */
		avgVol += vol;
	}

	/* Finish average. */
	avgVol /= (double)nEls;

	/* Calculate nodal weights. */
	avgVol = 1.0 / avgVol;
	for( n_i = 0; n_i < nNodes; n_i++ ) {
		unsigned	nIncEls = self->mesh->nodeElementCountTbl[n_i];
		unsigned*	incEls = self->mesh->nodeElementTbl[n_i];
		unsigned	nRealIncEls = 0;
		unsigned	ne_i;

		self->weights[n_i] = 0.0;
		for( ne_i = 0; ne_i < nIncEls; ne_i++ ) {
			if( incEls[ne_i] < nEls ) {
				self->weights[n_i] += volumes[incEls[ne_i]] * avgVol;
				nRealIncEls++;
			}
		}
		self->weights[n_i] /= (double)nRealIncEls;
	}

	/* Free volumes array. */
	FreeArray( volumes );
}
#endif
