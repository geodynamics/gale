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

#include "StGermain/StGermain.h"
#include "StgFEM/Discretisation/Discretisation.h"
#include "StgFEM/SLE/SLE.h"
#include "StgFEM/Assembly/Assembly.h"

#include "types.h"
#include "Polynomial.h"
#include "Quadrature.h"
#include "Context.h"
#include "LevelSet.h"


/* Textual name of this class */
const Type LevelSet_Type = "LevelSet";


unsigned short extractBit( unsigned n, unsigned b ) {
	return ((n >> b) & 1) ? 1 : 0;
}




/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

LevelSet* LevelSet_DefaultNew( Name name ) {
	return _LevelSet_New( sizeof(LevelSet), 
			      LevelSet_Type, 
			      _LevelSet_Delete, 
			      _LevelSet_Print, 
			      NULL, 
			      (void*(*)(Name))LevelSet_DefaultNew,
			      _LevelSet_Construct,
			      _LevelSet_Build, 
			      _LevelSet_Initialise, 
			      _LevelSet_Execute, 
			      _LevelSet_Destroy, 
			      name,
			      False,
			      _LevelSet_InterpValAt, 
			      _LevelSet_GetMinField, 
			      _LevelSet_GetMaxField, 
			      _LevelSet_GetLocalCrds, 
			      _LevelSet_GetGlobalCrds, 
			      1,
			      0,
				  False,
			      MPI_COMM_WORLD,
			      NULL );
}


LevelSet* LevelSet_New( Name			name, 
			FieldVariable_Register*	fieldVariable_Register )
{
	return _LevelSet_New( sizeof(LevelSet), 
			      LevelSet_Type, 
			      _LevelSet_Delete, 
			      _LevelSet_Print, 
			      NULL, 
			      (void*(*)(Name))FieldVariable_DefaultNew,
			      _LevelSet_Construct,
			      _LevelSet_Build, 
			      _LevelSet_Initialise, 
			      _LevelSet_Execute, 
			      _LevelSet_Destroy,
			      name,
			      True,
			      _LevelSet_InterpValAt, 
			      _LevelSet_GetMinField, 
			      _LevelSet_GetMaxField, 
			      _LevelSet_GetLocalCrds, 
			      _LevelSet_GetGlobalCrds, 
			      1, 
			      0,
				  False,	
			      MPI_COMM_WORLD, 
			      fieldVariable_Register );
}


LevelSet* _LevelSet_New( SizeT						_sizeOfSelf, 
			 Type						type,
			 Stg_Class_DeleteFunction*			_delete,
			 Stg_Class_PrintFunction*			_print, 
			 Stg_Class_CopyFunction*			_copy, 
			 Stg_Component_DefaultConstructorFunction*	_defaultConstructor, 
			 Stg_Component_ConstructFunction*		_construct,
			 Stg_Component_BuildFunction*			_build,
			 Stg_Component_InitialiseFunction*		_initialise,
			 Stg_Component_ExecuteFunction*			_execute,
			 Stg_Component_DestroyFunction*			_destroy,
			 Name						name,
			 Bool						initFlag,
			 FieldVariable_InterpolateValueAtFunction*	_interpolateValueAt,
			 FieldVariable_GetValueFunction*		_getMinGlobalFieldMagnitude,
			 FieldVariable_GetValueFunction*		_getMaxGlobalFieldMagnitude,		
			 FieldVariable_GetCoordFunction*		_getMinAndMaxLocalCoords,
			 FieldVariable_GetCoordFunction*		_getMinAndMaxGlobalCoords,
			 Index						            fieldComponentCount,
			 Dimension_Index           				dim,
			 Bool                                   isCheckpointedAndReloaded,
			 MPI_Comm					           communicator,
			 FieldVariable_Register*			fieldVariable_Register )
{
	LevelSet*	self;

	/* Allocate memory. */
	self = (LevelSet*)_FieldVariable_New( _sizeOfSelf,
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
					      _interpolateValueAt, 
					      _getMinGlobalFieldMagnitude, 
					      _getMaxGlobalFieldMagnitude, 
					      _getMinAndMaxLocalCoords, 
					      _getMinAndMaxGlobalCoords, 
					      fieldComponentCount, 
					      dim, 
						  isCheckpointedAndReloaded,
					      communicator, 
					      fieldVariable_Register );
	
	/* General info */
	
	/* Virtual info */
	
	/* LevelSet info */
	if( initFlag ) {
		_LevelSet_Init( self );
	}
	
	return self;
}


void LevelSet_Init( LevelSet*			self, 
		    FieldVariable_Register*	fV_Register )
{
	assert( 0 ); /* TODO */
#if 0
	/* General info */
	self->type = LevelSet_Type;
	self->_sizeOfSelf = sizeof(LevelSet);
	self->_deleteSelf = False;
	
	/* Virtual info */
	self->_delete = _LevelSet_Delete;
	self->_print = _LevelSet_Print;
	self->_copy = NULL;
	_Stg_Class_Init( (Stg_Class*)self );
	
	/* LevelSet info */
	_LevelSet_Init( self );
#endif
}


void _LevelSet_Init( LevelSet* self ) {
	/* General and Virtual info should already be set */
	
	/* LevelSet info */
	memset( &self->velVarName, 0, (size_t)&self->solver - (size_t)&self->grm + sizeof(MatrixSolver*) );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _LevelSet_Delete( void* levelSet ) {
	LevelSet*	self = (LevelSet*)levelSet;
	
	/* Delete the class itself */
	_LevelSet_Free( self );
	
	/* Delete parent */
	_FieldVariable_Delete( levelSet );
}


void _LevelSet_Print( void* levelSet, Stream* stream ) {
	LevelSet*	self = (LevelSet*)levelSet;
	Stream*		myStream;
	
	/* Set the Journal for printing informations */
	myStream = Journal_Register( InfoStream_Type, "LevelSetStream" );
	
	/* Print parent */
	_Stg_Class_Print( self, stream );
	
	/* General info */
	Journal_Printf( myStream, "LevelSet (ptr): (%p)\n", self );
	
	/* Virtual info */
	
	/* LevelSet info */
	Journal_Printf( myStream, "\tJacobian determinant: %g\n", self->jacDet );

	printf( "\n\n\n" );
	Matrix_View( self->aMat, stream );
	printf( "\n\n\n" );
	Vector_View( self->rhsVec, stream );
	printf( "\n\n\n" );
	Vector_View( self->solVec, stream );
}


void _LevelSet_Construct( void* levelSet, Stg_ComponentFactory* cf, void* data ) {
	LevelSet*		self = (LevelSet*)levelSet;
	Dictionary*		dict;

	assert( self );
	assert( cf );
	assert( cf->componentDict );

	_FieldVariable_Construct( levelSet, cf, data );
	_LevelSet_Init( levelSet );
	dict = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( cf->componentDict, self->name ) );
	assert( dict );

	self->velVarName = StG_Strdup( Dictionary_Entry_Value_AsString( Dictionary_Get( dict, "velocityVariable" ) ) );
	self->res = Dictionary_Entry_Value_AsInt( Dictionary_Get( dict, "resolution" ) );
	self->fac = Dictionary_Entry_Value_AsInt( Dictionary_Get( dict, "scaleFactor" ) );
	self->initDistDict = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( dict, "initialDistance" ) );
	self->cutoff = Dictionary_Entry_Value_AsDouble( Dictionary_Get( dict, "cutoff" ) );
}


void _LevelSet_Build( void* levelSet, void* data ) {
	_FieldVariable_Build( levelSet, data );
}


void _LevelSet_Execute( void* levelSet, void* data ) {
	LevelSet*		self = (LevelSet*)levelSet;
	FiniteElementContext*	feCtx = (FiniteElementContext*)data;

	assert( self );

	_FieldVariable_Execute( levelSet, data );
	LevelSet_Advance( levelSet, feCtx->dt );
}


void _LevelSet_Destroy( void* levelSet, void* data ) {
	_LevelSet_Free( levelSet );
	_FieldVariable_Destroy( levelSet, data );
}


void _LevelSet_Initialise( void* levelSet, void* data ) {
	LevelSet*		self = (LevelSet*)levelSet;
	LevelSetPlugin_Context*	lsCtx;
	FiniteElementContext*	feCtx;
	FieldVariable*		velVar;

	_FieldVariable_Initialise( levelSet, data );

	/* Store the context. */
	feCtx = (FiniteElementContext*)data;

	/* Locate velocity variable and initialise. */
	velVar = FieldVariable_Register_GetByName( self->fieldVariable_Register, self->velVarName );
	assert( velVar );
	LevelSet_SetVelocityField( self, velVar );

	/* Fill distance values. */
	if( self->initDistDict ) {
		char*	distType;

		distType = Dictionary_Entry_Value_AsString( Dictionary_Get( self->initDistDict, "type" ) );
		if( !strcmp( distType, "circular" ) ) {
			unsigned	nDims;
			Coord		origin;
			double		radius;

			nDims = Dictionary_Entry_Value_AsInt( Dictionary_Get( self->initDistDict, "dimensions" ) );
			assert( nDims > 0 && nDims <= 3 );
			origin[0] = Dictionary_Entry_Value_AsDouble( 
				Dictionary_GetDefault( self->initDistDict, "originX", Dictionary_Entry_Value_FromDouble( 0 ) ) );
			origin[1] = Dictionary_Entry_Value_AsDouble( 
				Dictionary_GetDefault( self->initDistDict, "originY", Dictionary_Entry_Value_FromDouble( 0 ) ) );
			origin[2] = Dictionary_Entry_Value_AsDouble( 
				Dictionary_GetDefault( self->initDistDict, "originZ", Dictionary_Entry_Value_FromDouble( 0 ) ) );
			radius = Dictionary_Entry_Value_AsDouble( 
				Dictionary_GetDefault( self->initDistDict, "radius", Dictionary_Entry_Value_FromDouble( 1 ) ) );
			assert( radius > 0.0 );

			/* Fill the distances. */
			if( nDims == 1 ) {
				assert( 0 );
			}
			else if( nDims == 2 ) {
				unsigned	lNodeInd;
				double		dist;
				Coord		crd;
				unsigned	inds[2];

				assert( self->grm.nDims == 2 );
				for( inds[1] = 0; inds[1] < self->grm.nNodes[1]; inds[1]++ ) {
					crd[1] = self->min[1] + (double)inds[1] * self->step[1] - origin[1];
					for( inds[0] = 0; inds[0] < self->grm.nNodes[0]; inds[0]++ ) {
						crd[0] = self->min[0] + (double)inds[0] * self->step[0] - origin[0];
						dist = sqrt( crd[0] * crd[0] + crd[1] * crd[1] ) - radius;
						GRM_Project( &self->grm, inds, &lNodeInd );
						self->dists[lNodeInd] = dist;
					}
				}
			}
			else {
				assert( 0 );
			}
			
		}
	}

	/* Setup initial active elements. */
	_LevelSet_InitActiveEls( self );

	/* Calc min and max. */
	_LevelSet_CalcMinMax( self );

	/* If the plugin is loaded, add to extension. */
	lsCtx = ExtensionManager_Get( feCtx->extensionMgr,
				      feCtx, 
				      LevelSetPlugin_ContextHandle );
	if( lsCtx ) {
		lsCtx->lSets = Memory_Realloc_Array( lsCtx->lSets, LevelSet*, lsCtx->nLSets + 1 );
		lsCtx->lSets[lsCtx->nLSets] = self;
		lsCtx->nLSets++;
	}
}


InterpolationResult _LevelSet_InterpValAt( void* levelSet, Coord gCrd, double* vecRes ) {
	LevelSet*	self = (LevelSet*)levelSet;
	unsigned	gNode, lNode, lEl;
	unsigned*	inds;
	double		res;
	double*		lCrd;
	unsigned	dim_i, elNode_i;

	assert( self );
	assert( vecRes );


	/*
	** For the moment we will use simple linear interpolation, however as we are using a strictly
	** regular mesh for level sets, we could use a smoother interpolation method if desired.
	**
	** Locate the element in which the global coordinate resides.
	*/

	gNode = 0;
	inds = Memory_Alloc_Array_Unnamed( unsigned, self->grm.nDims );
	for( dim_i = 0; dim_i < self->grm.nDims; dim_i++ ) {
		if( gCrd[dim_i] < self->min[dim_i] || gCrd[dim_i] > self->max[dim_i] ) {
			return OUTSIDE_GLOBAL;
		}
		inds[dim_i] = (unsigned)floor( (gCrd[dim_i] - self->min[dim_i]) / self->step[dim_i] );
		if( inds[dim_i] == self->grm.nNodes[dim_i] ) {
			inds[dim_i]--;
		}
		gNode += inds[dim_i] * self->grm.basis[dim_i];
	}

	if( !Sync_IsLocal( self->nodeSync, gNode ) ) {
		return OTHER_PROC;
	}
	lNode = Sync_MapGlobal( self->nodeSync, gNode );

	/*
	** Calculate local coordinates.
	*/

	lCrd = Memory_Alloc_Array_Unnamed( double, self->grm.nDims );
	for( dim_i = 0; dim_i < self->grm.nDims; dim_i++ ) {
		lCrd[dim_i] = 2.0 * (gCrd[dim_i] - self->min[dim_i] - (double)inds[dim_i] * self->step[dim_i]) / self->step[dim_i] - 1.0;
	}
	FreeArray( inds );


	/*
	** Use the shape functions to calculate the interpolation.
	*/

	res = 0.0;
	lEl = self->baseNodeEl[lNode];
	for( elNode_i = 0; elNode_i < self->nElNodes; elNode_i++ ) {
		/* Add value's contribution. */
		res += Function_Eval( self->shapeFuncs + elNode_i, lCrd ) * self->dists[self->elNodes[lEl][elNode_i]];
	}
	*vecRes = res;
	FreeArray( lCrd );

	return LOCAL;
}


double _LevelSet_GetMinField( void* levelSet ) {
	assert( levelSet );
	return ((LevelSet*)levelSet)->minDist;
}


double _LevelSet_GetMaxField( void* levelSet ) {
	assert( levelSet );
	return ((LevelSet*)levelSet)->maxDist;
}


void _LevelSet_GetLocalCrds( void* levelSet, Coord min, Coord max ) {
	FieldVariable_GetMinAndMaxLocalCoords( ((LevelSet*)levelSet)->velVar, min, max );
}


void _LevelSet_GetGlobalCrds( void* levelSet, Coord min, Coord max ) {
	FieldVariable_GetMinAndMaxGlobalCoords( ((LevelSet*)levelSet)->velVar, min, max );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void LevelSet_SetResolution( void* levelSet, double res, unsigned fac ) {
	LevelSet*	self = (LevelSet*)levelSet;

	assert( self );
	assert( res > 0.0 );
	assert( fac > 0 );

	self->res = res;
	self->fac = fac;
	if( self->velVar ) {
		LevelSet_SetVelocityField( self, self->velVar );
	}
}


void LevelSet_SetVelocityField( void* levelSet, FieldVariable* velVar )
{
	LevelSet*		self = (LevelSet*)levelSet;
	unsigned*		base;
	unsigned*		elNodes;
	Coord			minCrd, maxCrd, pos;
	unsigned		nProcs, rank, len;
	unsigned*		local;
	unsigned		nLocal;
	unsigned*		shadow;
	unsigned		nShadow;
	unsigned		nGNodes;
	unsigned*		gInds;
	InterpolationResult	erpRes;
	unsigned*		remNShad;
	unsigned**		remShad;
	unsigned*		intr;
	unsigned		nIntr;
	unsigned		step;
	double*			dummyVel;
	Function		func1, func2, func3;
	Function*		funcVec;
	double*			crds;
	unsigned		qdr_i;
	unsigned		node_i, node_j, node_k;
	unsigned		sf_i;
	unsigned		gNode_i, lNode_i, intr_i, dim_i;

	void setIntersect( unsigned** sets, unsigned* nEls, 
			   unsigned** intr, unsigned* nIntr, 
			   unsigned nProcs, unsigned rank );

	assert( self );
	assert( velVar );

	_LevelSet_Free( self );


	/*
	** Setup standard values. There is a slight hack in here to set up element mappings.
	*/

	self->communicator = velVar->communicator;
	self->dim = velVar->dim;
	self->grm.nDims = self->dim;
	self->min = Memory_Alloc_Array( double, self->dim, "LevelSet" );
	self->max = Memory_Alloc_Array( double, self->dim, "LevelSet" );
	self->step = Memory_Alloc_Array( double, self->dim, "LevelSet" );
	FieldVariable_GetMinAndMaxGlobalCoords( velVar, minCrd, maxCrd );
	nGNodes = 1;
	for( dim_i = 0; dim_i < self->dim; dim_i++ ) {
		self->min[dim_i] = minCrd[dim_i];
		self->max[dim_i] = maxCrd[dim_i];
		len = (self->max[dim_i] - self->min[dim_i]);
		self->grm.nNodes[dim_i] = self->res + 1;
		if( self->grm.nNodes[dim_i] < 2 ) {
			self->grm.nNodes[dim_i] = 2;
		}
		self->step[dim_i] = len / (double)(self->grm.nNodes[dim_i] - 1);
		self->grm.basis[dim_i] = (dim_i > 0) ? self->grm.basis[dim_i - 1] * self->grm.nNodes[dim_i - 1] : 1;
		nGNodes *= self->grm.nNodes[dim_i];
	}
	self->velVar = velVar;


	/*
	** Set up the sync, will require a bit of messing around.
	*/

	/* Figure out who owns what. */
	local = NULL;	nLocal = 0;
	shadow = NULL;	nShadow = 0;
	gInds = Memory_Alloc_Array_Unnamed( unsigned, self->dim );
	dummyVel = Memory_Alloc_Array_Unnamed( double, self->dim );
	for( gNode_i = 0; gNode_i < nGNodes; gNode_i++ ) {
		/* Calculate position. */
		GRM_Lift( &self->grm, gNode_i, gInds );
		for( dim_i = 0; dim_i < self->dim; dim_i++ ) {
			pos[dim_i] = self->min[dim_i] + self->step[dim_i] * gInds[dim_i];
		}

		/* Dummy interpolation. */
		erpRes = FieldVariable_InterpolateValueAt( self->velVar, pos, dummyVel );
		if( erpRes == LOCAL ) {
			local = Memory_Realloc_Array( local, unsigned, nLocal + 1 );
			local[nLocal++] = gNode_i;
		}
		else if( erpRes == SHADOW ) {
			shadow = Memory_Realloc_Array( shadow, unsigned, nShadow + 1 );
			shadow[nShadow++] = gNode_i;
		}
		else if( erpRes != OTHER_PROC ) {
			assert( 0 );
		}
	}
	FreeArray( gInds );
	FreeArray( dummyVel );

	/* There is a good chance there will be some nodes that noone owns, for the moment the lowest order
	   proc for each shared node will get it. */
	MPI_ArrayAllgather( nShadow, sizeof(unsigned), shadow, 
			    &remNShad, (void***)&remShad, 
			    self->communicator );
	MPI_Comm_rank( self->communicator, (int*)&rank );
	MPI_Comm_size( self->communicator, (int*)&nProcs );
	setIntersect( remShad, remNShad, 
		      &intr, &nIntr, 
		      nProcs, rank );
	step = nIntr / nProcs;
	local = Memory_Realloc_Array( local, unsigned, nLocal + step + (((nIntr % nProcs) > rank) ? 1 : 0) );
	for( intr_i = rank; intr_i < nIntr; intr_i += step ) {
		local[nLocal++] = intr[intr_i];
	}
	FreeArray( intr );

	self->nodeSync = Sync_New( "LevelSet" );
	Sync_Negotiate( self->nodeSync, 
			nGNodes, 
			local, nLocal, 
			NULL, 0, 
			shadow, nShadow, 
			self->communicator );


	/*
	** Create the mapping between elements and nodes, ensuring the ordering required for linear
	** finite elements is observed.
	*/

	base = Memory_Alloc_Array( unsigned, self->grm.nDims, "LevelSet" );
	self->nElNodes = (unsigned)pow( 2.0, (double)self->grm.nDims );
	elNodes = Memory_Alloc_Array( unsigned, self->nElNodes, "LevelSet" );
	self->baseNodeEl = Memory_Alloc_Array( unsigned, Sync_GetNLocal( self->nodeSync ), "LevelSet" );
	for( lNode_i = 0; lNode_i < Sync_GetNLocal( self->nodeSync ); lNode_i++ ) {
		unsigned	curElNode;
		unsigned	gNodeInd;

		gNodeInd = Sync_MapLocal( self->nodeSync, lNode_i );
		GRM_Lift( &self->grm, gNodeInd, base );
		curElNode = 0;
		if( _LevelSet_ValidateElNodes( &self->grm, self->nodeSync, 
					       base, self->grm.nDims, 
					       elNodes, &curElNode ) )
		{
			/* This is a local element, store it. */
			self->elNodes = Memory_Realloc_Array( self->elNodes, unsigned*, self->nEls + 1 );
			self->elNodes[self->nEls] = Memory_Alloc_Array_Unnamed( unsigned, self->nElNodes );
			memcpy( self->elNodes[self->nEls], elNodes, sizeof(unsigned) * self->nElNodes );
			self->baseNodeEl[lNode_i] = self->nEls;
			self->nEls++;
		}
	}


	/*
	** Calculate the jacobian constant.  We can do this because the level set mesh is always geometrically
	** regular and uses linear shape functions.
	*/

	self->jacDet = 1.0;
	for( dim_i = 0; dim_i < self->dim; dim_i++ ) {
		self->jacDet *= 0.5 * self->step[dim_i];
	}


	/*
	** Calculate quadrature weights and abscissas and apply to shape functions.  
	** I'm going to use 3 points just to be sure.
	*/

	self->qdr = Memory_Alloc_Array_Unnamed( Quadrature, 1 );
	Quadrature_Init( self->qdr, self->dim, 4 );
	Quadrature_GaussLegendre( self->qdr, -1.0, 1.0 );


	/*
	** Build the shape functions and associated tables.
	*/

	/* Construct linear shape functions for each element node. */
	self->shapeFuncs = Memory_Alloc_Array_Unnamed( Function, self->nElNodes );
	for( sf_i = 0; sf_i < self->nElNodes; sf_i++ ) {
		Function*	sf = self->shapeFuncs + sf_i;
		unsigned	bit;

		sf->nTerms = 1;
		sf->nDims = self->dim;
		sf->polys = Memory_Alloc_Array_Unnamed( Polynomial*, self->shapeFuncs[sf_i].nTerms );
		sf->polys[0] = Memory_Alloc_Array_Unnamed( Polynomial, self->shapeFuncs[sf_i].nDims );
		for( dim_i = 0; dim_i < self->dim; dim_i++ ) {
			Polynomial*	poly = sf->polys[0] + dim_i;

			bit = extractBit( sf_i, dim_i );

			poly->order = 1;
			poly->coefs = Memory_Alloc_Array_Unnamed( double, poly->order + 1 );
			poly->coefs[0] = 0.5;
			poly->coefs[1] = (bit == 0) ? -0.5 : 0.5;
		}
	}

	/* Build a table of integrated psi values. */
	memset( &func1, 0, sizeof(Function) );
	crds = Memory_Alloc_Array_Unnamed( double, self->dim );
	self->psiTbl = Memory_Alloc_Array_Unnamed( double, self->nElNodes );
	for( node_i = 0; node_i < self->nElNodes; node_i++ ) {
		self->psiTbl[node_i] = Quadrature_Eval( self->qdr, self->shapeFuncs + node_i );
		self->psiTbl[node_i] *= self->jacDet;
	}

	/* Build a table of integrated psi-squared values. */
	memset( &func2, 0, sizeof(Function) );
	self->psiSqTbl = Memory_Alloc_Array_Unnamed( double*, self->nElNodes );
	for( node_i = 0; node_i < self->nElNodes; node_i++ ) {
		self->psiSqTbl[node_i] = Memory_Alloc_Array_Unnamed( double, self->nElNodes );
		for( node_j = 0; node_j < self->nElNodes; node_j++ ) {
			Function_Mult( self->shapeFuncs + node_i, self->shapeFuncs + node_j, &func1 );
			self->psiSqTbl[node_i][node_j] = Quadrature_Eval( self->qdr, &func1 );
			self->psiSqTbl[node_i][node_j] *= self->jacDet;
		}
	}

	/* Build a table of integrated psi-grad-psi terms. */
	memset( &func3, 0, sizeof(Function) );
	funcVec = Memory_Alloc_Array_Unnamed( Function, self->dim );
	memset( funcVec, 0, sizeof(Function) * self->dim );
	self->gradTbl = Memory_Alloc_Array_Unnamed( double***, self->nElNodes );
	for( node_i = 0; node_i < self->nElNodes; node_i++ ) {
		self->gradTbl[node_i] = Memory_Alloc_Array_Unnamed( double**, self->nElNodes );
		for( node_j = 0; node_j < self->nElNodes; node_j++ ) {
			self->gradTbl[node_i][node_j] = Memory_Alloc_Array_Unnamed( double*, self->nElNodes );
			Function_Grad( self->shapeFuncs + node_j, funcVec );
			for( node_k = 0; node_k < self->nElNodes; node_k++ ) {
				self->gradTbl[node_i][node_j][node_k] = Memory_Alloc_Array_Unnamed( double, self->dim );
				Function_Mult( self->shapeFuncs + node_i, self->shapeFuncs + node_k, &func1 );
				for( dim_i = 0; dim_i < self->dim; dim_i++ ) {
					Function_Mult( &func1, funcVec + dim_i, &func2 );
					self->gradTbl[node_i][node_j][node_k][dim_i] = Quadrature_Eval( self->qdr, &func2 );
					self->gradTbl[node_i][node_j][node_k][dim_i] *= self->jacDet;

					/* Due to changes in variable (to the parent domain)... */
					self->gradTbl[node_i][node_j][node_k][dim_i] *= (2.0 / self->step[dim_i]);
				}
			}
		}
	}

	/* Free stuff. */
	FreeArray( crds );
	Function_Delete( &func1 );
	Function_Delete( &func2 );

	/* Build quadrature values for all shape functions. */
	crds = Memory_Alloc_Array_Unnamed( double, self->dim );
	self->psiQdr = Memory_Alloc_Array_Unnamed( double*, self->nElNodes );
	for( node_i = 0; node_i < self->nElNodes; node_i++ ) {
		self->psiQdr[node_i] = Memory_Alloc_Array_Unnamed( double, self->qdr->nPnts );
		for( qdr_i = 0; qdr_i < self->qdr->nPnts; qdr_i++ ) {
			/* Evaluate function. */
			self->psiQdr[node_i][qdr_i] = Function_Eval( self->shapeFuncs + node_i, self->qdr->x[qdr_i] );
		}
	}

	/* Build quadrature values for all gradient shape functions. */
	self->gradQdr = Memory_Alloc_Array_Unnamed( double**, self->nElNodes );
	for( node_i = 0; node_i < self->nElNodes; node_i++ ) {
		self->gradQdr[node_i] = Memory_Alloc_Array_Unnamed( double*, self->dim );
		Function_Grad( self->shapeFuncs + node_i, funcVec );
		for( dim_i = 0; dim_i < self->dim; dim_i++ ) {
			self->gradQdr[node_i][dim_i] = Memory_Alloc_Array_Unnamed( double, self->qdr->nPnts );
			for( qdr_i = 0; qdr_i < self->qdr->nPnts; qdr_i++ ) {
				/* Evaluate function. */
				self->gradQdr[node_i][dim_i][qdr_i] = Function_Eval( funcVec + dim_i, self->qdr->x[qdr_i] );

				/* Due to changes in variable (to the parent domain)... */
				self->gradQdr[node_i][dim_i][qdr_i] *= (2.0 / self->step[dim_i]);
			}
		}
	}

	/* Free stuff. */
	for( dim_i = 0; dim_i < self->dim; dim_i++ ) {
		Function_Delete( funcVec + dim_i );
	}
	FreeArray( funcVec );


	/*
	** Setup internal arrays for distance functions and velocities.
	*/

	self->dists = Memory_Alloc_Array( double, Sync_GetNLocal( self->nodeSync ), "LevelSet" );
	self->vels = Memory_Alloc_Array( double, Sync_GetNLocal( self->nodeSync ) * self->grm.nDims, "LevelSet" );


	/*
	** Allocate space for the matrix and two vectors and create the solver.
	*/

	self->aMat = Matrix_New( self->nodeSync->comm, 
				 Sync_GetNLocal( self->nodeSync ), Sync_GetNLocal( self->nodeSync ), 
				 (unsigned)pow( 3.0, (double)self->grm.nDims ) );
	self->rhsVec = Vector_New_SpecifyLocalSize( self->nodeSync->comm, 
						    Sync_GetNLocal( self->nodeSync ) );
	self->solVec = Vector_New_SpecifyLocalSize( self->nodeSync->comm, 
						    Sync_GetNLocal( self->nodeSync ) );
	self->solver = MatrixSolver_Build( self->nodeSync->comm, self->aMat );
	MatrixSolver_Setup( self->solver, self->aMat );

	/* Calculate min and max distances. */
	_LevelSet_CalcMinMax( self );
}


void LevelSet_InitDistances( void* levelSet, 
			     double (func)( const unsigned* inds, const double* pos, void* ctx ), void* funcCtx )
{
	LevelSet*	self = (LevelSet*)levelSet;
	unsigned*	inds;
	double*		pos;
	double*		step;
	unsigned	lNode_i;
	unsigned	dim_i;

	assert( levelSet );
	assert( func );
	assert( self->min && self->max );
	assert( self->nodeSync );


	/*
	** Loop over all the nodes of the mesh, sending each to the callback function which will give us the distance from the
	** surface at that point.
	*/

	inds = Memory_Alloc_Array( unsigned, self->grm.nDims, "LevelSet" );
	pos = Memory_Alloc_Array( double, self->grm.nDims, "LevelSet" );
	step = Memory_Alloc_Array( double, self->grm.nDims, "LevelSet" );
	for( dim_i = 0; dim_i < self->grm.nDims; dim_i++ ) {
		step[dim_i] = (self->max[dim_i] - self->min[dim_i]) / (double)(self->grm.nNodes[dim_i] - 1);
	}

	for( lNode_i = 0; lNode_i < Sync_GetNLocal( self->nodeSync ); lNode_i++ ) {
		unsigned	gNodeInd;

		gNodeInd = Sync_MapLocal( self->nodeSync, lNode_i );
		GRM_Lift( &self->grm, gNodeInd, inds );
		for( dim_i = 0; dim_i < self->grm.nDims; dim_i++ ) {
			pos[dim_i] = self->min[dim_i] + (double)inds[dim_i] * step[dim_i];
		}

		self->dists[lNode_i] = func( inds, pos, funcCtx );
	}

	Memory_Free( inds );
	Memory_Free( pos );
}


void LevelSet_Advance( void* levelSet, double dt ) {
	LevelSet*	self = (LevelSet*)levelSet;
	double*		pos;
	double**	velQdr;
	double*		imv;
	double*		uwConst;
	unsigned*	inds;
	unsigned	nIts;
	unsigned	fac_i, lNode_i;

	assert( self );
	assert( self->nodeSync );
	assert( self->velVar );

	if( dt < 1e-20 ) {
		return;
	}


	/*
	** Interpolate velocity from the supplied velocity variable onto all nodes of the level set mesh.
	*/

	inds = Memory_Alloc_Array( unsigned, self->grm.nDims, "LevelSet" );
	pos = Memory_Alloc_Array( double, self->grm.nDims, "LevelSet" );

	for( lNode_i = 0; lNode_i < Sync_GetNLocal( self->nodeSync ); lNode_i++ ) {
		InterpolationResult	res;
		unsigned		gNodeInd;
		unsigned		dim_i;

		/* Calculate the node's position. */
		gNodeInd = Sync_MapLocal( self->nodeSync, lNode_i );
		GRM_Lift( &self->grm, gNodeInd, inds );
		for( dim_i = 0; dim_i < self->grm.nDims; dim_i++ ) {
			pos[dim_i] = self->min[dim_i] + inds[dim_i] * self->step[dim_i];
		}

		/* Interpolate velocity. */
		res = FieldVariable_InterpolateValueAt( self->velVar, pos, self->vels + lNode_i * self->grm.nDims );
		for( dim_i = 0; dim_i < self->grm.nDims; dim_i++ ) {
			self->vels[lNode_i * self->grm.nDims + dim_i] *= -1.0;
		}
		if( res == OTHER_PROC || res == OUTSIDE_GLOBAL ) {
			assert( 0 );
		}
	}


	/*
	** Construct the matrix and RHS of the problem. Here we also need to solve 'fac' times to
	** account for the difference in dt scale.
	*/

	/* Allocate some space for upwinding terms. */
	velQdr = Memory_Alloc_2DArray_Unnamed( double, self->qdr->nPnts, self->dim );
	imv = Memory_Alloc_Array_Unnamed( double, self->qdr->nPnts );
	uwConst = Memory_Alloc_Array_Unnamed( double, self->qdr->nPnts );

	dt /= (double)self->fac;
	for( fac_i = 0; fac_i < self->fac; fac_i++ ) {
		unsigned	aEl_i;

		Matrix_Zero( self->aMat );
		Vector_Zero( self->rhsVec );
		Vector_Zero( self->solVec );

		/* Loop over local active elements to build the system. */
		for( aEl_i = 0; aEl_i < self->nActEls; aEl_i++ ) {
			unsigned	lEl_i = self->actEls[aEl_i];
			unsigned	node_i, qdr_i, dim_i;

			/* For the upwinding, calculate velocity at quadrature points. */
			for( qdr_i = 0; qdr_i < self->qdr->nPnts; qdr_i++ ) {
				for( dim_i = 0; dim_i < self->dim; dim_i++ ) {
					velQdr[qdr_i][dim_i] = 0.0;
					for( node_i = 0; node_i < self->nElNodes; node_i++ ) {
						unsigned	lNodeI = self->elNodes[lEl_i][node_i];

						velQdr[qdr_i][dim_i] += self->psiQdr[node_i][qdr_i] * 
							self->vels[lNodeI * self->dim + dim_i];
					}
				}
			}

			/* For the upwinding, calculate the inverse velocity magnitude. */
			memset( imv, 0, sizeof(double) * self->qdr->nPnts );
			for( qdr_i = 0; qdr_i < self->qdr->nPnts; qdr_i++ ) {
				for( dim_i = 0; dim_i < self->dim; dim_i++ ) {
					imv[qdr_i] += velQdr[qdr_i][dim_i] * velQdr[qdr_i][dim_i];
				}
				imv[qdr_i] = 1.0 / imv[qdr_i];
			}

			/* For the upwinding, calculate the artificial diffusivity constant. */
			memset( uwConst, 0, sizeof(double) * self->qdr->nPnts );
			for( qdr_i = 0; qdr_i < self->qdr->nPnts; qdr_i++ ) {
				for( dim_i = 0; dim_i < self->dim; dim_i++ ) {
					double	mag = (velQdr[qdr_i][dim_i] < 0.0) ? -velQdr[qdr_i][dim_i] : velQdr[qdr_i][dim_i];

					uwConst[qdr_i] += mag * self->step[dim_i];
				}
				uwConst[qdr_i] *= 0.5;
			}

			/* Begin filling matrix and RHS. */
			for( node_i = 0; node_i < self->nElNodes; node_i++ ) {
				unsigned	lNodeI;
				unsigned	eqNumI;
				double		val;
				unsigned	node_j;

				/* Get the equation number. */
				lNodeI = self->elNodes[lEl_i][node_i];
				eqNumI = Sync_MapLocal( self->nodeSync, lNodeI );

				/* Here is where we build this element's part of the matrix. */
				for( node_j = 0; node_j < self->nElNodes; node_j++ ) {
					unsigned	lNodeJ = self->elNodes[lEl_i][node_j];
					unsigned	eqNumJ = Sync_MapLocal( self->nodeSync, lNodeJ );
					double*		vel = self->vels + lNodeJ * self->dim;
					double**	gradTbl = self->gradTbl[node_i][node_j];
					double		val2 = 0.0;
					unsigned	node_k;

					/* First integral. */
					val = self->psiSqTbl[node_i][node_j];

					/* Second integral. */
					for( node_k = 0; node_k < self->nElNodes; node_k++ ) {
						for( dim_i = 0; dim_i < self->dim; dim_i++ ) {
							val2 += vel[dim_i] * gradTbl[node_k][dim_i];
						}
					}
					val -= dt * val2;

					/* Upwinded part of the second integral (advection). */
					val2 = 0.0;
					for( qdr_i = 0; qdr_i < self->qdr->nPnts; qdr_i++ ) {
						double	tmp1 = 0.0;
						double	tmp2 = 0.0;

						for( dim_i = 0; dim_i < self->dim; dim_i++ ) {
							tmp1 += velQdr[qdr_i][dim_i] * self->gradQdr[node_i][dim_i][qdr_i];
							tmp2 += velQdr[qdr_i][dim_i] * self->gradQdr[node_j][dim_i][qdr_i];
						}

						val2 += self->qdr->w[qdr_i] * uwConst[qdr_i] * imv[qdr_i] * tmp1 * tmp2;
					}
					val -= dt * val2;

					/* Add to matrix. */
					Matrix_AddValue( self->aMat, eqNumI, eqNumJ, val );

					/* RHS. */
					val = self->dists[lNodeJ] * self->psiSqTbl[node_i][node_j];
					Vector_AddEntry( self->rhsVec, eqNumI, val );
				}
			}
		}

		_LevelSet_InsertInactive( self );

		Matrix_AssemblyBegin( self->aMat );
		Vector_AssemblyBegin( self->rhsVec );
		Vector_AssemblyBegin( self->solVec );
		Matrix_AssemblyEnd( self->aMat );
		Vector_AssemblyEnd( self->rhsVec );
		Vector_AssemblyEnd( self->solVec );


		/*
		** Solve linear system and update distances.
		*/

		nIts = MatrixSolver_Solve( self->solver, self->solVec, self->rhsVec );
		Vector_GetValues( self->solVec, Sync_GetNLocal( self->nodeSync ), self->nodeSync->local, self->dists );


		/*
		** Redistance and extend.
		*/

		_LevelSet_Redistance( self, 10 );
		//_LevelSet_InitActiveEls( self );
	}

	/* Free all the quadrature memory. */
	FreeArray( velQdr );
	FreeArray( imv );
	FreeArray( uwConst );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void _LevelSet_Free( LevelSet* self ) {
	KillArray( self->min );
	KillArray( self->max );
	KillArray2D( self->nEls, self->elNodes );
	KillArray( self->baseNodeEl );
	KillArray( self->dists );
	KillArray( self->vels );
	KillArray( self->shapeFuncs );
}


unsigned ipow( unsigned base, unsigned exp ) {
	unsigned	res;
	unsigned	e_i;

	res = 1;
	for( e_i = 0; e_i < exp; e_i++ ) {
		res *= base;
	}

	return res;
}


void _LevelSet_Redistance( LevelSet* self, unsigned maxIts ) {
	double		dt;
	double*		vec;
	unsigned	it_i;

	dt = self->step[0] * self->step[0] * 0.5;
	vec = Memory_Alloc_Array_Unnamed( double, self->dim );
	for( it_i = 0; it_i < maxIts; it_i++ ) {
		unsigned	nIts;
		unsigned	aEl_i;

		Matrix_Zero( self->aMat );
		Vector_Zero( self->rhsVec );
		Vector_Zero( self->solVec );

		for( aEl_i = 0; aEl_i < self->nActEls; aEl_i++ ) {
			unsigned	lEl_i = self->actEls[aEl_i];
			unsigned	node_i;

			for( node_i = 0; node_i < self->nElNodes; node_i++ ) {
				unsigned	lNodeI, eqNumI;
				double		val;
				unsigned	qdr_i, node_j;

				/* Get the equation number. */
				lNodeI = self->elNodes[lEl_i][node_i];
				eqNumI = Sync_MapLocal( self->nodeSync, lNodeI );

				/* Update the matrix. */
				for( node_j = 0; node_j < self->nElNodes; node_j++ ) {
					unsigned	lNodeJ;
					unsigned	eqNumJ;

					lNodeJ = self->elNodes[lEl_i][node_j];
					eqNumJ = Sync_MapLocal( self->nodeSync, lNodeJ );
					Matrix_AddValue( self->aMat, eqNumI, eqNumJ, self->psiSqTbl[node_i][node_j] );

					/* Build one of the RHS integrals while we're here. */
					Vector_AddEntry( self->rhsVec, eqNumI, 
							 self->dists[lNodeJ] * self->psiSqTbl[node_i][node_j] );
				}

				/* Handle the RHS. */
				val = 0.0;
				for( qdr_i = 0; qdr_i < self->qdr->nPnts; qdr_i++ ) {
					double		sgn, mag;
					unsigned	dim_i, node_j;

					/* Evaluate nodal sums. */
					sgn = 0.0;
					memset( vec, 0, sizeof(double) * self->dim );
					for( node_j = 0; node_j < self->nElNodes; node_j++ ) {
						unsigned	lNodeJ = self->elNodes[lEl_i][node_j];
						double		dist = self->dists[lNodeJ];

						sgn += dist * self->psiQdr[node_j][qdr_i];
						for( dim_i = 0; dim_i < self->dim; dim_i++ ) {
							vec[dim_i] += dist * self->gradQdr[node_j][dim_i][qdr_i];
						}
					}
					sgn = (sgn < 0.0) ? -1.0 : 1.0;
					mag = 0.0;
					for( dim_i = 0; dim_i < self->dim; dim_i++ ) {
						mag += vec[dim_i] * vec[dim_i];
					}
					mag = sqrt( mag );

					/* Update value. */
					val += self->qdr->w[qdr_i] * self->psiQdr[node_i][qdr_i] * sgn * (1.0 - mag);
				}

				/* Update vector. */
				val *= self->jacDet;
				Vector_AddEntry( self->rhsVec, eqNumI, dt * val );
			}
		}

		_LevelSet_InsertInactive( self );

		/* Assemble. */
		Matrix_AssemblyBegin( self->aMat );
		Vector_AssemblyBegin( self->rhsVec );
		Vector_AssemblyBegin( self->solVec );
		Matrix_AssemblyEnd( self->aMat );
		Vector_AssemblyEnd( self->rhsVec );
		Vector_AssemblyEnd( self->solVec );

		/* Solve. */
		nIts = MatrixSolver_Solve( self->solver, self->solVec, self->rhsVec );
		Vector_GetValues( self->solVec, Sync_GetNLocal( self->nodeSync ), self->nodeSync->local, self->dists );
	}
	FreeArray( vec );
}


void _LevelSet_InsertInactive( LevelSet* self ) {
	unsigned	iNode_i;

	assert( self );

	/* Insert the inactive nodes. */
	for( iNode_i = 0; iNode_i < self->nInactNodes; iNode_i++ ) {
		unsigned	lNodeI = self->inactNodes[iNode_i];
		unsigned	eqNumI = Sync_MapLocal( self->nodeSync, lNodeI );

		Matrix_AddValue( self->aMat, eqNumI, eqNumI, 1.0 );
		Vector_AddEntry( self->rhsVec, eqNumI, self->dists[lNodeI] );
	}
}


void _LevelSet_InitActiveEls( LevelSet* self ) {
	unsigned	i;
	self->nInactNodes = 0;
	self->nActEls = self->nEls;
	self->actEls = Memory_Realloc_Array( self->actEls, unsigned, self->nActEls );
	for( i = 0; i < self->nActEls; i++ ) {
		self->actEls[i] = i;
	}
#if 0
	IndexSet*	actNodes;
	unsigned	lEl_i;

	assert( self );
	assert( self->dists );

	/*
	** Determine which elements are close enough to the surface to be involved in calculations.
	*/

	/* Clear existing arrays. */
	KillArray( self->actEls );
	KillArray( self->inactNodes );
	self->nActEls = 0;
	self->nInactNodes = 0;

	actNodes = IndexSet_New( Sync_GetNLocal( self->nodeSync ) );

	for( lEl_i = 0; lEl_i < self->nEls; lEl_i++ ) {
		unsigned	intNodes = 0;
		unsigned	node_i;

		for( node_i = 0; node_i < self->nElNodes; node_i++ ) {
			unsigned	lNodeI = self->elNodes[lEl_i][node_i];

			if( self->dists[lNodeI] <= self->cutoff ) {
				intNodes++;
			}
		}

		/* Both internal and boundary are added to the active list. */
		if( intNodes ) {
			self->actEls = Memory_Realloc_Array( self->actEls, unsigned, self->nActEls + 1 );
			self->actEls[self->nActEls] = lEl_i;
			self->nActEls++;

			/* Add the nodes to the active node set. */
			for( node_i = 0; node_i < self->nElNodes; node_i++ ) {
				IndexSet_Add( actNodes, self->elNodes[lEl_i][node_i] );
			}
		}
	}

	IndexSet_GetVacancies( actNodes, &self->nInactNodes, &self->inactNodes );
	Stg_Class_Delete( actNodes );
#endif
}


#if 0
Bool _LevelSet_Extend( LevelSet* self ) {
	Bool		checkInt = False;
	unsigned	nActEls = 0;
	unsigned	nBndEls = 0;
	unsigned	nIntEls = 0;
	unsigned*	actEls = NULL;
	unsigned*	bndEls = NULL;
	unsigned*	intEls = NULL;
	unsigned	el_i;

	assert( self );

	/*
	** Check if any of the elements of the level set have come within a certain range of
	** the active element boundary, and if so extend the active elements.
	*/

	for( el_i = 0; el_i < self->nBndEls; el_i++ ) {
		unsigned	intNodes = 0;
		unsigned	lElI = self->bndEls[el_i];

		for( node_i = 0; node_i < self->nElNodes; node_i++ ) {
			unsigned	lNodeI = self->elNodes[lElI][node_i];

			if( self->dists[lNodeI] <= self->cutoff ) {
				intNodes++;
			}
		}

		if( intNodes == self->nElNodes ) {
			/* Boundary element has become an internal element. */
			checkInt = True;
		}
		else if( intNodes ) {
			/* Boundary element remains so. */
		}
		else {
			/* Boundary element has become external. */
			checkInt = True;
		}
	}

	/*
	** Check the interior elements if needed.
	*/

	if( checkInt ) {
		for( el_i = 0; el_i < self->nIntEls; el_i++ ) {
			unsigned	intNodes = 0;
			unsigned	lElI = self->intEls[el_i];

			for( node_i = 0; node_i < self->nElNodes; node_i++ ) {
				unsigned	lNodeI = self->elNodes[lElI][node_i];

				if( self->dists[lNodeI] <= self->cutoff ) {
					intNodes++;
				}
			}

			if( intNodes == self->nElNodes ) {
				/* Interior element remains so. */
			}
			else if( intNodes ) {
				/* Interior element has become a boundary element. */
			}
			else {
				/* Interior element has become an external element. */
			}
		}
	}
}
#endif


Bool _LevelSet_ValidateElNodes( GRM* grm, Sync* nodeSync, 
				unsigned* base, unsigned short curDim, 
				unsigned* elNodes, unsigned* curElNode )
{
	if( curDim > 0 ) {
		/* Don't attempt this on maximum indices. */
		if( base[curDim - 1] == grm->nNodes[curDim - 1] - 1 ) {
			return False;
		}

		if( !_LevelSet_ValidateElNodes( grm, nodeSync, 
						base, curDim - 1, 
						elNodes, curElNode ) )
		{
			return False;
		}
		else {
			Bool	res;

			base[curDim - 1]++;
			res = _LevelSet_ValidateElNodes( grm, nodeSync, 
							 base, curDim - 1, 
							 elNodes, curElNode );
			base[curDim - 1]--;

			return res;
		}
	}
	else {
		GRM_Project( grm, base, elNodes + *curElNode );
		if( !Sync_IsLocal( nodeSync, elNodes[*curElNode] ) ) {
			return False;
		}
		else {
			elNodes[*curElNode] = Sync_MapGlobal( nodeSync, elNodes[*curElNode] );
			(*curElNode)++;
		}

		return True;
	}
}


void setIntersect( unsigned** sets, unsigned* nEls, 
		   unsigned** intr, unsigned* nIntr, 
		   unsigned nProcs, unsigned rank )
{
	unsigned*	own;
	unsigned	nOwn;
	unsigned*	rem;
	unsigned	nRem;
	unsigned	ownInd;
	unsigned	p_i, own_i, rem_i;

	assert( !nOwn || own );
	assert( sets );
	assert( nEls );
	assert( intr );
	assert( nIntr );
	assert( nProcs >= 1 );
	assert( rank < nProcs );

	own = sets[rank];
	nOwn = nEls[rank];
	*intr = NULL;
	*nIntr = 0;
	for( own_i = 0; own_i < nOwn; own_i++ ) {
		ownInd = own[own_i];
		for( p_i = 0; p_i < nProcs; p_i++ ) {
			assert( !nEls[p_i] || sets[p_i] );
			if( p_i == rank ) {
				continue;
			}

			rem = sets[p_i];
			nRem = nEls[p_i];
			for( rem_i = 0; rem_i < nRem; rem_i++ ) {
				if( ownInd == rem[rem_i] ) {
					break;
				}
			}
			if( rem_i == nEls[p_i] ) {
				break;
			}
		}
		if( p_i == nProcs ) {
			*intr = Memory_Realloc_Array( *intr, unsigned, *nIntr + 1 );
			(*intr)[(*nIntr)++] = ownInd;
		}
	}
}


void _LevelSet_CalcMinMax( LevelSet* self ) {
	unsigned	dim_i;

	assert( self );

	self->minDist = self->dists[0];
	self->maxDist = self->dists[0];
	for( dim_i = 1; dim_i < Sync_GetNLocal( self->nodeSync ); dim_i++ ) {
		if( self->dists[dim_i] < self->minDist ) {
			self->minDist = self->dists[dim_i];
		}

		if( self->dists[dim_i] > self->maxDist ) {
			self->maxDist = self->dists[dim_i];
		}
	}
}
