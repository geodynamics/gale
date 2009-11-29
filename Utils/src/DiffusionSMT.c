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
** $Id: DiffusionSMT.c 964 2007-10-11 08:03:06Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <assert.h>
#include <string.h>
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PopulationControl/PopulationControl.h>
#include <PICellerator/Weights/Weights.h>
#include <PICellerator/MaterialPoints/MaterialPoints.h>
#include "types.h"
#include "DiffusionSMT.h"
#include "MaterialSwarmVariable.h"


/* Textual name of this class */
const Type DiffusionSMT_Type = "DiffusionSMT";


DiffusionSMT* DiffusionSMT_New( 
    Name                                                name,
    FiniteElementContext*				                    context,
    StiffnessMatrix*                                    stiffnessMatrix,
    Swarm*                                              integrationSwarm )
{
    DiffusionSMT* self = (DiffusionSMT*) _DiffusionSMT_DefaultNew( name );

	self->isConstructed = True;
	_StiffnessMatrixTerm_Init( self, context, stiffnessMatrix, integrationSwarm, NULL );
	_DiffusionSMT_Init( self );

	return self;
}

/* Creation implementation / Virtual constructor */
DiffusionSMT* _DiffusionSMT_New(  DIFFUSIONSMT_DEFARGS  )
{
    DiffusionSMT* self;
	
    /* Allocate memory */
    assert( _sizeOfSelf >= sizeof(DiffusionSMT) );
    self = (DiffusionSMT*) _StiffnessMatrixTerm_New(  STIFFNESSMATRIXTERM_PASSARGS  );
	
    /* Virtual info */
	
    return self;
}

void _DiffusionSMT_Init( void* matrixTerm ) {
    DiffusionSMT* self = (DiffusionSMT*)matrixTerm;
}

void _DiffusionSMT_Delete( void* matrixTerm ) {
    DiffusionSMT* self = (DiffusionSMT*)matrixTerm;
    int ii;

    for(ii = 0; ii < self->materialSwarmCount; ii++)
	Stg_Class_Delete( self->diffusionSwarmVariables[ii] );

    _StiffnessMatrixTerm_Delete( self );
}

void _DiffusionSMT_Print( void* matrixTerm, Stream* stream ) {
    DiffusionSMT* self = (DiffusionSMT*)matrixTerm;
	
    _StiffnessMatrixTerm_Print( self, stream );

    /* General info */
}

void* _DiffusionSMT_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                 _sizeOfSelf = sizeof(DiffusionSMT);
	Type                                                         type = DiffusionSMT_Type;
	Stg_Class_DeleteFunction*                                 _delete = _DiffusionSMT_Delete;
	Stg_Class_PrintFunction*                                   _print = _DiffusionSMT_Print;
	Stg_Class_CopyFunction*                                     _copy = NULL;
	Stg_Component_DefaultConstructorFunction*     _defaultConstructor = _DiffusionSMT_DefaultNew;
	Stg_Component_ConstructFunction*                       _construct = _DiffusionSMT_AssignFromXML;
	Stg_Component_BuildFunction*                               _build = _DiffusionSMT_Build;
	Stg_Component_InitialiseFunction*                     _initialise = _DiffusionSMT_Initialise;
	Stg_Component_ExecuteFunction*                           _execute = _DiffusionSMT_Execute;
	Stg_Component_DestroyFunction*                           _destroy = _DiffusionSMT_Destroy;
	StiffnessMatrixTerm_AssembleElementFunction*     _assembleElement = _DiffusionSMT_AssembleElement;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;

    return (void*)_DiffusionSMT_New(  DIFFUSIONSMT_PASSARGS  );
}

void _DiffusionSMT_AssignFromXML( void* matrixTerm, Stg_ComponentFactory* cf, void* data ) {
    DiffusionSMT*            self             = (DiffusionSMT*)matrixTerm;
    PICelleratorContext*     context	      = (PICelleratorContext*)self->context;

    /* Construct Parent */
    _StiffnessMatrixTerm_AssignFromXML( self, cf, data );

    _DiffusionSMT_Init( self );

    assert( Stg_CheckType( context, PICelleratorContext ) );
    self->materials_Register = context->materials_Register;
}

void _DiffusionSMT_Build( void* matrixTerm, void* data ) {
    DiffusionSMT*             self             = (DiffusionSMT*)matrixTerm;
    AbstractContext*                 context;
    Stg_ComponentFactory*            cf;
    Materials_Register*              materials_Register = self->materials_Register;
    IntegrationPointsSwarm*          swarm = (IntegrationPointsSwarm*)self->integrationSwarm;
    MaterialPointsSwarm**            materialSwarms;
    DiffusionSMT_MaterialExt*   materialExt;
    Name                             name;
    Material*                        material;
    int ii;

    _StiffnessMatrixTerm_Build( self, data );

    /* Get Component Factory if we can */
    if ( Stg_Class_IsInstance( data, AbstractContext_Type ) ) {
	context = (AbstractContext*) data;
	cf = context->CF;
    }

    /* Sort out material extension stuff */
    self->materialExtHandle = Materials_Register_AddMaterialExtension(
	self->materials_Register, 
	self->type, 
	sizeof(DiffusionSMT_MaterialExt));

    for ( ii = 0 ; ii < Materials_Register_GetCount( materials_Register ) ; ii++) {
	material = Materials_Register_GetByIndex( materials_Register, ii );
	materialExt = ExtensionManager_GetFunc( material->extensionMgr, material,
						self->materialExtHandle );

	if ( cf ) {
	    materialExt->diffusion = Stg_ComponentFactory_GetDouble(
		cf, material->name, "diffusivity", 0.0 );
	}
	else {
	    materialExt->diffusion = Dictionary_GetDouble_WithDefault(
		material->dictionary, "diffusivity", 0.0 );
	}
    }

    /* Create Swarm Variables of each material swarm this ip swarm is mapped against */
    materialSwarms = IntegrationPointMapper_GetMaterialPointsSwarms(
	swarm->mapper, &(self->materialSwarmCount) );
    self->diffusionSwarmVariables = Memory_Alloc_Array(
	MaterialSwarmVariable*, self->materialSwarmCount, "DiffusionVariables");

    for ( ii = 0; ii < self->materialSwarmCount; ++ii ) {
	name = Stg_Object_AppendSuffix( materialSwarms[ii], "diffusivity" );
	self->diffusionSwarmVariables[ii] = MaterialSwarmVariable_New( 
	    name, 
	    (AbstractContext*) self->context, 
	    materialSwarms[ii], 
	    1, 
	    self->materials_Register, 
	    self->materialExtHandle, 
	    GetOffsetOfMember( *materialExt, diffusion ) );
	Memory_Free( name );
	
	/* Build new Swarm Variables */
	Stg_Component_Build( self->diffusionSwarmVariables[ii], data, False );
    }
}

void _DiffusionSMT_Initialise( void* matrixTerm, void* data ) {
    DiffusionSMT*             self             = (DiffusionSMT*)matrixTerm;
    int ii;

    _StiffnessMatrixTerm_Initialise( self, data );

    for ( ii = 0; ii < self->materialSwarmCount; ++ii )
	Stg_Component_Initialise( self->diffusionSwarmVariables[ii], data, False );
}

void _DiffusionSMT_Execute( void* matrixTerm, void* data ) {
    _StiffnessMatrixTerm_Execute( matrixTerm, data );
}

void _DiffusionSMT_Destroy( void* matrixTerm, void* data ) {
    _StiffnessMatrixTerm_Destroy( matrixTerm, data );
}


void _DiffusionSMT_AssembleElement( 
    void*                                              matrixTerm,
    StiffnessMatrix*                                   stiffnessMatrix, 
    Element_LocalIndex                                 lElement_I, 
    SystemLinearEquations*                             sle,
    FiniteElementContext*                              context,
    double**                                           elStiffMat ) 
{
    DiffusionSMT*       self         = Stg_CheckType( matrixTerm, DiffusionSMT );
    Swarm*                              swarm        = self->integrationSwarm;
    FeVariable*                         variable1    = stiffnessMatrix->rowVariable;
    Dimension_Index                     dim          = stiffnessMatrix->dim;
    IntegrationPoint*                   currIntegrationPoint;
    double*                             xi;
    double                              weight;
    Particle_InCellIndex                cParticle_I, cellParticleCount;
    Index                               nodesPerEl;
    Index                               i,j;
    Dimension_Index                     dim_I;
    double**                            GNx;
    double                              detJac;
	
    Cell_Index                          cell_I;
    ElementType*                        elementType;
    DiffusionSMT_MaterialExt*   materialExt;
    Material*                        material;
	
    /* Set the element type */
    elementType = FeMesh_GetElementType( variable1->feMesh, lElement_I );
    nodesPerEl = elementType->nodeCount;
	
    /* allocate */
    GNx = Memory_Alloc_2DArray( double, dim, nodesPerEl, "GNx" );
	
    cell_I = CellLayout_MapElementIdToCellId( swarm->cellLayout, lElement_I );
    cellParticleCount = swarm->cellParticleCountTbl[ cell_I ];

    /* Slap the laplacian together */
    for( cParticle_I = 0 ; cParticle_I < cellParticleCount ; cParticle_I++ ) {
		
	currIntegrationPoint = (IntegrationPoint*)Swarm_ParticleInCellAt(
	    swarm, cell_I, cParticle_I );

	material = IntegrationPointsSwarm_GetMaterialOn(
	    (IntegrationPointsSwarm*) swarm, currIntegrationPoint );
	materialExt = ExtensionManager_Get(
	    material->extensionMgr, material, self->materialExtHandle );

	xi = currIntegrationPoint->xi;
	weight = currIntegrationPoint->weight;
		
	ElementType_ShapeFunctionsGlobalDerivs( 
	    elementType,
	    variable1->feMesh, lElement_I,
	    xi, dim, &detJac, GNx );

	for( i=0; i<nodesPerEl; i++ ) {
	    for( j=0; j<nodesPerEl; j++ ) {
		for ( dim_I = 0; dim_I < dim ; dim_I++ ) { 
		    elStiffMat[i][j] += 
			materialExt->diffusion * detJac * weight * 
			GNx[dim_I][i] * GNx[dim_I][j];
		}
	    }
	}
    }
	
    Memory_Free(GNx); 
}


