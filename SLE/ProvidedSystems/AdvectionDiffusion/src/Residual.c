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
** $Id: Residual.c 920 2007-07-20 06:19:34Z RobertTurnbull $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/Discretisation/Discretisation.h>
#include <StgFEM/SLE/LinearAlgebra/LinearAlgebra.h>
#include <StgFEM/SLE/SystemSetup/SystemSetup.h>

#include "types.h"
#include "AdvectionDiffusionSLE.h"
#include "Residual.h"
#include "ShapeFunctions.h"
#include "UpwindParameter.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class */
const Type AdvDiffResidualForceTerm_Type = "AdvDiffResidualForceTerm";

AdvDiffResidualForceTerm* AdvDiffResidualForceTerm_New( 
		Name                                                name,
		ForceVector*                                        forceVector,
		Swarm*                                              integrationSwarm,
		Stg_Component*                                      sle, 
		FeVariable*                                         velocityField,
		Variable*                                           diffusivityVariable,
		double                                              defaultDiffusivity,
		AdvDiffResidualForceTerm_UpwindParamFuncType        upwindFuncType )
{
	AdvDiffResidualForceTerm* self = (AdvDiffResidualForceTerm*) _AdvDiffResidualForceTerm_DefaultNew( name );

	AdvDiffResidualForceTerm_InitAll( 
			self,
			forceVector,
			integrationSwarm,
			sle,
			velocityField,
			diffusivityVariable,
			defaultDiffusivity,
			upwindFuncType );

	return self;
}

/* Creation implementation / Virtual constructor */
AdvDiffResidualForceTerm* _AdvDiffResidualForceTerm_New( 
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
		ForceTerm_AssembleElementFunction*                  _assembleElement,		
		AdvDiffResidualForceTerm_UpwindParamFunction*       _upwindParam,
		Name                                                name )
{
	AdvDiffResidualForceTerm* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(AdvDiffResidualForceTerm) );
	self = (AdvDiffResidualForceTerm*) _ForceTerm_New( 
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
	self->_upwindParam = _upwindParam;
	
	return self;
}

void _AdvDiffResidualForceTerm_Init( 
		AdvDiffResidualForceTerm*                           self, 
		FeVariable*                                         velocityField,
		Variable*                                           diffusivityVariable,
		double                                              defaultDiffusivity,
		AdvDiffResidualForceTerm_UpwindParamFuncType        upwindFuncType ) //WHY IS THIS LINE HERE???
{
	self->velocityField       = velocityField;
	self->diffusivityVariable = diffusivityVariable;
	self->defaultDiffusivity  = defaultDiffusivity;
}

void AdvDiffResidualForceTerm_InitAll( 
		void*                                               residual,
		ForceVector*                                        forceVector,
		Swarm*                                              integrationSwarm,
		Stg_Component*                                      sle, 
		FeVariable*                                         velocityField,
		Variable*                                           diffusivityVariable,
		double                                              defaultDiffusivity,
		AdvDiffResidualForceTerm_UpwindParamFuncType        upwindFuncType )
{
	AdvDiffResidualForceTerm* self = (AdvDiffResidualForceTerm*) residual;

	ForceTerm_InitAll( self, forceVector, integrationSwarm, sle );
	_AdvDiffResidualForceTerm_Init( self, velocityField, diffusivityVariable, defaultDiffusivity, upwindFuncType );
}

void _AdvDiffResidualForceTerm_Delete( void* residual ) {
	AdvDiffResidualForceTerm* self = (AdvDiffResidualForceTerm*)residual;

	_ForceTerm_Delete( self );
}

void _AdvDiffResidualForceTerm_Print( void* residual, Stream* stream ) {
	AdvDiffResidualForceTerm* self = (AdvDiffResidualForceTerm*)residual;
	
	_ForceTerm_Print( self, stream );

	Journal_Printf( stream, "self->calculateUpwindParam = %s\n", 
			self->_upwindParam == AdvDiffResidualForceTerm_UpwindXiExact ? 
				"AdvDiffResidualForceTerm_UpwindXiExact" :
			self->_upwindParam == AdvDiffResidualForceTerm_UpwindXiDoublyAsymptoticAssumption ? 
				"AdvDiffResidualForceTerm_UpwindXiDoublyAsymptoticAssumption" :
			self->_upwindParam == AdvDiffResidualForceTerm_UpwindXiCriticalAssumption ? 
				"AdvDiffResidualForceTerm_UpwindXiCriticalAssumption" : "Unknown"  );

	/* General info */
	Journal_PrintPointer( stream, self->velocityField );
	Journal_PrintDouble( stream, self->defaultDiffusivity );
	Journal_Printf( stream, "self->diffusivityVariable = ");
	if ( self->diffusivityVariable )
		Journal_Printf( stream, "%s\n", self->diffusivityVariable->name );
	else
		Journal_Printf( stream, "<Unused>\n");

}

void* _AdvDiffResidualForceTerm_DefaultNew( Name name ) {
	return (void*)_AdvDiffResidualForceTerm_New( 
		sizeof(AdvDiffResidualForceTerm), 
		AdvDiffResidualForceTerm_Type,
		_AdvDiffResidualForceTerm_Delete,
		_AdvDiffResidualForceTerm_Print,
		NULL,
		_AdvDiffResidualForceTerm_DefaultNew,
		_AdvDiffResidualForceTerm_Construct,
		_AdvDiffResidualForceTerm_Build,
		_AdvDiffResidualForceTerm_Initialise,
		_AdvDiffResidualForceTerm_Execute,
		_AdvDiffResidualForceTerm_Destroy,
		_AdvDiffResidualForceTerm_AssembleElement,
		_AdvDiffResidualForceTerm_UpwindParam,
		name );
}

void _AdvDiffResidualForceTerm_Construct( void* residual, Stg_ComponentFactory* cf, void* data ) {
	AdvDiffResidualForceTerm*            self             = (AdvDiffResidualForceTerm*)residual;
	FeVariable*                          velocityField;
	Variable*                            diffusivityVariable;
	Name                                 upwindParamFuncName;
	double                               defaultDiffusivity;
	AdvDiffResidualForceTerm_UpwindParamFuncType  upwindFuncType       = 0;

	/* Construct Parent */
	_ForceTerm_Construct( self, cf, data );

	velocityField       = Stg_ComponentFactory_ConstructByKey( cf, self->name, "VelocityField",       FeVariable, True,  data );
	diffusivityVariable = Stg_ComponentFactory_ConstructByKey( cf, self->name, "DiffusivityVariable", Variable,   False, data );


	upwindParamFuncName = Stg_ComponentFactory_GetString( cf, self->name, "UpwindXiFunction", "Exact" );
	if ( strcasecmp( upwindParamFuncName, "DoublyAsymptoticAssumption" ) == 0 )
		upwindFuncType = DoublyAsymptoticAssumption;
	else if ( strcasecmp( upwindParamFuncName, "CriticalAssumption" ) == 0 )
		upwindFuncType = CriticalAssumption;
	else if ( strcasecmp( upwindParamFuncName, "Exact" ) == 0 )
		upwindFuncType = Exact;
	else 
		Journal_Firewall( False, Journal_Register( Error_Type, self->type ), 
				"Cannot understand '%s'\n", upwindParamFuncName );

	defaultDiffusivity = Stg_ComponentFactory_GetDouble( cf, self->name, "defaultDiffusivity", 1.0 );

	_AdvDiffResidualForceTerm_Init( self, velocityField, diffusivityVariable, defaultDiffusivity, upwindFuncType );
}

void _AdvDiffResidualForceTerm_Build( void* residual, void* data ) {
	AdvDiffResidualForceTerm*             self             = (AdvDiffResidualForceTerm*)residual;

	_ForceTerm_Build( self, data );

	Stg_Component_Build( self->velocityField, data, False );
	if ( self->diffusivityVariable )
		Stg_Component_Build( self->diffusivityVariable, data, False );
}

void _AdvDiffResidualForceTerm_Initialise( void* residual, void* data ) {
	AdvDiffResidualForceTerm*             self             = (AdvDiffResidualForceTerm*)residual;

	_ForceTerm_Initialise( self, data );

	Stg_Component_Initialise( self->velocityField, data, False );
	if ( self->diffusivityVariable )
		Stg_Component_Initialise( self->diffusivityVariable, data, False );
}

void _AdvDiffResidualForceTerm_Execute( void* residual, void* data ) {
	_ForceTerm_Execute( residual, data );
}

void _AdvDiffResidualForceTerm_Destroy( void* residual, void* data ) {
	_ForceTerm_Destroy( residual, data );
}


void _AdvDiffResidualForceTerm_AssembleElement( void* forceTerm, ForceVector* forceVector, Element_LocalIndex lElement_I, double* elementResidual ) {
	AdvDiffResidualForceTerm*  self               = Stg_CheckType( forceTerm, AdvDiffResidualForceTerm );
	AdvectionDiffusionSLE*     sle                = Stg_CheckType( self->extraInfo, AdvectionDiffusionSLE );
	Swarm*                     swarm              = self->integrationSwarm;
	Particle_Index             lParticle_I;
	Particle_Index             cParticle_I;
	Particle_Index             cellParticleCount;
	Cell_Index                 cell_I;    
	IntegrationPoint*          particle;
	FeVariable*                phiField           = sle->phiField;
	Dimension_Index            dim                = forceVector->dim;
	double                     velocity[3];
	double                     phi, phiDot;
	double*                    phiGrad;
	double                     detJac;
	double**                   GNx;
	double**                   elementUpwindShapeFunc;
	double*                    upwindShapeFunc;
	double*                    xi;
	double                     totalDerivative, diffusionTerm;
	double                     diffusivity         = self->defaultDiffusivity;
	Variable*                  diffusivityVariable = self->diffusivityVariable;
	ElementType*               elementType         = FeMesh_GetElementType( phiField->feMesh, lElement_I );
	Node_Index                 elementNodeCount    = elementType->nodeCount;
	Node_Index                 node_I;
	double                     factor;

	GNx     = Memory_Alloc_2DArray( double, dim, elementNodeCount, "Global Shape Function Derivatives" );
	phiGrad = Memory_Alloc_Array( double, dim, "Gradient of Phi" );

	elementUpwindShapeFunc = AdvDiffResidualForceTerm_BuildSUPGShapeFunctions( self, sle, swarm, lElement_I, dim );

	/* Determine number of particles in element */
	cell_I = CellLayout_MapElementIdToCellId( swarm->cellLayout, lElement_I );
	cellParticleCount = swarm->cellParticleCountTbl[ cell_I ];
	
	for ( cParticle_I = 0 ; cParticle_I < cellParticleCount ; cParticle_I++ ) {
		lParticle_I     = swarm->cellParticleTbl[cell_I][cParticle_I];

		particle        = (IntegrationPoint*) Swarm_ParticleAt( swarm, lParticle_I );
		xi              = particle->xi;
		
		/* Evalutate Shape Functions */
		upwindShapeFunc = elementUpwindShapeFunc[ cParticle_I ];

		/* Calculate Global Shape Function Derivatives */
		ElementType_ShapeFunctionsGlobalDerivs( 
			elementType,
			phiField->feMesh, lElement_I,
			xi, dim, &detJac, GNx );
		
		/* Calculate Velocity */
		if ( phiField->feMesh == self->velocityField->feMesh ) {
			/* If the phi field's mesh and the velocity field's mesh are identical - then we can assume we are in the same
			 * element and have the same local coordinate for this integration point */
			FeVariable_InterpolateWithinElement( self->velocityField, lElement_I, xi, velocity );	
		}
		else {
			Coord globalCoord;

			/* Since the phi field's mesha and the velocity field - we have to get the velocity the hard way */
			FeMesh_CoordLocalToGlobal( phiField->feMesh, lElement_I, xi, globalCoord );
			FieldVariable_InterpolateValueAt( self->velocityField, globalCoord, velocity );
		}

		/* Calculate phi on particle */
		_FeVariable_InterpolateNodeValuesToElLocalCoord( phiField, lElement_I, xi, &phi );

		/* Calculate Gradients of Phi */
		FeVariable_InterpolateDerivatives_WithGNx( phiField, lElement_I, GNx, phiGrad );

		/* Calculate time derivative of phi */
		_FeVariable_InterpolateNodeValuesToElLocalCoord( sle->phiDotField, lElement_I, xi, &phiDot );
		
		/* Calculate total derivative (i.e. Dphi/Dt = \dot \phi + u . \grad \phi) */
		totalDerivative = phiDot + StGermain_VectorDotProduct( velocity, phiGrad, dim );

		/* Get Diffusivity */
		/* diffusivityVariable will only be NOT NULL if:
		 * 1) The MaterialDiffusivityPlugin is used. It's in Underworld/Plugins
		 * 2) A special user defined DiffusivityVariable is given during the Construction phase
		 */  
		if ( diffusivityVariable != NULL )
			diffusivity = self->_getDiffusivityFromIntPoint( self, particle );

		/* Add to element residual */
		factor = particle->weight * detJac;
		for ( node_I = 0 ; node_I < elementNodeCount ; node_I++ ) {
			/* Calculate Diffusion Term */
			diffusionTerm = diffusivity * ( GNx[0][node_I] * phiGrad[0] + GNx[1][node_I] * phiGrad[1] );
			if (dim == 3)
				diffusionTerm += diffusivity * GNx[2][ node_I ] * phiGrad[2] ;
			
			elementResidual[ node_I ] -= factor * ( upwindShapeFunc[ node_I ] * totalDerivative + diffusionTerm );
		}
	}
	
	Memory_Free( elementUpwindShapeFunc );
	Memory_Free( phiGrad );
	Memory_Free( GNx );
}


/* Virtual Function Implementations */
double _AdvDiffResidualForceTerm_UpwindParam( void* residual, double pecletNumber ) {
	AdvDiffResidualForceTerm*             self             = (AdvDiffResidualForceTerm*)residual;

	switch ( self->upwindParamType ) {
		case Exact:
			self->_upwindParam = AdvDiffResidualForceTerm_UpwindXiExact; break;
		case DoublyAsymptoticAssumption:
			self->_upwindParam = AdvDiffResidualForceTerm_UpwindXiDoublyAsymptoticAssumption; break;
		case CriticalAssumption:
			self->_upwindParam = AdvDiffResidualForceTerm_UpwindXiCriticalAssumption; break;
	}

	return AdvDiffResidualForceTerm_UpwindParam( self, pecletNumber );
}
