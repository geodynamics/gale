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
** $Id: ThermalBuoyancyForceTerm.c 964 2007-10-11 08:03:06Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/Discretisation/Discretisation.h>
#include <StgFEM/SLE/SLE.h>

#include "types.h"
#include "ThermalBuoyancyForceTerm.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class */
const Type ThermalBuoyancyForceTerm_Type = "ThermalBuoyancyForceTerm";

ThermalBuoyancyForceTerm* ThermalBuoyancyForceTerm_New( 
	Name							name,
	FiniteElementContext*	context,
	ForceVector*				forceVector,
	Swarm*						integrationSwarm,
	FeVariable*					temperatureField,
	double						rayleighNumber )
{
	ThermalBuoyancyForceTerm* self = (ThermalBuoyancyForceTerm*) _ThermalBuoyancyForceTerm_DefaultNew( name );

	self->isConstructed = True;
	_ForceTerm_Init( self, context, forceVector, integrationSwarm, NULL );
	_ThermalBuoyancyForceTerm_Init( self, temperatureField, rayleighNumber );

	return self;
}

/* Creation implementation / Virtual constructor */
ThermalBuoyancyForceTerm* _ThermalBuoyancyForceTerm_New(  THERMALBUOYANCYFORCETERM_DEFARGS  )
{
	ThermalBuoyancyForceTerm* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(ThermalBuoyancyForceTerm) );
	/* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
	/* This means that any values of these parameters that are passed into this function are not passed onto the parent function
	   and so should be set to ZERO in any children of this class. */
	nameAllocationType = NON_GLOBAL;

	self = (ThermalBuoyancyForceTerm*) _ForceTerm_New(  FORCETERM_PASSARGS  );
	
	/* Virtual info */
	
	return self;
}

void _ThermalBuoyancyForceTerm_Init( void* forceTerm, FeVariable* temperatureField, double rayleighNumber ) {
	ThermalBuoyancyForceTerm* self = (ThermalBuoyancyForceTerm*)forceTerm;

	self->temperatureField = temperatureField;
	self->rayleighNumber = rayleighNumber;
}

void _ThermalBuoyancyForceTerm_Delete( void* forceTerm ) {
	ThermalBuoyancyForceTerm* self = (ThermalBuoyancyForceTerm*)forceTerm;

	_ForceTerm_Delete( self );
}

void _ThermalBuoyancyForceTerm_Print( void* forceTerm, Stream* stream ) {
	ThermalBuoyancyForceTerm* self = (ThermalBuoyancyForceTerm*)forceTerm;
	
	_ForceTerm_Print( self, stream );

	/* General info */
	Journal_PrintPointer( stream, self->temperatureField );
	Journal_PrintDouble( stream, self->rayleighNumber );
}

void* _ThermalBuoyancyForceTerm_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(ThermalBuoyancyForceTerm);
	Type                                                      type = ThermalBuoyancyForceTerm_Type;
	Stg_Class_DeleteFunction*                              _delete = _ThermalBuoyancyForceTerm_Delete;
	Stg_Class_PrintFunction*                                _print = _ThermalBuoyancyForceTerm_Print;
	Stg_Class_CopyFunction*                                  _copy = NULL;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _ThermalBuoyancyForceTerm_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _ThermalBuoyancyForceTerm_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _ThermalBuoyancyForceTerm_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _ThermalBuoyancyForceTerm_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _ThermalBuoyancyForceTerm_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _ThermalBuoyancyForceTerm_Destroy;
	ForceTerm_AssembleElementFunction*            _assembleElement = _ThermalBuoyancyForceTerm_AssembleElement;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;

	return (void*)_ThermalBuoyancyForceTerm_New(  THERMALBUOYANCYFORCETERM_PASSARGS  );
}

void _ThermalBuoyancyForceTerm_AssignFromXML( void* forceTerm, Stg_ComponentFactory* cf, void* data ) {
	ThermalBuoyancyForceTerm*	self = (ThermalBuoyancyForceTerm*)forceTerm;
	FeVariable*						temperatureField;
	double							rayleighNumber;

	/* Construct Parent */
	_ForceTerm_AssignFromXML( self, cf, data );

	temperatureField = Stg_ComponentFactory_ConstructByKey( cf, self->name, "TemperatureField", FeVariable, True, data ) ;
	rayleighNumber   = Stg_ComponentFactory_GetDouble( cf, self->name, "Ra", 0.0 );

	_ThermalBuoyancyForceTerm_Init( self, temperatureField, rayleighNumber );
}

void _ThermalBuoyancyForceTerm_Build( void* forceTerm, void* data ) {
	ThermalBuoyancyForceTerm* self = (ThermalBuoyancyForceTerm*)forceTerm;

	Stg_Component_Build( self->temperatureField, data, False );
	_ForceTerm_Build( self, data );
}

void _ThermalBuoyancyForceTerm_Initialise( void* forceTerm, void* data ) {
	ThermalBuoyancyForceTerm* self = (ThermalBuoyancyForceTerm*)forceTerm;

	Stg_Component_Initialise( self->temperatureField, data, False );
	_ForceTerm_Initialise( self, data );
}

void _ThermalBuoyancyForceTerm_Execute( void* forceTerm, void* data ) {
	ThermalBuoyancyForceTerm* self = (ThermalBuoyancyForceTerm*)forceTerm;

	_ForceTerm_Execute( self, data );
}

void _ThermalBuoyancyForceTerm_Destroy( void* forceTerm, void* data ) {
	ThermalBuoyancyForceTerm* self = (ThermalBuoyancyForceTerm*)forceTerm;

   Stg_Component_Destroy( self->temperatureField, data, False );
	_ForceTerm_Destroy( self, data );
}


void _ThermalBuoyancyForceTerm_AssembleElement( void* forceTerm, ForceVector* forceVector, Element_LocalIndex lElement_I, double* elForceVec ) {
	ThermalBuoyancyForceTerm*	self = Stg_CheckType( forceTerm, ThermalBuoyancyForceTerm );
	Swarm*							swarm = self->integrationSwarm;
	Dimension_Index				dim = forceVector->dim;
	IntegrationPoint*				particle;
	FeVariable*						temperatureField;
	FeMesh*							mesh;
	FeMesh*							temperatureMesh;
	double*							xi;
	Particle_InCellIndex			cParticle_I;
	Particle_InCellIndex			cellParticleCount;
	Element_NodeIndex				elementNodeCount;
	Node_ElementLocalIndex		node_I;
	ElementType*					elementType;
	Dof_Index						dofsPerNode;
	Cell_Index						cell_I;
	double							detJac;
	double							factor;
	/*double							Ni[8];*/
	double							Ni[27];
	double							force;
	double							rayleighNumber;
	double							temperature;
	double							tmpGlobalCoord[3];

	/* Get context extension */
	rayleighNumber   = self->rayleighNumber;
	temperatureField = self->temperatureField;
	temperatureMesh  = temperatureField->feMesh;

	/* Since we are integrating over the velocity mesh - we want the velocity mesh here and not the temperature mesh */
	mesh             = forceVector->feVariable->feMesh;
	
	/* Set the element type */
	elementType      = FeMesh_GetElementType( mesh, lElement_I ); 
	elementNodeCount = elementType->nodeCount;

	/* assumes constant number of dofs per element */
	dofsPerNode = dim;
	
	cell_I = CellLayout_MapElementIdToCellId( swarm->cellLayout, lElement_I );
	cellParticleCount = swarm->cellParticleCountTbl[ cell_I ];
	
	for ( cParticle_I = 0 ; cParticle_I < cellParticleCount ; cParticle_I++ ) {
		particle = (IntegrationPoint*) Swarm_ParticleInCellAt( swarm, cell_I, cParticle_I );
		xi       = particle->xi;
		
		/* Calculate Determinant of Jacobian and Shape Functions */
		detJac = ElementType_JacobianDeterminant( elementType, mesh, lElement_I, xi, dim );
		ElementType_EvaluateShapeFunctionsAt( elementType, xi, Ni );

		/* Field Get Temperature from Field Variable */
		FeVariable_InterpolateFromMeshLocalCoord( temperatureField, mesh, lElement_I, xi, &temperature );

		force = rayleighNumber * temperature;

		factor = detJac * particle->weight * force;
		for( node_I = 0 ; node_I < elementNodeCount ; node_I++ ) 
			elForceVec[node_I * dofsPerNode + J_AXIS ] += factor * Ni[ node_I ] ;
		
	}
}



