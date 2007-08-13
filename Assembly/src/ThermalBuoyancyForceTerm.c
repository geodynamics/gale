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
** $Id: ThermalBuoyancyForceTerm.c 936 2007-08-13 00:47:42Z DavidLee $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/Discretisation/Discretisation.h>
#include <StgFEM/SLE/SLE.h>

#include "types.h"
#include "ThermalBuoyancyForceTerm.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class */
const Type ThermalBuoyancyForceTerm_Type = "ThermalBuoyancyForceTerm";

ThermalBuoyancyForceTerm* ThermalBuoyancyForceTerm_New( 
		Name                                                name,
		ForceVector*                                        forceVector,
		Swarm*                                              integrationSwarm,
		FeVariable*                                         temperatureField,
		double                                              rayleighNumber )
{
	ThermalBuoyancyForceTerm* self = (ThermalBuoyancyForceTerm*) _ThermalBuoyancyForceTerm_DefaultNew( name );

	ThermalBuoyancyForceTerm_InitAll( 
			self,
			forceVector,
			integrationSwarm,
			temperatureField,
			rayleighNumber );

	return self;
}

/* Creation implementation / Virtual constructor */
ThermalBuoyancyForceTerm* _ThermalBuoyancyForceTerm_New( 
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
		ForceTerm_AssembleElementFunction*                   _assembleElement,		
		Name                                                name )
{
	ThermalBuoyancyForceTerm* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(ThermalBuoyancyForceTerm) );
	self = (ThermalBuoyancyForceTerm*) _ForceTerm_New( 
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

void _ThermalBuoyancyForceTerm_Init( 
		ThermalBuoyancyForceTerm*                                    self, 
		FeVariable*                                         temperatureField,
		double                                              rayleighNumber )
{
	self->temperatureField    = temperatureField;
	self->rayleighNumber      = rayleighNumber;
}

void ThermalBuoyancyForceTerm_InitAll( 
		void*                                               forceTerm,
		ForceVector*                                        forceVector,
		Swarm*                                              integrationSwarm,
		FeVariable*                                         temperatureField,
		double                                              rayleighNumber )
{
	ThermalBuoyancyForceTerm* self = (ThermalBuoyancyForceTerm*) forceTerm;

	ForceTerm_InitAll( self, forceVector, integrationSwarm, NULL );
	_ThermalBuoyancyForceTerm_Init( self, temperatureField, rayleighNumber );
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
	return (void*)_ThermalBuoyancyForceTerm_New( 
		sizeof(ThermalBuoyancyForceTerm), 
		ThermalBuoyancyForceTerm_Type,
		_ThermalBuoyancyForceTerm_Delete,
		_ThermalBuoyancyForceTerm_Print,
		NULL,
		_ThermalBuoyancyForceTerm_DefaultNew,
		_ThermalBuoyancyForceTerm_Construct,
		_ThermalBuoyancyForceTerm_Build,
		_ThermalBuoyancyForceTerm_Initialise,
		_ThermalBuoyancyForceTerm_Execute,
		_ThermalBuoyancyForceTerm_Destroy,
		_ThermalBuoyancyForceTerm_AssembleElement,
		name );
}

void _ThermalBuoyancyForceTerm_Construct( void* forceTerm, Stg_ComponentFactory* cf, void* data ) {
	ThermalBuoyancyForceTerm*            self             = (ThermalBuoyancyForceTerm*)forceTerm;
	FeVariable*                 temperatureField;
	double                      rayleighNumber;

	/* Construct Parent */
	_ForceTerm_Construct( self, cf, data );

	temperatureField = Stg_ComponentFactory_ConstructByKey( cf, self->name, "TemperatureField", FeVariable, True, data ) ;
	rayleighNumber   = Stg_ComponentFactory_GetDouble( cf, self->name, "Ra", 0.0 );

	_ThermalBuoyancyForceTerm_Init( self, temperatureField, rayleighNumber );
}

void _ThermalBuoyancyForceTerm_Build( void* forceTerm, void* data ) {
	ThermalBuoyancyForceTerm*             self             = (ThermalBuoyancyForceTerm*)forceTerm;

	_ForceTerm_Build( self, data );

	Stg_Component_Build( self->temperatureField, data, False );
}

void _ThermalBuoyancyForceTerm_Initialise( void* forceTerm, void* data ) {
	ThermalBuoyancyForceTerm*             self             = (ThermalBuoyancyForceTerm*)forceTerm;

	_ForceTerm_Initialise( self, data );

	Stg_Component_Initialise( self->temperatureField, data, False );
}

void _ThermalBuoyancyForceTerm_Execute( void* forceTerm, void* data ) {
	_ForceTerm_Execute( forceTerm, data );
}

void _ThermalBuoyancyForceTerm_Destroy( void* forceTerm, void* data ) {
	_ForceTerm_Destroy( forceTerm, data );
}


void _ThermalBuoyancyForceTerm_AssembleElement( void* forceTerm, ForceVector* forceVector, Element_LocalIndex lElement_I, double* elForceVec ) {
	ThermalBuoyancyForceTerm*        self               = Stg_CheckType( forceTerm, ThermalBuoyancyForceTerm );
	Swarm*                  swarm              = self->integrationSwarm;
	Dimension_Index         dim                = forceVector->dim;
	IntegrationPoint*       particle;
	FeVariable*             temperatureField;
	FeMesh*                 mesh;
	FeMesh*                 temperatureMesh;
	double*                 xi;
	Particle_InCellIndex    cParticle_I;
	Particle_InCellIndex    cellParticleCount;
	Element_NodeIndex       elementNodeCount;
	Node_ElementLocalIndex  node_I;
	ElementType*            elementType;
	Dof_Index               dofsPerNode;
	Cell_Index              cell_I;
	double                  detJac;
	double                  factor;
	/*double                  Ni[8];*/
	double                  Ni[27];
	double                  force;
	double                  rayleighNumber;
	double                  temperature;
	double                  tmpGlobalCoord[3];

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
		if ( temperatureMesh == mesh ) {
			/* If meshes are identical - then we can use a shortcut for working out the temperature */
			FeVariable_InterpolateWithinElement( temperatureField, lElement_I, xi, &temperature );
		}
		else {
			FeMesh_CoordLocalToGlobal( mesh, lElement_I, xi, tmpGlobalCoord );
			FieldVariable_InterpolateValueAt( temperatureField, tmpGlobalCoord, &temperature );
		}

		force = rayleighNumber * temperature;

		factor = detJac * particle->weight * force;
		for( node_I = 0 ; node_I < elementNodeCount ; node_I++ ) 
			elForceVec[node_I * dofsPerNode + J_AXIS ] += factor * Ni[ node_I ] ;
		
	}
}

