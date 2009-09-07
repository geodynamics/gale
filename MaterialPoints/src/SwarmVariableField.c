/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
** $Id: SwarmVariableField.c 376 2006-10-18 06:58:41Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PopulationControl/PopulationControl.h>
#include <PICellerator/Weights/Weights.h>

#include "types.h"
#include "ParticleFeVariable.h"
#include "SwarmVariableField.h"
#include "IntegrationPointsSwarm.h"
#include "MaterialPointsSwarm.h"
#include "IntegrationPointMapper.h"

#include <assert.h>
#include <string.h>

const Type SwarmVariableField_Type = "SwarmVariableField";

/* Creation implementation / Virtual constructor */
SwarmVariableField* _SwarmVariableField_New( 
		SizeT                                             _sizeOfSelf,  
		Type                                              type,
		Stg_Class_DeleteFunction*                         _delete,
		Stg_Class_PrintFunction*                          _print,
		Stg_Class_CopyFunction*                           _copy, 
		Stg_Component_DefaultConstructorFunction*         _defaultConstructor,
		Stg_Component_ConstructFunction*                  _construct,
		Stg_Component_BuildFunction*                      _build,
		Stg_Component_InitialiseFunction*                 _initialise,
		Stg_Component_ExecuteFunction*                    _execute,
		Stg_Component_DestroyFunction*                    _destroy,
		FieldVariable_InterpolateValueAtFunction*         _interpolateValueAt,
		FieldVariable_GetValueFunction*	                  _getMinGlobalFeMagnitude,
		FieldVariable_GetValueFunction*                   _getMaxGlobalFeMagnitude,
		FieldVariable_GetCoordFunction*             	  _getMinAndMaxLocalCoords,
		FieldVariable_GetCoordFunction*             	  _getMinAndMaxGlobalCoords,
		FeVariable_InterpolateWithinElementFunction*      _interpolateWithinElement,	
		FeVariable_GetValueAtNodeFunction*	  	  _getValueAtNode,
		/*SwarmVariableField_GetValueAtNodeFunction*	  _getValueAtNode,*/
		/*SwarmVariableField_ValueAtParticleFunction*       _valueAtParticle,*/
		ParticleFeVariable_ValueAtParticleFunction*       _valueAtParticle,
		Name                                              name ) 
{
	SwarmVariableField* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(SwarmVariableField) );

	self = (SwarmVariableField*) _ParticleFeVariable_New(
 						_sizeOfSelf,
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
						_interpolateValueAt,
						_getMinGlobalFeMagnitude,
						_getMaxGlobalFeMagnitude,
						_getMinAndMaxLocalCoords,
						_getMinAndMaxGlobalCoords,		
						_interpolateWithinElement,	
						_getValueAtNode,
						_valueAtParticle,
						name );
	/* Virtual info */
	
	return self;
}

void _SwarmVariableField_Init( SwarmVariableField* swarmVariableField, SwarmVariable* swarmVar, Variable_Register* variable_Register ) {
	SwarmVariableField* 	self 	= (SwarmVariableField*)swarmVariableField;
	self->swarmVar = swarmVar;
	self->variable_Register = variable_Register;
}

void _SwarmVariableField_Delete( void* swarmVariableField ) {
	SwarmVariableField* self = (SwarmVariableField*)swarmVariableField;

	_ParticleFeVariable_Delete( self );
}

void _SwarmVariableField_Print( void* swarmVariableField, Stream* stream ) {
	SwarmVariableField* self = (SwarmVariableField*)swarmVariableField;
	
	_ParticleFeVariable_Print( self, stream );
}

void* _SwarmVariableField_DefaultNew( Name name ) {
	return (void*)_SwarmVariableField_New( 
			sizeof(SwarmVariableField), 
			SwarmVariableField_Type,
			_SwarmVariableField_Delete,
			_SwarmVariableField_Print,
			NULL,
			_SwarmVariableField_DefaultNew,
			_SwarmVariableField_Construct,
			_SwarmVariableField_Build,
			_SwarmVariableField_Initialise,
			_SwarmVariableField_Execute,
			_SwarmVariableField_Destroy,
			_FeVariable_InterpolateValueAt,
			_FeVariable_GetMinGlobalFieldMagnitude,
			_FeVariable_GetMaxGlobalFieldMagnitude,
			_FeVariable_GetMinAndMaxLocalCoords,
			_FeVariable_GetMinAndMaxGlobalCoords,
			_FeVariable_InterpolateNodeValuesToElLocalCoord,
			_FeVariable_GetValueAtNode,
			/*_SwarmVariableField_GetValueAtNode,*/
			_SwarmVariableField_ValueAtParticle,
			name );
}

void _SwarmVariableField_Construct( void* swarmVariableField, Stg_ComponentFactory* cf, void* data ) {
	SwarmVariableField*     	self		= (SwarmVariableField*)swarmVariableField;
	SwarmVariable*			swarmVar;
	IntegrationPointsSwarm* 	integrationSwarm;
	Variable_Register*		variable_Register;

	variable_Register = self->context->variable_Register; 

	/* Construct Parent */
	_ParticleFeVariable_Construct( self, cf, data );

	// TODO: just get the textual name here - see gLucifer's SwarmPlotter DrawignObject 
	self->swarmVarName = Stg_ComponentFactory_GetString( cf, self->name, "swarmVariable", "" );
	assert( swarmVar );

	self->materialSwarm = Stg_ComponentFactory_ConstructByKey( cf, self->name, "MaterialSwarm", MaterialPointsSwarm, True, data );

	integrationSwarm = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Swarm", IntegrationPointsSwarm, True, NULL );
	assert( integrationSwarm );

	_SwarmVariableField_Init( self, swarmVar, variable_Register );
	_ParticleFeVariable_Init( self, integrationSwarm, self->context );
}

void _SwarmVariableField_Build( void* swarmVariableField, void* data ) {
	SwarmVariableField*	self	= (SwarmVariableField*)swarmVariableField;
	Name			tmpName;
	unsigned		nDomainVerts	= Mesh_GetDomainSize( self->feMesh, MT_VERTEX );

	/* make this more flexible to handle vector values at each node - will have to get the num dofs from the XML
	 * as other components are not necessarily built yet... dave. 03.10.07 */
	assert( Class_IsSuper( self->feMesh->topo, IGraph ) );
	tmpName = Stg_Object_AppendSuffix( self, "DataVariable" );
	self->dataVariable = Variable_NewScalar( tmpName,
		       				 Variable_DataType_Double,
						 &nDomainVerts,
						 NULL,
						 (void**)&self->data,
						 self->variable_Register );
	Memory_Free( tmpName );
	self->fieldComponentCount = 1;

	tmpName = Stg_Object_AppendSuffix( self, "DofLayout" );
	self->dofLayout = DofLayout_New( tmpName, self->variable_Register, 0, self->feMesh );

	/* must build before adding the variable to the dof layout, dave. 04.10.07 */
	//Stg_Component_Build( self->dofLayout, data, False );
	Stg_Component_Build( self->dofLayout->mesh, data, False );
	self->dofLayout->_numItemsInLayout = Mesh_GetDomainSize( self->dofLayout->mesh, MT_VERTEX );
	DofLayout_AddAllFromVariableArray( self->dofLayout, 1, &self->dataVariable );
	Memory_Free( tmpName );
	
	self->eqNum->dofLayout = self->dofLayout;

	_ParticleFeVariable_Build( self, data );
	/* TODO: build self->swarmVar */
	// TODO: granb self->SwarmVariableName out of the swarm vart register, save reference as self->swarmVar
	//_SwarmVariable_Build( self->swarmVar, data ); /* dave, 18.09.07 */
}

void _SwarmVariableField_Initialise( void* swarmVariableField, void* data ) {
	SwarmVariableField*		self			= (SwarmVariableField*)swarmVariableField;
	SwarmVariable_Register*		swarmVar_Register	= self->materialSwarm->swarmVariable_Register;
	Stream*				errorStream		= Journal_Register( Error_Type, self->type );

	/* assign the swarm variable, assume its already built */
	if( 0 != strcmp( self->swarmVarName, "" ) ) {
		self->swarmVar = SwarmVariable_Register_GetByName( swarmVar_Register, self->swarmVarName );
		Journal_Firewall( self->swarmVar != NULL, errorStream, "Error - cannot find swarm variable \"%s\" in %s() for swarm \"%s\".\n",
			          self->swarmVarName, __func__, self->materialSwarm->name );

		Stg_Component_Build( self->swarmVar, data, False );
		Stg_Component_Initialise( self->swarmVar, data, False );
	}

	_ParticleFeVariable_Initialise( self, data );

	_SwarmVariable_Initialise( self->swarmVar, data );
}

void _SwarmVariableField_Execute( void* swarmVariableField, void* data ) {
	_ParticleFeVariable_Execute( swarmVariableField, data );
}

void _SwarmVariableField_Destroy( void* swarmVariableField, void* data ) {
	_ParticleFeVariable_Destroy( swarmVariableField, data );
}

void _SwarmVariableField_ValueAtParticle( void* swarmVariableField, 
					  IntegrationPointsSwarm* swarm, 
					  Element_LocalIndex lElement_I, 
					  IntegrationPoint* particle,
					  double* value ) 
{
	SwarmVariableField*	self            = (SwarmVariableField*)swarmVariableField;
	GlobalParticle*		matParticle;
	double			distance;
	Cell_Index		cell_I;
	Particle_InCellIndex	cParticle_I;
	Particle_Index		lParticle_I;

	cell_I = CellLayout_MapElementIdToCellId( swarm->cellLayout, lElement_I );
	cParticle_I = Swarm_FindClosestParticleInCell( swarm,
		       				       cell_I,
						       Mesh_GetDimSize( self->dofLayout->mesh ),
						       particle->xi,
						       &distance );
	// this function doesn't seem to be doing its joob properly!
	//lParticle_I = IntegrationPointMapper_GetMaterialIndexAt( swarm->mapper, swarm->cellParticleTbl[cell_I][cParticle_I] );
	
	/* assume that the material and intergation points swarms map 1:1 */
	lParticle_I = swarm->cellParticleTbl[cell_I][cParticle_I];
	
	
	SwarmVariable_ValueAt( self->swarmVar, lParticle_I, value ); /* does the copy inside this func. dave, 18.09.07 */
}

void _SwarmVariableField_GetValueAtNode( void* swarmVariableField, Node_DomainIndex dNode_I, double* value ) {
	FeVariable_GetValueAtNode( swarmVariableField, dNode_I, value );
}

/* implement these two functions later... */
/*
double _SwarmVariableField_GetMinGlobalMagnitude( void* swarmVariableField ) {
	return _SwarmVariable_GetMinGlobalMagnitude( swarmVariableField );
}
double _SwarmVariableField_GetMaxGlobalMagnitude( void* swarmVariableField ) {
	return _SwarmVariable_GetMaxGlobalMagnitude( swarmVariableField );
}
*/



