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
** $Id: MaterialFeVariable.c 462 2007-05-16 01:13:21Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/Voronoi/Voronoi.h>
#include <PICellerator/PopulationControl/PopulationControl.h>
#include <PICellerator/Weights/Weights.h>

#include "MaterialPoints.h"

#include <assert.h>

const Type MaterialFeVariable_Type = "MaterialFeVariable";

MaterialFeVariable* _MaterialFeVariable_New(
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
		FieldVariable_GetCoordFunction*                   _getMinAndMaxLocalCoords,
		FieldVariable_GetCoordFunction*                   _getMinAndMaxGlobalCoords,		
		FeVariable_InterpolateWithinElementFunction*      _interpolateWithinElement,	
		FeVariable_GetValueAtNodeFunction*                _getValueAtNode,
		ParticleFeVariable_ValueAtParticleFunction*       _valueAtParticle,
		Name                                              name )
{
	MaterialFeVariable*		self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(MaterialFeVariable) );
	self = (MaterialFeVariable*)
		_ParticleFeVariable_New(
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
	
	return self;
}

void _MaterialFeVariable_Init( MaterialFeVariable* self, Material* material ) {
	IntegrationPointsSwarm* swarm;

	/* Assign Pointers */
	swarm = Stg_CheckType( self->assemblyTerm->integrationSwarm, IntegrationPointsSwarm );
	self->picIntegrationPoints = swarm;
	self->material       = material;
}

/* --- Virtual Function Implementations --- */
void _MaterialFeVariable_Delete( void* materialFeVariable ) {
	MaterialFeVariable* self = (MaterialFeVariable*) materialFeVariable;

	_ParticleFeVariable_Delete( self );
}

void _MaterialFeVariable_Print( void* materialFeVariable, Stream* stream ) {
	MaterialFeVariable* self = (MaterialFeVariable*) materialFeVariable;
	
	/* General info */
	Journal_Printf( stream, "MaterialFeVariable (ptr): %p\n", self );
	
	/* Print parent */
	_ParticleFeVariable_Print( self, stream );
	
	/* MaterialFeVariable info */
	Journal_PrintPointer( stream, self->material );
}


void* _MaterialFeVariable_Copy( void* feVariable, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	MaterialFeVariable*	self = (MaterialFeVariable*)feVariable;
	MaterialFeVariable*	newMaterialFeVariable;
	
	newMaterialFeVariable = (MaterialFeVariable*) _ParticleFeVariable_Copy( feVariable, dest, deep, nameExt, ptrMap );

	newMaterialFeVariable->picIntegrationPoints = self->picIntegrationPoints;
	newMaterialFeVariable->material = self->material;
	
	return (void*)newMaterialFeVariable;
}

void* _MaterialFeVariable_DefaultNew( Name name ) {
	return (void*) _MaterialFeVariable_New(
		sizeof(MaterialFeVariable),
		MaterialFeVariable_Type,
		_MaterialFeVariable_Delete,
		_MaterialFeVariable_Print,
		_MaterialFeVariable_Copy,
		_MaterialFeVariable_DefaultNew,
		_MaterialFeVariable_Construct,
		_MaterialFeVariable_Build, 
		_MaterialFeVariable_Initialise,
		_MaterialFeVariable_Execute,
		_MaterialFeVariable_Destroy,
		_FeVariable_InterpolateValueAt,
		_FeVariable_GetMinGlobalFieldMagnitude,
		_FeVariable_GetMaxGlobalFieldMagnitude,
		_FeVariable_GetMinAndMaxLocalCoords,
		_FeVariable_GetMinAndMaxGlobalCoords,
		_FeVariable_InterpolateNodeValuesToElLocalCoord,
		_FeVariable_GetValueAtNode,
		_MaterialFeVariable_ValueAtParticle,
		name );
}

void _MaterialFeVariable_Construct( void* materialFeVariable, Stg_ComponentFactory* cf, void* data ){
	MaterialFeVariable*   self              = (MaterialFeVariable*) materialFeVariable;
	Material*             material;
	
	material = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Material", Material, True, data );

	/* Construct Parent */
	_ParticleFeVariable_Construct( self, cf, data );

	_FieldVariable_Construct( self, cf, data );
	_MaterialFeVariable_Init( self, material );
}

void _MaterialFeVariable_Build( void* materialFeVariable, void* data ) {
	MaterialFeVariable* self = (MaterialFeVariable*) materialFeVariable;
	IntegrationPointsSwarm* swarm;
	Name tmpName;
	Variable_Register* variable_Register = NULL;

	Build( self->feMesh, data, False );

	/* Create Dof Layout */
	swarm = self->picIntegrationPoints;
	if ( swarm->swarmVariable_Register )
		variable_Register = swarm->swarmVariable_Register->variable_Register;

	tmpName = Stg_Object_AppendSuffix( self, "DataVariable" );
	self->dataVariable = Variable_NewScalar( 
			tmpName,
			Variable_DataType_Double, 
			&self->feMesh->topo->remotes[MT_VERTEX]->nDomains, 
			(void**)&self->data, 
			variable_Register );
	Memory_Free( tmpName );
	self->fieldComponentCount = 1;
	
	tmpName = Stg_Object_AppendSuffix( self, "DofLayout" );
	self->dofLayout = DofLayout_New( tmpName, variable_Register, self->feMesh->topo->remotes[MT_VERTEX]->nDomains, NULL );
	DofLayout_AddAllFromVariableArray( self->dofLayout, 1, &self->dataVariable );
	Memory_Free( tmpName );
	self->eqNum->dofLayout = self->dofLayout;
	
	_ParticleFeVariable_Build( self, data );
}
void _MaterialFeVariable_Initialise( void* materialFeVariable, void* data ) {
	MaterialFeVariable* self = (MaterialFeVariable*) materialFeVariable;

	_ParticleFeVariable_Initialise( self, data );
}
void _MaterialFeVariable_Execute( void* materialFeVariable, void* data ) {
	MaterialFeVariable* self = (MaterialFeVariable*) materialFeVariable;

	_ParticleFeVariable_Execute( self, data );
}
void _MaterialFeVariable_Destroy( void* materialFeVariable, void* data ) {
	MaterialFeVariable* self = (MaterialFeVariable*) materialFeVariable;

	_ParticleFeVariable_Destroy( self, data );
}
void _MaterialFeVariable_ValueAtParticle( 
		void*                   materialFeVariable,
		IntegrationPointsSwarm* swarm,
		Element_LocalIndex      lElement_I,
		void*                   particle,
		double*                 particleValue )
{
	MaterialFeVariable* self = (MaterialFeVariable*) materialFeVariable;
	*particleValue = (double) ( self->material->index == IntegrationPointsSwarm_GetMaterialIndexOn( swarm, particle ) );
}

