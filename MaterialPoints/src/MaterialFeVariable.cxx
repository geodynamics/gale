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
** $Id: MaterialFeVariable.c 518 2007-10-11 08:07:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PopulationControl/PopulationControl.h>
#include <PICellerator/Weights/Weights.h>

#include "MaterialPoints.h"

#include <assert.h>

const Type MaterialFeVariable_Type = "MaterialFeVariable";

MaterialFeVariable* _MaterialFeVariable_New(  MATERIALFEVARIABLE_DEFARGS  ) {
	MaterialFeVariable* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(MaterialFeVariable) );
	self = (MaterialFeVariable*) _ParticleFeVariable_New(  PARTICLEFEVARIABLE_PASSARGS  );
	
	return self;
}

void _MaterialFeVariable_Init( MaterialFeVariable* self, Material* material ) {
	IntegrationPointsSwarm* swarm;

	/* Assign Pointers */
	swarm = Stg_CheckType( self->assemblyTerm->integrationSwarm, IntegrationPointsSwarm );
	self->picIntegrationPoints = swarm;
	self->material = material;
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

void* _MaterialFeVariable_Copy( const void* feVariable, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	MaterialFeVariable*	self = (MaterialFeVariable*)feVariable;
	MaterialFeVariable*	newMaterialFeVariable;
	
	newMaterialFeVariable = (MaterialFeVariable*) _ParticleFeVariable_Copy( feVariable, dest, deep, nameExt, ptrMap );

	newMaterialFeVariable->picIntegrationPoints = self->picIntegrationPoints;
	newMaterialFeVariable->material = self->material;
	
	return (void*)newMaterialFeVariable;
}

void* _MaterialFeVariable_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                       _sizeOfSelf = sizeof(MaterialFeVariable);
	Type                                                               type = MaterialFeVariable_Type;
	Stg_Class_DeleteFunction*                                       _delete = _MaterialFeVariable_Delete;
	Stg_Class_PrintFunction*                                         _print = _MaterialFeVariable_Print;
	Stg_Class_CopyFunction*                                           _copy = _MaterialFeVariable_Copy;
	Stg_Component_DefaultConstructorFunction*           _defaultConstructor = _MaterialFeVariable_DefaultNew;
	Stg_Component_ConstructFunction*                             _construct = _MaterialFeVariable_AssignFromXML;
	Stg_Component_BuildFunction*                                     _build = _MaterialFeVariable_Build;
	Stg_Component_InitialiseFunction*                           _initialise = _MaterialFeVariable_Initialise;
	Stg_Component_ExecuteFunction*                                 _execute = _MaterialFeVariable_Execute;
	Stg_Component_DestroyFunction*                                 _destroy = _MaterialFeVariable_Destroy;
	FieldVariable_InterpolateValueAtFunction*           _interpolateValueAt = _FeVariable_InterpolateValueAt;
	FieldVariable_GetCoordFunction*                _getMinAndMaxLocalCoords = _FeVariable_GetMinAndMaxLocalCoords;
	FieldVariable_GetCoordFunction*               _getMinAndMaxGlobalCoords = _FeVariable_GetMinAndMaxGlobalCoords;
	FeVariable_InterpolateWithinElementFunction*  _interpolateWithinElement = _FeVariable_InterpolateNodeValuesToElLocalCoord;
	FeVariable_GetValueAtNodeFunction*                      _getValueAtNode = _FeVariable_GetValueAtNode;
	ParticleFeVariable_ValueAtParticleFunction*            _valueAtParticle = _MaterialFeVariable_ValueAtParticle;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType                             nameAllocationType = (AllocationType)ZERO;
	FieldVariable_GetValueFunction*   _getMinGlobalFieldMagnitude = NULL;
	FieldVariable_GetValueFunction*   _getMaxGlobalFieldMagnitude = NULL;
	FeVariable_SyncShadowValuesFunc*            _syncShadowValues = NULL;

	return (void*) _MaterialFeVariable_New(  MATERIALFEVARIABLE_PASSARGS  );
}

void _MaterialFeVariable_AssignFromXML( void* materialFeVariable, Stg_ComponentFactory* cf, void* data ){
	MaterialFeVariable*	self = (MaterialFeVariable*) materialFeVariable;
	Material*				material;
	
	material = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"Material", Material, True, data  );

	/* Construct Parent */
	_ParticleFeVariable_AssignFromXML( self, cf, data );

	_MaterialFeVariable_Init( self, material );
}

void _MaterialFeVariable_Build( void* materialFeVariable, void* data ) {
	MaterialFeVariable* self = (MaterialFeVariable*) materialFeVariable;
	IntegrationPointsSwarm* swarm;
	char *tmpName;
	Variable_Register* variable_Register = NULL;

	Stg_Component_Build( self->feMesh, data, False );

	/* Create Dof Layout */
	swarm = self->picIntegrationPoints;
	if ( swarm->swarmVariable_Register )
		variable_Register = swarm->swarmVariable_Register->variable_Register;

	tmpName = Stg_Object_AppendSuffix( self, "DataVariable"  );
	assert( Class_IsSuper( self->feMesh->topo, IGraph ) );
	self->dataVariable = Variable_NewScalar( 
		tmpName,
		(AbstractContext*)self->context,
		Variable_DataType_Double, 
		(Index*)(&((IGraph*)self->feMesh->topo)->remotes[MT_VERTEX]->nDomains),
		NULL,
		(void**)&self->data, 
		variable_Register );
	Memory_Free( tmpName );
	self->fieldComponentCount = 1;
	
	tmpName = Stg_Object_AppendSuffix( self, (Name)"DofLayout"  );
	self->dofLayout = DofLayout_New( tmpName, self->context, variable_Register, ((IGraph*)self->feMesh->topo)->remotes[MT_VERTEX]->nDomains, NULL );
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
	*particleValue = (double) ( self->material->index == IntegrationPointsSwarm_GetMaterialIndexOn( swarm, (IntegrationPoint*)particle ) );
}



