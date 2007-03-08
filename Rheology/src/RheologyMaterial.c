/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2005, Monash Cluster Computing 
** All rights reserved.
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
** 		* Redistributions of source code must retain the above copyright notice, 
** 			this list of conditions and the following disclaimer.
** 		* Redistributions in binary form must reproduce the above copyright 
**			notice, this list of conditions and the following disclaimer in the 
**			documentation and/or other materials provided with the distribution.
** 		* Neither the name of the Monash University nor the names of its contributors 
**			may be used to endorse or promote products derived from this software 
**			without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
** THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
** PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS 
** BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
** CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
** SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
** HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
** LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
** OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
**
** Contact:
*%		Louis Moresi - Louis.Moresi@sci.monash.edu.au
*%
** Contributors:
*+		Robert Turnbull
*+		Vincent Lemiale
*+		Louis Moresi
*+		David May
*+		David Stegman
*+		Mirko Velic
*+		Patrick Sunter
*+		Julian Giordani
*+
** $Id: RheologyMaterial.c 446 2007-03-04 09:55:46Z PatrickSunter $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "RheologyMaterial.h"
#include "Rheology_Register.h"
#include "RheologyClass.h"
#include "Compressible.h"

#include <assert.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type RheologyMaterial_Type = "RheologyMaterial";

RheologyMaterial* RheologyMaterial_New( 
		Name                             name,
		Stg_Shape*                       shape,
		Dictionary*                      materialDictionary,
		Materials_Register*              materialRegister,
		Rheology**                       rheologyList,
		Rheology_Index                   rheologyCount,
		Compressible*                    compressible )
{
	return (void*) _RheologyMaterial_New( 
		sizeof(RheologyMaterial), 
		RheologyMaterial_Type, 
		_RheologyMaterial_Delete,
		_RheologyMaterial_Print, 
		_RheologyMaterial_Copy,
		_RheologyMaterial_DefaultNew,
		_RheologyMaterial_Construct,
		_RheologyMaterial_Build, 
		_RheologyMaterial_Initialise,
		_RheologyMaterial_Execute,
		_RheologyMaterial_Destroy,
		_RheologyMaterial_RunRheologies,
		name,
		shape,
		materialDictionary,
		materialRegister,
		rheologyList,
		rheologyCount,
		compressible );
}


void* _RheologyMaterial_DefaultNew( Name name ) {
	return (void*) _RheologyMaterial_New( 
		sizeof(RheologyMaterial), 
		RheologyMaterial_Type, 
		_RheologyMaterial_Delete,
		_RheologyMaterial_Print, 
		_RheologyMaterial_Copy,
		_RheologyMaterial_DefaultNew,
		_RheologyMaterial_Construct,
		_RheologyMaterial_Build, 
		_RheologyMaterial_Initialise,
		_RheologyMaterial_Execute,
		_RheologyMaterial_Destroy,
		_RheologyMaterial_RunRheologies,
		name,
		NULL,
		NULL,
		NULL,
		NULL,
		0,
		NULL
		);
}


/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
RheologyMaterial* _RheologyMaterial_New( 
		SizeT                                              sizeOfSelf,
		Type                                               type,
		Stg_Class_DeleteFunction*                          _delete,
		Stg_Class_PrintFunction*                           _print,
		Stg_Class_CopyFunction*                            _copy, 
		Stg_Component_DefaultConstructorFunction*          _defaultConstructor,
		Stg_Component_ConstructFunction*                   _construct,
		Stg_Component_BuildFunction*                       _build,
		Stg_Component_InitialiseFunction*                  _initialise,
		Stg_Component_ExecuteFunction*                     _execute,
		Stg_Component_DestroyFunction*                     _destroy,
		RheologyMaterial_RunRheologiesFunction*            _runRheologies,
		Name                                               name,
		Stg_Shape*                                         shape,
		Dictionary*                                        materialDictionary,
		Materials_Register*                                materialRegister,
		Rheology**                                         rheologyList,
		Rheology_Index                                     rheologyCount,
		Compressible*                                      compressible ) 
{
	RheologyMaterial*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the
 	 *  hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and
 	 *  initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(RheologyMaterial) );
	self = (RheologyMaterial*) _Material_New( 
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
			name,
			shape,
			materialDictionary,
			materialRegister
			);
	
	/* Function pointers for this class that are not on the parent class should be set here */
	self->_runRheologies = _runRheologies;
	
	/* Now call the member-assigner, Init() */
	_RheologyMaterial_Init( self, rheologyList, rheologyCount, compressible );

	return self;
}


void _RheologyMaterial_Construct( void* rheologyMaterial, Stg_ComponentFactory* cf, void* data ){
	RheologyMaterial*           self                 = (RheologyMaterial*)rheologyMaterial;
	Rheology**                  rheologyList;
	Rheology_Index              rheologyCount;
	Compressible*               compressible;

	_Material_Construct( self, cf, data );

	rheologyList = Stg_ComponentFactory_ConstructByList( 
		cf, 
		self->name, 
		"Rheology", 
		Stg_ComponentFactory_Unlimited, 
		Rheology, 
		True, 
		&rheologyCount,
		data );

	compressible = Stg_ComponentFactory_ConstructByKey( 
		cf, 
		self->name, 
		"Compressible", 
		Compressible, 
		False,
		data ) ;

	_RheologyMaterial_Init( self, rheologyList, rheologyCount, compressible );

	Memory_Free( rheologyList );
}


void _RheologyMaterial_Init(
		void*                                              rheologyMaterial,
		Rheology**                                         rheologyList,
		Rheology_Index                                     rheologyCount,
		Compressible*                                      compressible )
{
	RheologyMaterial* self = (RheologyMaterial*)rheologyMaterial;
	Rheology_Index    rheology_I;

	self->compressible      = compressible;
	self->rheology_Register = Rheology_Register_New();

	/* Add rheologies */
	for ( rheology_I = 0 ; rheology_I < rheologyCount ; rheology_I++ ) 
		Rheology_Register_Add( self->rheology_Register, rheologyList[ rheology_I ] );

//	self->debug = Journal_Register( Debug_Type, self->type ); /* TODO make child of Underworld_Debug */
}


void _RheologyMaterial_Delete( void* rheologyMaterial ) {
	RheologyMaterial*					self = (RheologyMaterial*)rheologyMaterial;

	Stg_Class_Delete( self->rheology_Register );
	_Stg_Component_Delete( self );
}


void _RheologyMaterial_Print( void* rheologyMaterial, Stream* stream ) {}


void* _RheologyMaterial_Copy( void* rheologyMaterial, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	RheologyMaterial*	self = (RheologyMaterial*)rheologyMaterial;

	/* TODO */ abort();
	return (void*) self;
}


void _RheologyMaterial_Build( void* rheologyMaterial, void* data ) {
	RheologyMaterial*   self                 = (RheologyMaterial*)rheologyMaterial;

	_Material_Build( self, data );
}
void _RheologyMaterial_Initialise( void* rheologyMaterial, void* data ) {
	RheologyMaterial*   self                 = (RheologyMaterial*)rheologyMaterial;

	_Material_Initialise( self, data );
}
void _RheologyMaterial_Execute( void* rheologyMaterial, void* data ) {
	RheologyMaterial*   self                 = (RheologyMaterial*)rheologyMaterial;
	_Material_Execute( self, data );
}
void _RheologyMaterial_Destroy( void* rheologyMaterial, void* data ) {
	RheologyMaterial*   self                 = (RheologyMaterial*)rheologyMaterial;
	_Material_Destroy( self, data );

}

void RheologyMaterial_RunRheologies( 	
		void*                                              rheologyMaterial,
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               swarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi ) 
{
	RheologyMaterial*   self                 = (RheologyMaterial*)rheologyMaterial;

	self->_runRheologies( self, constitutiveMatrix, swarm, lElement_I, materialPoint, xi );
}

void _RheologyMaterial_RunRheologies( 	
		void*                                              rheologyMaterial,
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               swarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi )
{
	RheologyMaterial*   self                 = (RheologyMaterial*)rheologyMaterial;

	Rheology_Register_RunRheologies( 
			self->rheology_Register,
			constitutiveMatrix,
			swarm, 
			lElement_I,
			materialPoint,
			xi );
}

Bool RheologyMaterial_IsNonLinear( void* rheologyMaterial ) {
	RheologyMaterial*   self                 = (RheologyMaterial*)rheologyMaterial;
	Rheology_Index      rheology_I; 
	Rheology_Index      rheologyCount = Rheology_Register_GetCount( self->rheology_Register ); 
	Rheology*           rheology; 
	
	Journal_DFirewall( 
			rheologyCount > 0, 
			Journal_Register( Error_Type, self->type ), 
			"No rheologies registered on %s '%s'.\n", self->type, self->name ); 
	
	for( rheology_I = 0; rheology_I < rheologyCount ; rheology_I++ ) { 
		rheology = Rheology_Register_GetByIndex( self->rheology_Register, rheology_I ); 
		if ( rheology->nonLinear )
			return True;
	}

	return False;
}

Rheology* RheologyMaterial_GetRheologyByType( void* rheologyMaterial, Type type ) {
	RheologyMaterial*   self                 = (RheologyMaterial*)rheologyMaterial;
	Rheology_Index      rheology_I; 
	Rheology_Index      rheologyCount = Rheology_Register_GetCount( self->rheology_Register ); 
	Rheology*           rheology;

	for( rheology_I = 0; rheology_I < rheologyCount ; rheology_I++ ) { 
		rheology = Rheology_Register_GetByIndex( self->rheology_Register, rheology_I ); 

		if ( Stg_Class_IsInstance( rheology, type ) )
			return rheology;
	}

	return NULL;
}

Bool RheologyMaterial_HasRheology( void* rheologyMaterial, void* rheology ) {
	RheologyMaterial*   self                 = (RheologyMaterial*)rheologyMaterial;
	Rheology_Index      rheology_I; 
	Rheology_Index      rheologyCount = Rheology_Register_GetCount( self->rheology_Register ); 
	Rheology*           currRheology;
	
	for( rheology_I = 0; rheology_I < rheologyCount ; rheology_I++ ) { 
		currRheology = Rheology_Register_GetByIndex( self->rheology_Register, rheology_I ); 

		if ( rheology == currRheology )
			return True;
	}
	return False;
}
