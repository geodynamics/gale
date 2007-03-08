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
** $Id: MultiRheologyMaterial.c 446 2007-03-04 09:55:46Z PatrickSunter $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "RheologyMaterial.h"
#include "Rheology_Register.h"
#include "MultiRheologyMaterial.h"
#include "ConstitutiveMatrix.h"
#include "Compressible.h"
#include "RheologyClass.h"

#include <assert.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type MultiRheologyMaterial_Type = "MultiRheologyMaterial";


MultiRheologyMaterial* MultiRheologyMaterial_New( 
		Name                             name,
		Stg_Shape*                       shape,
		Dictionary*                      materialDictionary,
		Materials_Register*              materialRegister,
		Rheology**                       rheologyList,
		Rheology_Index                   rheologyCount,
		Compressible*                    compressible,
		Rheology***                      rheologyListList,
		Rheology_Index*                  rheologyCountList, 
		Index                            rheologyListCount ) 
{
	return (void*) _MultiRheologyMaterial_New(
		sizeof(MultiRheologyMaterial),
		MultiRheologyMaterial_Type,
		_MultiRheologyMaterial_Delete,
		_RheologyMaterial_Print,
		_RheologyMaterial_Copy,
		_MultiRheologyMaterial_DefaultNew,
		_MultiRheologyMaterial_Construct,
		_RheologyMaterial_Build,
		_RheologyMaterial_Initialise,
		_RheologyMaterial_Execute,
		_RheologyMaterial_Destroy,
		_MultiRheologyMaterial_RunRheologies,
		name,
		shape,
		materialDictionary,
		materialRegister,
		rheologyList,
		rheologyCount,
		compressible,
		rheologyListList,
		rheologyCountList,
		rheologyListCount );
}


void* _MultiRheologyMaterial_DefaultNew( Name name ) {
	return (void*) _MultiRheologyMaterial_New(
		sizeof(MultiRheologyMaterial),
		MultiRheologyMaterial_Type,
		_MultiRheologyMaterial_Delete,
		_RheologyMaterial_Print,
		_RheologyMaterial_Copy,
		_MultiRheologyMaterial_DefaultNew,
		_MultiRheologyMaterial_Construct,
		_RheologyMaterial_Build,
		_RheologyMaterial_Initialise,
		_RheologyMaterial_Execute,
		_RheologyMaterial_Destroy,
		_MultiRheologyMaterial_RunRheologies,
		name,
		NULL,
		NULL,
		NULL,
		NULL,
		0,
		NULL,
		NULL,
		NULL,
		0
		);
}


/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
MultiRheologyMaterial* _MultiRheologyMaterial_New( 
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
		Compressible*                                      compressible,
		Rheology***                                        rheologyListList,
		Rheology_Index*                                    rheologyCountList, 
		Index                                              rheologyListCount
		) 
{
	MultiRheologyMaterial*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up
	 *  the hierarchy tree. At the beginning of the tree it will allocate memory of the size of
	 *  object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(MultiRheologyMaterial) );
	self = (MultiRheologyMaterial*) _RheologyMaterial_New( 
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
			_runRheologies,
			name,
			shape,
			materialDictionary,
			materialRegister,
			rheologyList,
			rheologyCount,
			compressible );

	_MultiRheologyMaterial_Init(
		self,
		rheologyListList,
		rheologyCountList,
		rheologyListCount );
	
	return self;
}


void _MultiRheologyMaterial_Construct( void* material, Stg_ComponentFactory* cf, void* data ){
	MultiRheologyMaterial*  self = (MultiRheologyMaterial*)material;
	Rheology***             rheologyListList;
	Rheology_Index*         rheologyCountList;
	Index                   rheologyListCount;
	Index                   rheologyList_I;
	Dictionary*             currDictionary     = Dictionary_GetDictionary( cf->componentDict, self->name );
	Name                    rheologyName;
	Dictionary_Entry_Value* multiRheologyList;
	Dictionary_Entry_Value* rheologyList;
	Dictionary_Entry_Value* rheologyEntry;
	Rheology_Index          rheology_I;
	Stream*                 stream             = cf->infoStream;

	/* Construct Parent */
	_Material_Construct( self, cf, data );
	_RheologyMaterial_Init( 
		self, 
		NULL, 
		0, 
		Stg_ComponentFactory_ConstructByKey( cf, self->name, "Compressible", Compressible, False, data ) );

	multiRheologyList = Dictionary_Get( currDictionary, "MultiRheologies" );
	assert( multiRheologyList );
	rheologyListCount = Dictionary_Entry_Value_GetCount( multiRheologyList );

	rheologyCountList = Memory_Alloc_Array( Rheology_Index, rheologyListCount, "rheologyCountList" );
	rheologyListList = Memory_Alloc_Array( Rheology**, rheologyListCount, "rheologyListList" );
			
	Journal_PrintfL( stream, 2, "%s '%s' has %u rheological components.\n", self->type, self->name, rheologyListCount );
	Stream_Indent( stream );

	for ( rheologyList_I = 0 ; rheologyList_I < rheologyListCount ; rheologyList_I++ ){
		rheologyList = Dictionary_Entry_Value_GetElement( multiRheologyList, rheologyList_I );
		rheologyCountList[ rheologyList_I ] = Dictionary_Entry_Value_GetCount( rheologyList );
		rheologyListList[ rheologyList_I ] = Memory_Alloc_Array( Rheology*, rheologyCountList[ rheologyList_I ], "rheologyList" );
	
		Journal_PrintfL( stream, 2, "%s \"%s's\" rheological component '%u' has %u rheologies.\n", 
				self->type, self->name, rheologyList_I, rheologyCountList[ rheologyList_I ] );

		Stream_Indent( stream );
		for ( rheology_I = 0 ; rheology_I < rheologyCountList[ rheologyList_I ] ; rheology_I++ ) {
			rheologyEntry = Dictionary_Entry_Value_GetElement( rheologyList, rheology_I );

			rheologyName = Dictionary_Entry_Value_AsString( rheologyEntry );
		
			Journal_PrintfL( stream, 2, "Rheological component '%u' for %s '%s' needs rheology: '%s'.\n", 
				rheologyList_I, self->type, self->name, rheologyName );

			rheologyListList[ rheologyList_I ][ rheology_I ] = 
	Stg_ComponentFactory_ConstructByName( cf, rheologyName, Rheology, True, data ); 
		}
		Stream_UnIndent( stream );
	}
	Stream_UnIndent( stream );
	
	_MultiRheologyMaterial_Init( 
			self,
			rheologyListList,
			rheologyCountList,
			rheologyListCount );

	/* Clean up */
	for ( rheologyList_I = 0 ; rheologyList_I < rheologyListCount ; rheologyList_I++ ){
		Memory_Free( rheologyListList[ rheologyList_I ] );
	}
	Memory_Free( rheologyListList );
	Memory_Free( rheologyCountList );
}


void _MultiRheologyMaterial_Init( 
		MultiRheologyMaterial*  self, 
		Rheology***             rheologyListList,
		Rheology_Index*         rheologyCountList, 
		Index                   rheologyListCount ) 
{
	Index          rheology_Register_I;
	Rheology_Index rheology_I;
	Rheology*      rheology;

	self->rheology_RegisterCount = rheologyListCount;
	self->rheology_RegisterList = Memory_Alloc_Array( Rheology_Register*, rheologyListCount, "rheology_RegisterList" );
	for ( rheology_Register_I = 0 ; rheology_Register_I < rheologyListCount ; rheology_Register_I++ ) {
		self->rheology_RegisterList[ rheology_Register_I ] = Rheology_Register_New();

		for ( rheology_I = 0 ; rheology_I < rheologyCountList[ rheology_Register_I ] ; rheology_I++ ) {
			rheology = rheologyListList[ rheology_Register_I ][ rheology_I ];
			Rheology_Register_Add( self->rheology_Register, rheology );
			Rheology_Register_Add( self->rheology_RegisterList[ rheology_Register_I ], rheology );
		}
	}
}


void _MultiRheologyMaterial_Delete( void* material ){
	MultiRheologyMaterial*  self = (MultiRheologyMaterial*)material;
	Index                   rheology_Register_I;

	for ( rheology_Register_I = 0 ; rheology_Register_I < self->rheology_RegisterCount ; rheology_Register_I++ ) {
		Stg_Class_Delete( self->rheology_RegisterList[ rheology_Register_I ] );
	}
	Memory_Free( self->rheology_RegisterList );

	_RheologyMaterial_Delete( self );
}


void _MultiRheologyMaterial_RunRheologies( 	
		void*                                              material,
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               swarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi )
{
	MultiRheologyMaterial*   self                   = (MultiRheologyMaterial*) material;
	Index                    rheology_RegisterCount = self->rheology_RegisterCount;
	Index                    rheology_Register_I;
	double                   oneOnViscosity         = 0.0;

	for ( rheology_Register_I = 0 ; rheology_Register_I < rheology_RegisterCount ; rheology_Register_I++ ) {
		/* Calculate Isotropic Viscosity for this rheology register */
		Rheology_Register_RunRheologies( 
				self->rheology_RegisterList[ rheology_Register_I ],
				constitutiveMatrix,
				swarm, 
				lElement_I,
				materialPoint,
				xi );

		oneOnViscosity += 1.0/ConstitutiveMatrix_GetIsotropicViscosity( constitutiveMatrix );
		
		/* Initialise for next run through */
		ConstitutiveMatrix_ZeroMatrix( constitutiveMatrix );
	}

	ConstitutiveMatrix_SetIsotropicViscosity( constitutiveMatrix, 1.0/oneOnViscosity );
}
