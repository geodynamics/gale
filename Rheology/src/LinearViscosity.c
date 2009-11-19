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
** $Id: LinearViscosity.c 743 2008-06-23 01:49:43Z JulianGiordani $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "RheologyClass.h"
#include "ConstitutiveMatrix.h"
#include "LinearViscosity.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type LinearViscosity_Type = "LinearViscosity";

/* Public Constructor */
LinearViscosity* LinearViscosity_New(
      Name              name,
      AbstractContext*  context,
      double            C,
      double            X,
      double            Y,
      double            Z,
      double            XY,
      double            XZ,
      double            YZ,
      double            XYZ )
{
   LinearViscosity* self = (LinearViscosity*) _LinearViscosity_DefaultNew( name );

   _Rheology_Init( self, context );
   _LinearViscosity_Init( self,
	     C,
	     X,
	     Y,
	     Z,
	     XY,
	     XZ,
	     YZ,
	     XYZ );

   self->isConstructed = True;
   return self;
}

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
LinearViscosity* _LinearViscosity_New( 
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
		Rheology_ModifyConstitutiveMatrixFunction*         _modifyConstitutiveMatrix,
		Name                                               name ) 
{
	LinearViscosity*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(LinearViscosity) );
	self = (LinearViscosity*) _Rheology_New( 
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
			_modifyConstitutiveMatrix,
			name );
	
	return self;
}

void _LinearViscosity_Init( 
	LinearViscosity* self,
	double           C,
	double           X,
	double           Y,
	double           Z,
	double           XY,
	double           XZ,
	double           YZ,
	double           XYZ )
{
	self->C          = C;
	self->X          = X;
	self->Y          = Y;
	self->Z          = Z;
	self->XY         = XY;
	self->XZ         = XZ;
	self->YZ         = YZ;
	self->XYZ        = XYZ;
}

void* _LinearViscosity_DefaultNew( Name name ) {
	return (void*) _LinearViscosity_New(
		sizeof(LinearViscosity),
		LinearViscosity_Type,
		_Rheology_Delete,
		_Rheology_Print,
		_Rheology_Copy,
		_LinearViscosity_DefaultNew,
		_LinearViscosity_AssignFromXML,
		_Rheology_Build,
		_Rheology_Initialise,
		_Rheology_Execute,
		_Rheology_Destroy,
		_LinearViscosity_ModifyConstitutiveMatrix,
		name );
}

void _LinearViscosity_AssignFromXML( void* rheology, Stg_ComponentFactory* cf, void* data ){
	LinearViscosity*  self                   = (LinearViscosity*)rheology;

	/* Construct Parent */
	_Rheology_AssignFromXML( self, cf, data );
	
	_LinearViscosity_Init( 
			self, 
			Stg_ComponentFactory_GetDouble( cf, self->name, "C", 1.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "X", 1.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "Y", 1.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "Z", 1.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "XY", 1.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "XZ", 1.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "YZ", 1.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "XYZ", 1.0 ) );
}

void _LinearViscosity_ModifyConstitutiveMatrix( 
		void*                                              rheology, 
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               swarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi )
{
	LinearViscosity*  self              = (LinearViscosity*) rheology;
	double            viscosity;
	Coord             coord;
	double            x, y, z;
	double            C                 = self->C;
	double            X                 = self->X;
	double            Y                 = self->Y;
	double            Z                 = self->Z;
	double            XY                = self->XY;
	double            XZ                = self->XZ;
	double            YZ                = self->YZ;
	double            XYZ               = self->XYZ;
	FeMesh*           mesh              = ConstitutiveMatrix_GetMesh( constitutiveMatrix );
	

	/* This rheology assumes particle is an integration points thats can be mapped to a particle
	 * that has no meaningful coord. Best to re-calc the global from local */
	FeMesh_CoordLocalToGlobal( mesh, lElement_I, xi, coord );
	x = coord[ I_AXIS ];
	y = coord[ J_AXIS ];
	z = coord[ K_AXIS ];

	/* Calculate New Viscosity */
	viscosity = X * x + Y * y + XY * x * y + C;
	if ( swarm->dim == 3 )
		viscosity += Z * z + XZ * x * z + YZ * y * z + XYZ * x * y * z;
		
	ConstitutiveMatrix_SetIsotropicViscosity( constitutiveMatrix, viscosity );
}
