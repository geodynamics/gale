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
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "RheologyClass.h"
#include "ConstitutiveMatrix.h"
#include "OrthotropicAligned.h"

#include <assert.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type OrthotropicAligned_Type = "OrthotropicAligned";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
OrthotropicAligned* _OrthotropicAligned_New( 
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
	OrthotropicAligned*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(OrthotropicAligned) );
	self = (OrthotropicAligned*) _Rheology_New( 
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

void _OrthotropicAligned_Init( OrthotropicAligned* self,
                             /*  MaterialPointsSwarm* materialPointsSwarm, */
			       double viscosity1, double viscosity2, 
			       double viscosity3, double viscosity4, 
			       double viscosity5, double viscosity6 ) {
//	self->director         = director;
//	self->viscosityRatio  = viscosityRatio;
/*	OrthotropicAligned_ParticleExt*   particleExt;
	StandardParticle                  materialPoint;
	
	self->particleExtHandle       = ExtensionManager_Add( materialPointsSwarm->particleExtensionMgr,
	                                                	OrthotropicAligned_Type, sizeof(OrthotropicAligned_ParticleExt) ); */
        self->viscosity1  = viscosity1;
	self->viscosity2  = viscosity2;
	self->viscosity3  = viscosity3;
	self->viscosity4  = viscosity4;
	self->viscosity5  = viscosity5;  
	self->viscosity6  = viscosity6;
}

void* _OrthotropicAligned_DefaultNew( Name name ) {
	return (void*) _OrthotropicAligned_New(
		sizeof(OrthotropicAligned),
		OrthotropicAligned_Type,
		_Rheology_Delete,
		_Rheology_Print,
		_Rheology_Copy,
		_OrthotropicAligned_DefaultNew,
		_OrthotropicAligned_Construct,
		_Rheology_Build,
		_Rheology_Initialise,
		_Rheology_Execute,
		_Rheology_Destroy,
		_OrthotropicAligned_ModifyConstitutiveMatrix,
		name );
}

void _OrthotropicAligned_Construct( void* rheology, Stg_ComponentFactory* cf, void* data ){
	OrthotropicAligned*     self = (OrthotropicAligned*)rheology;
//	Director*        director;
	double viscosity1;
	double viscosity2;
	double viscosity3;
	double viscosity4;
	double viscosity5;
	double viscosity6;

	/* Construct Parent */
	_Rheology_Construct( self, cf, data );
	
//	director =  Stg_ComponentFactory_ConstructByKey(  cf,  self->name,  "Director", Director,  True  ) ;
	viscosity1 = Stg_ComponentFactory_GetDouble( cf, self->name, "viscosity1",  True );
	viscosity2 = Stg_ComponentFactory_GetDouble( cf, self->name, "viscosity2",  True );
	viscosity3 = Stg_ComponentFactory_GetDouble( cf, self->name, "viscosity3",  True );
	viscosity4 = Stg_ComponentFactory_GetDouble( cf, self->name, "viscosity4",  True );
	viscosity5 = Stg_ComponentFactory_GetDouble( cf, self->name, "viscosity5",  True );
	viscosity6 = Stg_ComponentFactory_GetDouble( cf, self->name, "viscosity6",  True );

	_OrthotropicAligned_Init( 
			self, 
			viscosity1,
			viscosity2,
			viscosity3,
			viscosity4,
			viscosity5,
			viscosity6 );
}

/* how does this get called eventually? */
void _OrthotropicAligned_ModifyConstitutiveMatrix( 
		void*                                              rheology, 
		ConstitutiveMatrix*                                constitutiveMatrix, // constitutive matrix
		MaterialPointsSwarm*                               swarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi )
{
	OrthotropicAligned*	                self = (OrthotropicAligned*) rheology;
	Dimension_Index                   dim  = swarm->dim;
//	double                          isotropicViscosity = ConstitutiveMatrix_GetIsotropicViscosity( constitutiveMatrix );
//	double                          deltaViscosity;
//	XYZ                             normal;
	int i,j;
	double**   D  = constitutiveMatrix->matrixData;
//	static int flag = 0;
//	deltaViscosity = isotropicViscosity * (1.0 - self->viscosityRatio);
//	Director_GetNormal( self->director, materialPoint, normal );

//	ConstitutiveMatrix_SetSecondViscosity( constitutiveMatrix, deltaViscosity, normal );
	
//	if(!flag){/* if not visited modify matrix else no need to update */
        /* Snark dies if I only allow this to be called once.. */
	/* ahh need to allow it to be called once for every particle */
	   for(i=0;i<dim*(dim+1)/2;i++){
	      for(j=0;j<dim*(dim+1)/2;j++){
		 D[i][j] = 0.0;
	      }
	   }
	   
	   D[0][0] = self->viscosity1;
	   D[1][1] = self->viscosity2;
	   D[2][2] = self->viscosity3;
	   if(dim == 3){
	      D[3][3] = self->viscosity4;
	      D[4][4] = self->viscosity5;
	      D[5][5] = self->viscosity6;
	   }
	   constitutiveMatrix->isDiagonal = True;	   
//	   printf("In %s OK\n\n",__func__);
//	   flag = 1;
//	}
/* 	for(i=0;i<dim*(dim+1)/2;i++){ */
/* 	   for(j=0;j<dim*(dim+1)/2;j++){ */
/* 	      printf("Matrix Data = %g [%d %d]\n",constitutiveMatrix->matrixData[i][j],i,j); */
/* 	   } */
/* 	} */
}

#if 0
void _OrthotropicAligned_UpdateDrawParameters( void* rheology ) 
{
	OrthotropicAligned*                   self               = (OrthotropicAligned*) rheology;
	Particle_Index                   lParticle_I;
	Particle_Index                   particleLocalCount;
	StandardParticle*                	materialPoint;
	Dimension_Index                 dim                   = self->materialPointsSwarm->dim;
	OrthotropicAligned_ParticleExt*            particleExt;   /* this new type is defined in OrthotropicAligned.h */

	/* do stuff */
}
#endif
