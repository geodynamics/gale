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
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "RheologyClass.h"
#include "ConstitutiveMatrix.h"
#include "Orthotropic.h"

#include <assert.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type Orthotropic_Type = "Orthotropic";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
Orthotropic* _Orthotropic_New(  ORTHOTROPIC_DEFARGS  ) 
{
	Orthotropic*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(Orthotropic) );
	self = (Orthotropic*) _Rheology_New(  RHEOLOGY_PASSARGS  );
	
	return self;
}

void _Orthotropic_Init( Orthotropic* self,
			/*  MaterialPointsSwarm* materialPointsSwarm, */
			double C11, double C22, double C33, 
			double C44, double C55, double C66,
			double C12, double C13, double C23,
			double *n, double *m, double *q ) {
      int i;
      
      self->C11  = C11;
      self->C22  = C22;
      self->C33  = C33;
      self->C44  = C44;
      self->C55  = C55;  
      self->C66  = C66;
      self->C12  = C12;
      self->C13  = C13;
      self->C23  = C23;

      for(i=0;i<3;i++){
	    self->n[i] = n[i];
	    self->m[i] = m[i];
	    self->q[i] = q[i];
      }
}

void* _Orthotropic_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof(Orthotropic);
	Type                                                             type = Orthotropic_Type;
	Stg_Class_DeleteFunction*                                     _delete = _Rheology_Delete;
	Stg_Class_PrintFunction*                                       _print = _Rheology_Print;
	Stg_Class_CopyFunction*                                         _copy = _Rheology_Copy;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _Orthotropic_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _Orthotropic_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _Rheology_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _Rheology_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _Rheology_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _Rheology_Destroy;
	Rheology_ModifyConstitutiveMatrixFunction*  _modifyConstitutiveMatrix = _Orthotropic_ModifyConstitutiveMatrix;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _Orthotropic_New(  ORTHOTROPIC_PASSARGS  );
}

void _Orthotropic_AssignFromXML( void* rheology, Stg_ComponentFactory* cf, void* data ){

	Orthotropic*     self = (Orthotropic*)rheology;
	double C11;
	double C22;
	double C33;
	double C44;
	double C55;
	double C66;
	double C12;
	double C13;
	double C23;
	double n[3], m[3], q[3];

	/* Construct Parent */
	_Rheology_AssignFromXML( self, cf, data );
	
	/* get parameters */
	/*******************************************************************
         The input parameters Cij correspond to the components of the matrix
         aligned with the directions n,m,q:
         With the following correspondence between indices and vectors:
         1 = n,   2 = m, 3 = q.
         The indices on the vectors correspond to the coordinate system
         as follows:
         1 = x 2 = y 3 = z.
         *******************************************************************/
	C11 = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"C11", True  );
	C22 = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"C22", True  );
	C33 = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"C33", True  );
	C44 = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"C44", True  );
	C55 = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"C55", True  );
	C66 = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"C66", True  );
	C12 = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"C12", True  );
	C13 = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"C13", True  );
	C23 = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"C23", True  );
	n[0] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"n1", True  );
	n[1] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"n2", True  );
	n[2] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"n3", True  );
	m[0] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"m1", True  );
	m[1] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"m2", True  );
	m[2] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"m3", True  );
	q[0] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"q1", True  );
	q[1] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"q2", True  );
	q[2] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"q3", True  );

	_Orthotropic_Init( 
			self, 
			C11, C22, C33, C44, C55, C66,
			C12, C13, C23,
			n,m,q);
}

/* how does this get called eventually? */
void _Orthotropic_ModifyConstitutiveMatrix( 
		void*                                              rheology, 
		ConstitutiveMatrix*                                constitutiveMatrix, /* constitutive matrix */
		MaterialPointsSwarm*                               swarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi )
{
	Orthotropic*	                self = (Orthotropic*) rheology;
	Dimension_Index                   dim  = swarm->dim;

	double**   C  = constitutiveMatrix->matrixData;
	double n1,n2,n3;
	double m1,m2,m3;
	double q1,q2,q3;
	double C_h11, C_h22, C_h33, C_h44, C_h55, C_h66, C_h12, C_h13, C_h23;

	C_h11 = self->C11;
	C_h22 = self->C22;
	C_h33 = self->C33;
	C_h12 = self->C12;
	C_h13 = self->C13;
	C_h23 = self->C23;

	/* the vectors n, m and q must form an orthonormal triad in 3D */
	n1 = self->n[0];
	n2 = self->n[1];

	for(unsigned int i=0;i<dim*(dim+1)/2;i++){
	      for(unsigned int j=0;j<dim*(dim+1)/2;j++){
		    C[i][j] = 0.0;
	      }
	}
	

	if(2==dim) {/* vector m appears implicitly here. m=(-n2,n1) i.e. n dot m = 0 */
	      C[0][0] = n1*n1*n1*n1*C_h11+2*n1*n1*n2*n2*C_h12-2*n1*n1*n1*n2*C_h12+n2*n2*n2*n2*C_h22-4*n2*n2*n2*n1*C_h23-2*n1*n1*n1*n2*C_h13+4*n1*n1*n2*n2*C_h33;
	      C[0][1] = C[1][0] = n1*n1*n2*n2*C_h11+n2*n2*n2*n2*C_h12-2*n1*n2*n2*n2*C_h12+n1*n1*n1*n1*C_h12+n2*n2*n1*n1*C_h22-2*n1*n1*n1*n2*C_h23+2*n1*n1*n1*n2*C_h13+2*n2*n2*n2*n1*C_h23-4*n1*n1*n2*n2*C_h33;
	      C[0][2] = C[2][0] = n1*n1*n1*n2*C_h11+n1*n2*n2*n2*C_h12-2*n1*n1*n2*n2*C_h12-n1*n1*n1*n2*C_h12-n1*n2*n2*n2*C_h22+3*n1*n1*n2*n2*C_h23+n1*n1*n1*n1*C_h13-n1*n1*C_h13*n2*n2-n2*n2*n2*n2*C_h23-2*n1*n1*n1*n2*C_h33+2*n1*n2*n2*n2*C_h33;
	      C[1][1] = n2*n2*n2*n2*C_h11+2*n1*n1*n2*n2*C_h12+2*n1*n2*n2*n2*C_h12+n1*n1*n1*n1*C_h22+4*n1*n1*n1*n2*C_h23+2*n1*n2*n2*n2*C_h13+4*n1*n1*n2*n2*C_h33;
	      C[1][2] = C[2][1] = n1*n2*n2*n2*C_h11+n1*n1*n1*n2*C_h12+2*n1*n1*n2*n2*C_h12-n1*n2*n2*n2*C_h12-n1*n1*n1*n2*C_h22-3*n1*n1*n2*n2*C_h23+n1*n1*C_h13*n2*n2-n2*n2*n2*n2*C_h13+n1*n1*n1*n1*C_h23+2*n1*n1*n1*n2*C_h33-2*n1*n2*n2*n2*C_h33;
	      C[2][2] = n1*n1*n2*n2*C_h11-2*n1*n1*n2*n2*C_h12+n1*n1*n1*n2*C_h12-n1*n2*n2*n2*C_h12+n2*n2*n1*n1*C_h22-2*n1*n1*n1*n2*C_h23+2*n2*n2*n2*n1*C_h23+n1*n1*n1*n2*C_h13-n1*n2*n2*n2*C_h13+C_h33*n1*n1*n1*n1-2*n1*n1*n2*n2*C_h33+C_h33*n2*n2*n2*n2;
	}


	if(3==dim){
	      m1 = self->m[0];
	      m2 = self->m[1];
	      m3 = self->m[2];
	      n3 = self->n[2];
	      q1 = self->q[0];
	      q2 = self->q[1];
	      q3 = self->q[2];
	      C_h44 = self->C44;
	      C_h55 = self->C55;
	      C_h66 = self->C66;

	      C[0][0] = C_h11*n1*n1*n1*n1+C_h22*m1*m1*m1*m1+C_h33*q1*q1*q1*q1+4*C_h44*n1*n1*m1*m1+4*C_h55*n1*n1*q1*q1+4*C_h66*q1*q1*m1*m1+2*C_h12*n1*n1*m1*m1+2*C_h12*n1*n1*q1*q1+2*C_h23*q1*q1*m1*m1;
	
	      C[1][1] = C_h11*n2*n2*n2*n2+C_h22*m2*m2*m2*m2+C_h33*q2*q2*q2*q2+4*C_h44*n2*n2*m2*m2+4*C_h55*n2*n2*q2*q2+4*C_h66*q2*q2*m2*m2+2*C_h12*n2*n2*m2*m2+2*C_h12*n2*n2*q2*q2+2*C_h23*q2*q2*m2*m2;

	      C[1][0]=C[0][1] = C_h11*n1*n1*n2*n2+C_h22*m1*m1*m2*m2+C_h33*q1*q1*q2*q2+4*C_h44*m1*n2*n1*m2+4*C_h55*q1*n2*n1*q2+4*C_h66*m1*q2*q1*m2+C_h12*(n1*n1*m2*m2+m1*m1*n2*n2)+C_h12*(n1*n1*q2*q2+q1*q1*n2*n2)+C_h23*(q1*q1*m2*m2+m1*m1*q2*q2);

	      C[2][2] = C_h11*n3*n3*n3*n3+C_h22*m3*m3*m3*m3+C_h33*q3*q3*q3*q3+4*C_h44*n3*n3*m3*m3+4*C_h55*n3*n3*q3*q3+4*C_h66*q3*q3*m3*m3+2*C_h12*n3*n3*m3*m3+2*C_h12*n3*n3*q3*q3+2*C_h23*q3*q3*m3*m3;

	      C[2][0]=C[0][2] = C_h11*n1*n1*n3*n3+C_h22*m1*m1*m3*m3+C_h33*q1*q1*q3*q3+4*C_h44*n1*m1*n3*m3+4*C_h55*n1*q1*n3*q3+4*C_h66*q1*m1*q3*m3+C_h12*(n1*n1*m3*m3+m1*m1*n3*n3)+C_h12*(n1*n1*q3*q3+q1*q1*n3*n3)+C_h23*(q1*q1*m3*m3+m1*m1*q3*q3);

	      C[2][1]=C[1][2] =  C_h11*n2*n2*n3*n3+C_h22*m2*m2*m3*m3+C_h33*q2*q2*q3*q3+4*C_h44*n2*m2*n3*m3+4*C_h55*n2*q2*n3*q3+4*C_h66*q2*m2*q3*m3+C_h12*(n2*n2*m3*m3+m2*m2*n3*n3)+C_h12*(n2*n2*q3*q3+q2*q2*n3*n3)+C_h23*(q2*q2*m3*m3+m2*m2*q3*q3);

	      C[5][0]=C[0][5] = C_h11*n1*n1*n3*n2+C_h22*m1*m1*m3*m2+C_h33*q1*q1*q3*q2+2*C_h44*m1*n1*(n2*m3+m2*n3)+2*C_h55*q1*n1*(n2*q3+q2*n3)+2*C_h66*m1*q1*(m2*q3+q2*m3)+C_h12*(n1*n1*m3*m2+m1*m1*n3*n2)+C_h12*(n1*n1*q3*q2+q1*q1*n3*n2)+C_h23*(q1*q1*m3*m2+m1*m1*q3*q2);

	      C[4][0]=C[0][4] = C_h11*n1*n1*n1*n3+C_h22*m1*m1*m1*m3+C_h33*q1*q1*q1*q3+2*C_h44*n1*m1*(n3*m1+n1*m3)+2*C_h55*n1*q1*(n3*q1+n1*q3)+2*C_h66*q1*m1*(m1*q3+q1*m3)+C_h12*n1*m1*(n3*m1+n1*m3)+C_h12*n1*q1*(n3*q1+n1*q3)+C_h23*q1*m1*(m1*q3+q1*m3);

	      C[3][0]=C[0][3] =  C_h11*n1*n1*n1*n2+C_h22*m1*m1*m1*m2+C_h33*q1*q1*q1*q2+2*C_h44*m1*n1*(n2*m1+m2*n1)+2*C_h55*q1*n1*(n2*q1+q2*n1)+2*C_h66*m1*q1*(q1*m2+q2*m1)+C_h12*m1*n1*(n2*m1+m2*n1)+C_h12*q1*n1*(n2*q1+q2*n1)+C_h23*m1*q1*(q1*m2+q2*m1);

	      C[5][1]=C[1][5] =  C_h11*n2*n2*n2*n3+C_h22*m2*m2*m2*m3+C_h33*q2*q2*q2*q3+2*C_h44*m2*n2*(n2*m3+m2*n3)+2*C_h55*q2*n2*(n2*q3+q2*n3)+2*C_h66*q2*m2*(m2*q3+q2*m3)+C_h12*m2*n2*(n2*m3+m2*n3)+C_h12*q2*n2*(n2*q3+q2*n3)+C_h23*q2*m2*(m2*q3+q2*m3);

	      C[4][1]=C[1][4] =  C_h11*n2*n2*n1*n3+C_h22*m2*m2*m1*m3+C_h33*q2*q2*q1*q3+2*C_h44*n2*m2*(n3*m1+n1*m3)+2*C_h55*n2*q2*(n3*q1+n1*q3)+2*C_h66*m2*q2*(m1*q3+q1*m3)+C_h12*(n2*n2*m1*m3+m2*m2*n1*n3)+C_h12*(n2*n2*q1*q3+q2*q2*n1*n3)+C_h23*(q2*q2*m1*m3+m2*m2*q1*q3);

	      C[3][1]=C[1][3] =  C_h11*n2*n2*n2*n1+C_h22*m2*m2*m2*m1+C_h33*q2*q2*q2*q1+2*C_h44*m2*n2*(n2*m1+m2*n1)+2*C_h55*q2*n2*(n2*q1+q2*n1)+2*C_h66*q2*m2*(q1*m2+q2*m1)+C_h12*m2*n2*(n2*m1+m2*n1)+C_h12*q2*n2*(n2*q1+q2*n1)+C_h23*q2*m2*(q1*m2+q2*m1);

	      C[5][2]=C[2][5] =  C_h11*n3*n3*n3*n2+C_h22*m3*m3*m3*m2+C_h33*q3*q3*q3*q2+2*C_h44*m3*n3*(n2*m3+m2*n3)+2*C_h55*q3*n3*(n2*q3+q2*n3)+2*C_h66*q3*m3*(m2*q3+q2*m3)+C_h12*m3*n3*(n2*m3+m2*n3)+C_h12*q3*n3*(n2*q3+q2*n3)+C_h23*q3*m3*(m2*q3+q2*m3);

	      C[4][2]=C[2][4] =  C_h11*n3*n3*n3*n1+C_h22*m3*m3*m3*m1+C_h33*q3*q3*q3*q1+2*C_h44*m3*n3*(n3*m1+n1*m3)+2*C_h55*q3*n3*(n3*q1+n1*q3)+2*C_h66*q3*m3*(m1*q3+q1*m3)+C_h12*m3*n3*(n3*m1+n1*m3)+C_h12*q3*n3*(n3*q1+n1*q3)+C_h23*q3*m3*(m1*q3+q1*m3);

	      C[3][2]=C[2][3] =  C_h11*n3*n3*n1*n2+C_h22*m3*m3*m1*m2+C_h33*q3*q3*q1*q2+2*C_h44*n3*m3*(n2*m1+m2*n1)+2*C_h55*n3*q3*(n2*q1+q2*n1)+2*C_h66*m3*q3*(q1*m2+q2*m1)+C_h12*(n3*n3*m1*m2+m3*m3*n1*n2)+C_h12*(n3*n3*q1*q2+q3*q3*n1*n2)+C_h23*(q3*q3*m1*m2+m3*m3*q1*q2);

	      C[5][5] =  C_h11*n2*n2*n3*n3+C_h22*m2*m2*m3*m3+C_h33*q2*q2*q3*q3+C_h44*(n2*m3+m2*n3)*(n2*m3+m2*n3)+C_h55*(n2*q3+q2*n3)*(n2*q3+q2*n3)+C_h66*(m2*q3+q2*m3)*(m2*q3+q2*m3)+2*C_h12*n2*m2*n3*m3+2*C_h12*n2*q2*n3*q3+2*C_h23*q2*m2*q3*m3;

	      C[4][5]=C[5][4] = C_h11*n3*n3*n1*n2+C_h22*m3*m3*m1*m2+C_h33*q3*q3*q1*q2+C_h44*(n3*m1+n1*m3)*(n2*m3+m2*n3)+C_h55*(n3*q1+n1*q3)*(n2*q3+q2*n3)+C_h66*(m1*q3+q1*m3)*(m2*q3+q2*m3)+C_h12*n3*m3*(n2*m1+m2*n1)+C_h12*n3*q3*(n2*q1+q2*n1)+C_h23*m3*q3*(q1*m2+q2*m1);

	      C[3][5]=C[5][3] =  C_h11*n2*n2*n1*n3+C_h22*m2*m2*m1*m3+C_h33*q2*q2*q1*q3+C_h44*(n2*m1+m2*n1)*(n2*m3+m2*n3)+C_h55*(n2*q1+q2*n1)*(n2*q3+q2*n3)+C_h66*(q1*m2+q2*m1)*(m2*q3+q2*m3)+C_h12*n2*m2*(n3*m1+n1*m3)+C_h12*n2*q2*(n3*q1+n1*q3)+C_h23*m2*q2*(m1*q3+q1*m3);

	      C[4][4] = C_h11*n1*n1*n3*n3+C_h22*m1*m1*m3*m3+C_h33*q1*q1*q3*q3+C_h44*(n3*m1+n1*m3)*(n3*m1+n1*m3)+C_h55*(n3*q1+n1*q3)*(n3*q1+n1*q3)+C_h66*(m1*q3+q1*m3)*(m1*q3+q1*m3)+2*C_h12*n1*m1*n3*m3+2*C_h12*n1*q1*n3*q3+2*C_h23*q1*m1*q3*m3;

	      C[3][4]=C[4][3] = C_h11*n1*n1*n3*n2+C_h22*m1*m1*m3*m2+C_h33*q1*q1*q3*q2+C_h44*(n2*m1+m2*n1)*(n3*m1+n1*m3)+C_h55*(n2*q1+q2*n1)*(n3*q1+n1*q3)+C_h66*(q1*m2+q2*m1)*(m1*q3+q1*m3)+C_h12*m1*n1*(n2*m3+m2*n3)+C_h12*q1*n1*(n2*q3+q2*n3)+C_h23*m1*q1*(m2*q3+q2*m3);

	      C[3][3] = C_h11*n1*n1*n2*n2+C_h22*m1*m1*m2*m2+C_h33*q1*q1*q2*q2+C_h44*(n2*m1+m2*n1)*(n2*m1+m2*n1)+C_h55*(n2*q1+q2*n1)*(n2*q1+q2*n1)+C_h66*(q1*m2+q2*m1)*(q1*m2+q2*m1)+2*C_h12*m1*n2*n1*m2+2*C_h12*q1*n2*n1*q2+2*C_h23*m1*q2*q1*m2;
	}

 
	constitutiveMatrix->isDiagonal = False;	   
           /*	   printf("In %s OK\n\n",__func__); */
           /*	   flag = 1; */
           /*	} */
	/*for(i=0;i<dim*(dim+1)/2;i++){ 
	 	   for(j=0;j<dim*(dim+1)/2;j++){
			printf("Matrix Data = %g [%d %d]\n",constitutiveMatrix->matrixData[i][j],i,j); 
 	   } 
 	}
	exit(0);*/ 
}

#if 0
void _Orthotropic_UpdateDrawParameters( void* rheology ) 
{
	Orthotropic*                   self               = (Orthotropic*) rheology;
	Particle_Index                   lParticle_I;
	Particle_Index                   particleLocalCount;
	StandardParticle*                	materialPoint;
	Dimension_Index                 dim                   = self->materialPointsSwarm->dim;
	Orthotropic_ParticleExt*            particleExt;   /* this new type is defined in Orthotropic.h */

	/* do stuff */
}
#endif


