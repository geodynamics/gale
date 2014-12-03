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
** $Id: BuoyancyDampingTerm.c 964 2007-10-11 08:03:06Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <assert.h>
#include <string.h>
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PopulationControl/PopulationControl.h>
#include <PICellerator/Weights/Weights.h>
#include <PICellerator/MaterialPoints/MaterialPoints.h>
#include "types.h"
#include "BuoyancyForceTerm.h"
#include "BuoyancyDampingTerm.h"
#include "MaterialSwarmVariable.h"


/* Textual name of this class */
const Type BuoyancyDampingTerm_Type = "BuoyancyDampingTerm";


BuoyancyDampingTerm*
BuoyancyDampingTerm_New(Name name,
                        FiniteElementContext* context,
                        StiffnessMatrix* stiffnessMatrix,
                        Swarm* integrationSwarm,
                        BuoyancyForceTerm* buoyancyForceTerm)
{
  BuoyancyDampingTerm* self =
    (BuoyancyDampingTerm*) _BuoyancyDampingTerm_DefaultNew( name );

  self->isConstructed = True;
  _StiffnessMatrixTerm_Init(self,context,stiffnessMatrix,integrationSwarm,NULL);
                            
  _BuoyancyDampingTerm_Init(self, buoyancyForceTerm);

  return self;
}

/* Creation implementation / Virtual constructor */
BuoyancyDampingTerm* _BuoyancyDampingTerm_New(  BUOYANCYDAMPINGTERM_DEFARGS  )
{
  BuoyancyDampingTerm* self;
	
  /* Allocate memory */
  assert( _sizeOfSelf >= sizeof(BuoyancyDampingTerm) );
  self=(BuoyancyDampingTerm*)
    _StiffnessMatrixTerm_New(  STIFFNESSMATRIXTERM_PASSARGS  );

  /* Virtual info */
	
  return self;
}

void _BuoyancyDampingTerm_Init(void* matrixTerm,
                               BuoyancyForceTerm* buoyancyForceTerm)
{
  BuoyancyDampingTerm* self = (BuoyancyDampingTerm*)matrixTerm;
  self->buoyancyForceTerm=buoyancyForceTerm;
}

void _BuoyancyDampingTerm_Delete( void* matrixTerm ) {
  BuoyancyDampingTerm* self = (BuoyancyDampingTerm*)matrixTerm;
  _StiffnessMatrixTerm_Delete( self );
}

void _BuoyancyDampingTerm_Print( void* matrixTerm, Stream* stream ) {
  BuoyancyDampingTerm* self = (BuoyancyDampingTerm*)matrixTerm;
	
  _StiffnessMatrixTerm_Print( self, stream );

  /* General info */
}

void* _BuoyancyDampingTerm_DefaultNew(Name name)
{
  /* Variables set in this function */
  SizeT _sizeOfSelf = sizeof(BuoyancyDampingTerm);
  Type type = BuoyancyDampingTerm_Type;
  Stg_Class_DeleteFunction* _delete= _BuoyancyDampingTerm_Delete;
  Stg_Class_PrintFunction* _print= _BuoyancyDampingTerm_Print;
  Stg_Class_CopyFunction* _copy= NULL;
  Stg_Component_DefaultConstructorFunction* _defaultConstructor=
    _BuoyancyDampingTerm_DefaultNew;
  Stg_Component_ConstructFunction* _construct=_BuoyancyDampingTerm_AssignFromXML;
  Stg_Component_BuildFunction* _build=_BuoyancyDampingTerm_Build;
  Stg_Component_InitialiseFunction* _initialise=_BuoyancyDampingTerm_Initialise;
  Stg_Component_ExecuteFunction* _execute=_BuoyancyDampingTerm_Execute;
  Stg_Component_DestroyFunction* _destroy=_BuoyancyDampingTerm_Destroy;
  StiffnessMatrixTerm_AssembleElementFunction* _assembleElement=
    _BuoyancyDampingTerm_AssembleElement;

  /* Variables that are set to ZERO are variables that will be
     set either by the current _New function or another parent
     _New function further up the hierachy */
  AllocationType  nameAllocationType = NON_GLOBAL /* default
                                                     value
                                                     NON_GLOBAL */;

  return (void*)_BuoyancyDampingTerm_New(  BUOYANCYDAMPINGTERM_PASSARGS  );
}

void _BuoyancyDampingTerm_AssignFromXML(void* matrixTerm,
                                        Stg_ComponentFactory* cf, void* data)
{
  BuoyancyDampingTerm* self=(BuoyancyDampingTerm*)matrixTerm;
  /* Construct Parent */
  _StiffnessMatrixTerm_AssignFromXML( self, cf, data );

  BuoyancyForceTerm *force=(BuoyancyForceTerm*)
    Stg_ComponentFactory_ConstructByKey(cf,self->name,(Dictionary_Entry_Key)"BuoyancyForceTerm",
                                        BuoyancyForceTerm,False,data);
  _BuoyancyDampingTerm_Init(self, force);
}

void _BuoyancyDampingTerm_Build( void* matrixTerm, void* data )
{
  BuoyancyDampingTerm* self = (BuoyancyDampingTerm*)matrixTerm;
  _StiffnessMatrixTerm_Build( self, data );
  Stg_Component_Build( self->buoyancyForceTerm, data, False);
}

void _BuoyancyDampingTerm_Initialise( void* matrixTerm, void* data )
{
  BuoyancyDampingTerm* self = (BuoyancyDampingTerm*)matrixTerm;

  _StiffnessMatrixTerm_Initialise( self, data );
  Stg_Component_Initialise( self->buoyancyForceTerm, data, False );
}

void _BuoyancyDampingTerm_Execute( void* matrixTerm, void* data )
{
  _StiffnessMatrixTerm_Execute( matrixTerm, data );
}

void _BuoyancyDampingTerm_Destroy( void* matrixTerm, void* data )
{
  _StiffnessMatrixTerm_Destroy( matrixTerm, data );
}


void _BuoyancyDampingTerm_AssembleElement(void* matrixTerm,
                                          StiffnessMatrix* stiffnessMatrix, 
                                          Element_LocalIndex lElement_I, 
                                          SystemLinearEquations* sle,
                                          FiniteElementContext* context,
                                          double** elStiffMat)
{
  BuoyancyDampingTerm* self=Stg_CheckType(matrixTerm,BuoyancyDampingTerm);
  if(!self->buoyancyForceTerm || !self->buoyancyForceTerm->damping)
    return;

  FeVariable* variable1=stiffnessMatrix->rowVariable;
  Dimension_Index dim=stiffnessMatrix->dim;
  Dof_Index nodeDofCount=dim;

  /* Set the element type */
  ElementType* elementType=FeMesh_GetElementType(variable1->feMesh,lElement_I);
	
  double density, alpha;

  BuoyancyForceTerm_average_density_alpha
    (self->buoyancyForceTerm,
     (IntegrationPointsSwarm*)self->buoyancyForceTerm->integrationSwarm,
     variable1->feMesh,lElement_I,&density,&alpha);

  double **jac=Memory_Alloc_2DArray(double,dim,dim,(Name)"Temporary Jacobian");
                                    
  /* Gauss-Legendre integration */
  const double legendre_points[]={-sqrt(3/5.0),0,sqrt(3/5.0)};
  const double weights[]={5/9.0, 8/9.0, 5/9.0};

  for(unsigned int local_norm=0;local_norm<dim; ++local_norm)
    {
      const int tangent((local_norm+1)%dim), tangent2((tangent+1)%dim);

      for(int sgn=-1;sgn<=1;sgn+=2)
        {
          const unsigned int face((local_norm==0 ? 1 :
                                   (local_norm==1 ? 0 : 2))*2+(sgn+1)/2);
          IJK ijk;
          RegularMeshUtils_Node_1DTo3D
            (variable1->feMesh,
             Mesh_DomainToGlobal(variable1->feMesh,MT_VERTEX,
                                 elementType->faceNodes[face][0]),ijk);

          /* Do not apply a damping force to the bottom boundary.
             Assume that it is counteracted by whatever is below
             it. */
          if(local_norm==1 && sgn==-1 && ijk[1]==0)
            continue;
          
          for(int i=0; i<3; ++i)
            for(int j=0;j<(dim==2 ? 1 : 3); ++j)
              {
                double xi[dim];
                xi[local_norm]=sgn;
                xi[tangent]=legendre_points[i];
                if(dim==3)
                  {
                    xi[tangent2]=legendre_points[j];
                  }
                const int num_nodes(dim==2 ? 9 : 27);
                double Ni[num_nodes];
                ElementType_EvaluateShapeFunctionsAt(elementType,xi,Ni);
                ElementType_Jacobian(elementType,variable1->feMesh,
                                     lElement_I,xi,dim,jac,NULL);

                double gravity[dim];
                BuoyancyForceTerm_CalcGravity(self->buoyancyForceTerm,
                                              variable1->feMesh,lElement_I,
                                              xi,gravity);
                double temperature(0);
                if(alpha!=0 && self->buoyancyForceTerm->temperatureField)
                  FeVariable_InterpolateFromMeshLocalCoord
                    (self->buoyancyForceTerm->temperatureField,
                     variable1->feMesh, lElement_I, xi, &temperature);

                double rho(density*(1.0-alpha*temperature));
                for(unsigned int d=0;d<dim;++d)
                  {
                    double damping(-gravity[d]*rho*context->dt);
                    double geometric_factor;

                    /* For a top boundary, consider y=a*x.  So
                       dx/dxi>0, and dx/deta=1/a.  If a>0, Fy<0, so
                       need to change the sign with the d==local_norm
                       term. */

                    if(dim==2)
                      {
                        geometric_factor=jac[tangent][(d+1)%dim]
                          *sgn*(d==local_norm ? 1 : -1);
                      }
                    else
                      {
                        geometric_factor=
                          (jac[tangent][(d+1)%dim] * jac[tangent2][(d+2)%dim]
                           -jac[tangent][(d+2)%dim]*jac[tangent2][(d+1)%dim])
                          *sgn*(d==local_norm ? 1 : -1);
                      }
                    /* Now loop over the nodes of the element */
                    const int num_face_nodes(dim==2 ? (num_nodes==4 ? 2 : 3)
                                             : (num_nodes==8 ? 4 : 9));
                    for(int el_0=0;el_0<num_face_nodes;++el_0)
                      {
                        const int node_0(elementType->faceNodes[face][el_0]);
                        for(int el_1=0;el_1<num_face_nodes;++el_1)
                          {
                            const int node_1
                              (elementType->faceNodes[face][el_1]);

                            double total=damping*weights[i]*geometric_factor
                              *Ni[node_0]*Ni[node_1];
                            if(dim==3)
                              {
                                total*=weights[j];
                              }
                            elStiffMat[node_0*nodeDofCount+d]
                              [node_1*nodeDofCount+d] +=total;
                          }
                      }
                  }
              }
        }
    }
  Memory_Free( jac );
}
