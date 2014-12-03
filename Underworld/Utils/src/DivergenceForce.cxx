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
** Copyright (C) 2008, California Institute of Technology
** Modified for DivergenceForce by Walter Landry
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
** $Id: /cig/src/Gale/Utils/src/DivergenceForce.c 1691 2007-03-13T18:13:42.248551Z boo  $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "DivergenceForce.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class */
const Type DivergenceForce_Type = "DivergenceForce";

DivergenceForce* DivergenceForce_New(Name name,
                                     FiniteElementContext*	context,
                                     ForceVector* forceVector,
                                     Swarm* integrationSwarm,
                                     Stg_Shape* domainShape,
                                     StressBC_Entry force)
{
	DivergenceForce* self = (DivergenceForce*) _DivergenceForce_DefaultNew( name );

	DivergenceForce_InitAll( 
			self,
                        context,
			forceVector,
			integrationSwarm,
                        domainShape,
                        force);

	return self;
}

/* Creation implementation / Virtual constructor */
DivergenceForce* _DivergenceForce_New( DIVERGENCEFORCE_DEFARGS)
{
	DivergenceForce* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(DivergenceForce) );
	self = (DivergenceForce*) _ForceTerm_New( FORCETERM_PASSARGS);
	
	return self;
}

void _DivergenceForce_Init(DivergenceForce* self,
                           Stg_Shape* domainShape,
                           StressBC_Entry force)
{
  self->isConstructed    = True;

  self->domainShape=domainShape;
  self->force=force;
}

void DivergenceForce_InitAll( 
		void*                                               forceTerm,
                FiniteElementContext*	context,
		ForceVector*                                        forceVector,
		Swarm*                                              integrationSwarm,
                Stg_Shape* domainShape,
                StressBC_Entry force)
{
	DivergenceForce* self = (DivergenceForce*) forceTerm;

	_ForceTerm_Init( self, context, forceVector, integrationSwarm, NULL );
	_DivergenceForce_Init( self, domainShape, force);
}

void _DivergenceForce_Delete( void* forceTerm ) {
	DivergenceForce* self = (DivergenceForce*)forceTerm;

	_ForceTerm_Delete( self );
}

void _DivergenceForce_Print( void* forceTerm, Stream* stream ) {
	DivergenceForce* self = (DivergenceForce*)forceTerm;
	
	_ForceTerm_Print( self, stream );

	/* General info */
}

void* _DivergenceForce_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(DivergenceForce);
	Type                                                      type = DivergenceForce_Type;
	Stg_Class_DeleteFunction*                              _delete = _DivergenceForce_Delete;
	Stg_Class_PrintFunction*                                _print = _DivergenceForce_Print;
	Stg_Class_CopyFunction*                                  _copy = NULL;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _DivergenceForce_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _DivergenceForce_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _DivergenceForce_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _DivergenceForce_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _DivergenceForce_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _DivergenceForce_Destroy;
	ForceTerm_AssembleElementFunction*            _assembleElement = _DivergenceForce_AssembleElement;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*)_DivergenceForce_New( DIVERGENCEFORCE_PASSARGS);
}

void _DivergenceForce_AssignFromXML( void* forceTerm, Stg_ComponentFactory* cf, void* data ) {
	DivergenceForce*          self             = (DivergenceForce*)forceTerm;
        Stg_Shape* domainShape=NULL;
        StressBC_Entry force;
        char *type;

	/* Construct Parent */
	_ForceTerm_AssignFromXML( self, cf, data );

	domainShape =  Stg_ComponentFactory_ConstructByKey(  cf,  self->name,  "DomainShape", Stg_Shape,  True, data  ) ;
        type = Stg_ComponentFactory_GetString(cf, self->name, "force_type",
                                              "equation");

        if(!strcasecmp(type,"double") || !strcasecmp(type,"float"))
          {
            force.type = StressBC_Double;
            force.DoubleValue =
              Stg_ComponentFactory_GetDouble( cf, self->name, "force_value", 0.0);
          }
        else if (strlen(type)==0 || 0==strcasecmp(type, "equation"))
          {
            force.type = StressBC_Equation;
            /* This leaks memory */
            force.equation =
              StG_Strdup(Stg_ComponentFactory_GetString(cf,self->name,
                                                        "force_value",""));
          }
        else if(!strcasecmp(type,"func"))
          {
            char *funcName = Stg_ComponentFactory_GetString
              ( cf, self->name, "force_value", "");
            
            Index cfIndex;
            cfIndex = ConditionFunction_Register_GetIndex
              ( condFunc_Register, funcName);
            force.type = StressBC_ConditionFunction;
            if ( cfIndex == (unsigned)-1 ) {	
              Stream*	errorStr = Journal_Register( Error_Type, self->type );
              
              Journal_Printf( errorStr, "Error- in %s: While parsing "
                              "definition of DivergenceForce, the cond. func. "
                              " \"%s\" - wasn't found in the c.f. register.\n",
                              __func__, funcName );
              Journal_Printf( errorStr, "(Available functions in the C.F. register are: ");	
              ConditionFunction_Register_PrintNameOfEachFunc
                ( condFunc_Register, errorStr );
              Journal_Printf( errorStr, ")\n");	
              assert(0);
            }
            force.CFIndex = cfIndex;
          }
        else
          {
            Stream* errorStr = Journal_Register( Error_Type, self->type );
            Journal_Printf( errorStr, "Error- in %s: While parsing "
                            "definition of DivergenceForce, the type of condition \"%s\"\nis not supported.  Supported types are \"equation\", \"double\" and \"func\".\n",
                            __func__, type );
            assert(0);
          }
        
	_DivergenceForce_Init( self, domainShape, force);
}

void _DivergenceForce_Build( void* forceTerm, void* data ) {
	DivergenceForce*               self               = (DivergenceForce*)forceTerm;
	_ForceTerm_Build( self, data );
}

void _DivergenceForce_Initialise( void* forceTerm, void* data ) {
	DivergenceForce*             self             = (DivergenceForce*)forceTerm;
	_ForceTerm_Initialise( self, data );
}

void _DivergenceForce_Execute( void* forceTerm, void* data ) {
	_ForceTerm_Execute( forceTerm, data );
}

void _DivergenceForce_Destroy( void* forceTerm, void* data ) {
	_ForceTerm_Destroy( forceTerm, data );
}


void _DivergenceForce_AssembleElement( void* forceTerm,
                                       ForceVector* forceVector, 
                                       Element_LocalIndex lElement_I, 
                                       double* elForceVec ) {
  DivergenceForce* self=(DivergenceForce*) forceTerm;
  FeMesh* mesh=forceVector->feVariable->feMesh;

  IArray *incidence=IArray_New();
  const MeshTopology_Dim dim(Mesh_GetDimSize(mesh));
  Mesh_GetIncidence(mesh,dim,lElement_I,MT_VERTEX,incidence);
  const int num_nodes=IArray_GetSize(incidence);
  Node_DomainIndex *nodes=(Node_DomainIndex*)IArray_GetPtr(incidence);
  
  /* Integrate the divergence force over the element using gaussian
     quadrature */
  double xi[]={-sqrt(3/5.0),0,sqrt(3/5.0)};
  double weights[]={5/9.0, 8/9.0, 5/9.0};

  for(int i=0;i<3;++i)
    {
      for(int j=0;j<3;++j)
        {
          for(int k=0;k<(dim==2 ? 1 : 3); ++k)
            {
              double local_coord[dim], global_coord[dim];
              local_coord[0]=xi[i];
              local_coord[1]=xi[j];
              if(dim==3)
                local_coord[2]=xi[k];
              double Ni[num_nodes];
              ElementType_EvaluateShapeFunctionsAt(mesh->feElType,local_coord,
                                                   Ni);
              double jac(ElementType_JacobianDeterminant(mesh->feElType,mesh,
                                                          lElement_I,
                                                          local_coord,dim));
              FeMesh_CoordLocalToGlobal(mesh,lElement_I,local_coord,
                                        global_coord);
              double force;
              switch(self->force.type)
                {
                case StressBC_Double:
                  force=self->force.DoubleValue;
                  break;
                case StressBC_Equation:
                  force=Equation_eval(global_coord,
                                      (DomainContext*)(self->context),
                                      self->force.equation);
                  break;
                case StressBC_ConditionFunction:
            
                  ConditionFunction_Apply
                    (condFunc_Register->_cf[self->force.CFIndex],
                     global_coord,self->context,&force);
                  break;
                default:
                  abort();
                  break;
                }
              for(int node=0;node<num_nodes;++node)
                {
                  if(Stg_Shape_IsCoordInside(self->domainShape,
                                             Mesh_GetVertex(mesh,nodes[node])))
                    {
                      double total_force=force*jac*weights[i]*weights[j]*Ni[node];
                      if(dim==3)
                        total_force*=weights[k];
                      elForceVec[node]+=total_force;
                    }
                }
            }
        }
    }
  NewClass_Delete(incidence);
}
