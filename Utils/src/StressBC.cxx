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
**  Modified 2006-2010 by Walter Landry (California Institute of Technology)
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
** $Id$
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <PICellerator/PopulationControl/PopulationControl.h>
#include <PICellerator/Weights/Weights.h>
#include <PICellerator/MaterialPoints/MaterialPoints.h>

#include "PICellerator/Utils/types.h"
#include "PICellerator/Utils/MaterialSwarmVariable.h"

#include "types.h"
#include "StressBC.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class */
const Type StressBC_Type = "StressBC";

extern Name WallEnumToStr[Wall_Size];

StressBC* StressBC_New( 
		Name                                      name,
                FiniteElementContext*	context,
		ForceVector*                              forceVector,
		Swarm*                                    integrationSwarm,
                ConditionFunction_Register*               conFunc_Register)
{
	StressBC* self = (StressBC*) _StressBC_DefaultNew( name );

	StressBC_InitAll( 
			self,
                        context,
			forceVector,
			integrationSwarm,
                        conFunc_Register);

	return self;
}

/* Creation implementation / Virtual constructor */
StressBC* _StressBC_New( STRESSBC_DEFARGS)
{
	StressBC* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(StressBC) );
	self = (StressBC*) _ForceTerm_New( FORCETERM_PASSARGS);
	
        self->conFunc_Register=condFunc_Register;
	
	return self;
}

void _StressBC_Init( 
		StressBC*                                  self, 
                ConditionFunction_Register*		   conFunc_Register)
{
  self->isConstructed    = True;
  self->numEntries        = 0;
  self->conFunc_Register=conFunc_Register;
}

void StressBC_InitAll( 
		void*                                               forceTerm,
                FiniteElementContext*	context,
		ForceVector*                                        forceVector,
		Swarm*                                              integrationSwarm,
                ConditionFunction_Register*			    conFunc_Register)
{
	StressBC* self = (StressBC*) forceTerm;

	_ForceTerm_Init( self, context, forceVector, integrationSwarm, NULL );
	_StressBC_Init( self, conFunc_Register );
}

void _StressBC_Delete( void* forceTerm ) {
	StressBC* self = (StressBC*)forceTerm;

	_ForceTerm_Delete( self );
}

void _StressBC_Print( void* forceTerm, Stream* stream ) {
	StressBC* self = (StressBC*)forceTerm;
	
	_ForceTerm_Print( self, stream );
}

void* _StressBC_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(StressBC);
	Type                                                      type = StressBC_Type;
	Stg_Class_DeleteFunction*                              _delete = _StressBC_Delete;
	Stg_Class_PrintFunction*                                _print = _StressBC_Print;
	Stg_Class_CopyFunction*                                  _copy = NULL;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _StressBC_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _StressBC_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _StressBC_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _StressBC_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _StressBC_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _StressBC_Destroy;
	ForceTerm_AssembleElementFunction*            _assembleElement = _StressBC_AssembleElement;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*)_StressBC_New( STRESSBC_PASSARGS);
}

void _StressBC_AssignFromXML( void* forceTerm, Stg_ComponentFactory* cf, void* data ) {
	StressBC*          self             = (StressBC*)forceTerm;

	/* Construct Parent */
	_ForceTerm_AssignFromXML( self, cf, data );

	_StressBC_Init( self, condFunc_Register );
        {
          char*	wallStr;
		
          /* Obtain which wall */
          wallStr = Stg_ComponentFactory_GetString( cf, self->name,
                                                    "wall", "");
          if (!strcasecmp(wallStr, "back"))
            self->_wall = Wall_Back;
          else if (!strcasecmp(wallStr, "left"))
            self->_wall = Wall_Left;
          else if (!strcasecmp(wallStr, "bottom"))
            self->_wall = Wall_Bottom;
          else if (!strcasecmp(wallStr, "right"))
            self->_wall = Wall_Right;
          else if (!strcasecmp(wallStr, "top"))
            self->_wall = Wall_Top;
          else if (!strcasecmp(wallStr, "front"))
            self->_wall = Wall_Front;
          else {
            assert( 0 );
            self->_wall = Wall_Size; /* invalid entry */
          }
        }
        _StressBC_GetValues(cf,self,'x',data);
        _StressBC_GetValues(cf,self,'y',data);
        _StressBC_GetValues(cf,self,'z',data);

        self->bottom_density=
          Stg_ComponentFactory_GetDouble(cf,self->name,"bottomDensity",0.0);
}

/* Gets the actual values used by the StressBC (e.g. a float or a function). */
void _StressBC_GetValues(Stg_ComponentFactory* cf, void *stressBC,
                         const char direction, void *data)
{
          StressBC*          self             = (StressBC*)stressBC;
          char *type;

          switch(direction)
            {
            case 'x':
              self->_entryTbl[self->numEntries].axis=I_AXIS;
              break;
            case 'y':
              self->_entryTbl[self->numEntries].axis=J_AXIS;
              break;
            case 'z':
              self->_entryTbl[self->numEntries].axis=K_AXIS;
              break;
            }

          std::string temp_str(direction + std::string("_type"));
          type=Stg_ComponentFactory_GetString( cf, self->name, temp_str.c_str(), "");
          temp_str=direction+std::string("_value");

          if(!strcasecmp(type,"double") || !strcasecmp(type,"float"))
            {
              self->_entryTbl[self->numEntries].type = StressBC_Double;
              self->_entryTbl[self->numEntries].DoubleValue =
                Stg_ComponentFactory_GetDouble( cf, self->name, temp_str.c_str(), 0.0);
              (self->numEntries)++;
            }
          else if(!strcasecmp(type,"func"))
            {
              char *funcName =
                Stg_ComponentFactory_GetString( cf, self->name, temp_str.c_str(), "");
              int cfIndex = ConditionFunction_Register_GetIndex
                ( self->conFunc_Register, funcName);
              self->_entryTbl[self->numEntries].type =
                StressBC_ConditionFunction;
              if ( cfIndex == -1 ) {	
                Stream*	errorStr = Journal_Register( Error_Type, self->type );
                
                Journal_Printf(errorStr, "Error- in %s: While parsing "
                               "definition of StressBC (applies to wall \"%s\"), the cond. func. applied to "
                               "direction \"%c\": \"%s\" - wasn't found in the c.f. register.\n",
                               __func__, WallEnumToStr[self->_wall],
                               direction, funcName);
                Journal_Printf(errorStr,
                               "(Available functions in the C.F. register are: ");	
                ConditionFunction_Register_PrintNameOfEachFunc
                  ( self->conFunc_Register, errorStr );
                Journal_Printf( errorStr, ")\n");	
                abort();
              }	
              self->_entryTbl[self->numEntries].CFIndex = cfIndex;
              (self->numEntries)++;
            }
          else if(!strcasecmp(type,"HydrostaticTerm"))
            {
              self->_entryTbl[self->numEntries].type = StressBC_HydrostaticTerm;
              self->_entryTbl[self->numEntries].hydrostaticTerm =
                Stg_ComponentFactory_ConstructByKey( cf, self->name, temp_str.c_str(),
                                                     HydrostaticTerm, True,
                                                     data ) ;
              (self->numEntries)++;
            }
          else if(strlen(type)!=0)
            {
              Stream* errorStr = Journal_Register( Error_Type, self->type );
              Journal_Firewall( 0, errorStr, "Error- in %s: While parsing "
                              "definition of StressBC (applies to wall \"%s\"), the type of condition \"%s\"\nis not supported.  Supported types are \"double\" and \"function\".",
                                __func__, WallEnumToStr[self->_wall],
                                type );
            }
}

void _StressBC_Build( void* forceTerm, void* data ) {
  StressBC*               self               = (StressBC*)forceTerm;
  _ForceTerm_Build( self, data );
}

void _StressBC_Initialise( void* forceTerm, void* data ) {
  StressBC*             self             = (StressBC*)forceTerm;
  _ForceTerm_Initialise( self, data );
}

void _StressBC_Execute( void* forceTerm, void* data ) {
  _ForceTerm_Execute( forceTerm, data );
}

void _StressBC_Destroy( void* forceTerm, void* data ) {
  _ForceTerm_Destroy( forceTerm, data );
}

namespace {
  bool on_wall(const Wall &wall, Grid *grid, IJK ijk);
  inline double square(const double &a)
  {
    return a*a;
  }
}

void _StressBC_AssembleElement( void* forceTerm, ForceVector* forceVector, Element_LocalIndex lElement_I, double* elForceVec ) {
  StressBC*               self               = (StressBC*) forceTerm;
  Dimension_Index                  dim                = forceVector->dim;
  FeMesh*              mesh               = forceVector->feVariable->feMesh;
  Grid*		grid;
  IArray *incidence;
  int *elementNodes;
  int nodeDofCount(dim);
  
  grid=*(Grid**)ExtensionManager_Get(mesh->info, mesh, 
                                     ExtensionManager_GetHandle(mesh->info,
                                                                "vertexGrid"));
  incidence=IArray_New();
  Mesh_GetIncidence(mesh, Mesh_GetDimSize(mesh), lElement_I,
                    MT_VERTEX,incidence);
  elementNodes=IArray_GetPtr(incidence);
  const int num_nodes=IArray_GetSize(incidence);

  /* Make sure that we are on the boundary */
  bool wall_element=false;
  for(int el=0;el<num_nodes;++el)
    {
      IJK ijk;
      RegularMeshUtils_Node_1DTo3D
        (mesh,Mesh_DomainToGlobal(mesh,MT_VERTEX,elementNodes[el]),ijk);
      if(on_wall(self->_wall,grid,ijk))
        {
          wall_element=true;
          break;
        }
    }
  if(!wall_element)
    return;

  /* Set up directions particular to each wall */
  int x,y,z;
  int local_bottom(-1);
  switch(self->_wall)
    {
      case Wall_Right:
        local_bottom=1;
      case Wall_Left:
        x=1;
        y=2;
        z=0;
        break;
      case Wall_Top:
        local_bottom=1;
      case Wall_Bottom:
        x=0;
        y=2;
        z=1;
        break;
      case Wall_Front:
        local_bottom=1;
      case Wall_Back:
        x=0;
        y=1;
        z=2;
        break;
      default:
        abort();
      }
  /* Fix up for 2D */
  if(dim==2)
    y=1-x;

  /* For each gaussian integration point, compute the value of the
     Jacobian */
  double **jacobian=Memory_Alloc_2DArray(double,dim,dim,
                                         (Name)"Temporary Jacobian");
  double xi[]={-sqrt(3/5.0),0,sqrt(3/5.0)};
  double weights[]={5/9.0, 8/9.0, 5/9.0};

  for(int i=0;i<3;++i)
    {
      for(int j=0;j<(dim==2 ? 1 : 3); ++j)
        {
          double local_coord[dim], global_coord[dim];
          local_coord[x]=xi[i];
          if(dim==2)
            {
              local_coord[y]=local_bottom;
            }
          else
            {
              local_coord[y]=xi[j];
              local_coord[z]=local_bottom;
            }
          double Ni[num_nodes];
          ElementType_EvaluateShapeFunctionsAt(mesh->feElType,local_coord,Ni);
          ElementType_Jacobian_AxisIndependent(mesh->feElType,mesh,lElement_I,
                                               local_coord,dim,jacobian,NULL,
                                               x,y,z);
          double geometric_factor;
          if(dim==2)
            {
              geometric_factor=sqrt(jacobian[0][x]*jacobian[0][x]
                                    + jacobian[1][x]*jacobian[1][x]);
            }
          else
            {
              geometric_factor=
                sqrt(square(jacobian[0][x]*jacobian[1][y]
                            - jacobian[1][x]*jacobian[0][y])
                     + square(jacobian[0][x]*jacobian[2][y]
                              - jacobian[2][x]*jacobian[0][y])
                     + square(jacobian[1][x]*jacobian[2][y]
                              - jacobian[2][x]*jacobian[1][y]));
            }
          

          FeMesh_CoordLocalToGlobal(mesh,lElement_I,local_coord,global_coord);

          /* Get the applied stress */
          ConditionFunction* cf;
          for(int entry_I=0; entry_I<self->numEntries; ++entry_I)
            {
              double stress;
              switch(self->_entryTbl[entry_I].type)
                {
                case StressBC_Double:
                  stress=self->_entryTbl[entry_I].DoubleValue;
                  break;
                case StressBC_ConditionFunction:
                  cf=self->conFunc_Register->
                    _cf[self->_entryTbl[entry_I].CFIndex];

                  /* We use a variable number of zero "0",
                     because we don't use the variable number
                     and that one is always going to exist. */
                  ConditionFunction_Apply(cf,global_coord,
                                          self->context,&stress);
                  break;
                default:
                  abort();
                  break;
                }

              /* Now loop over the nodes of the element */
              for(int el=0;el<num_nodes;++el)
                {
                  /* Make sure that node is on the boundary */
                  IJK ijk;
                  RegularMeshUtils_Node_1DTo3D
                    (mesh,Mesh_DomainToGlobal(mesh,MT_VERTEX,elementNodes[el]),
                     ijk);
                  if(on_wall(self->_wall,grid,ijk))
                    {
                      double total=stress*weights[i]*geometric_factor*Ni[el];
                      if(dim==3)
                        {
                          total*=weights[j];
                        }
                      elForceVec[el*nodeDofCount+self->_entryTbl[entry_I].axis]+=
                        total;
                    }
                }
            }
        }
    }
  Memory_Free( jacobian );
  NewClass_Delete(incidence);
}

namespace {
  bool on_wall(const Wall &wall, Grid *grid, IJK ijk)
  {
    bool result(false);
    switch(wall)
      {
      case Wall_Left:
        result=(ijk[0] == 0);
        break;
      case Wall_Right:
        result=(ijk[0] == ( grid->sizes[0] - 1 ));
        break;
      case Wall_Bottom:
        result=(ijk[1] == 0);
        break;
      case Wall_Top:
        result=(ijk[1] == ( grid->sizes[1] - 1 ));
        break;
      case Wall_Front:
        result=(ijk[2] == 0);
        break;
      case Wall_Back:
        result=(ijk[2] == ( grid->sizes[2] - 1 ));
        break;
      default:
        abort();
      }
    return result;
  }

}
