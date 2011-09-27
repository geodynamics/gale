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
        _StressBC_GetValues(cf,self,"x",data);
        _StressBC_GetValues(cf,self,"y",data);
        _StressBC_GetValues(cf,self,"z",data);

        self->bottom_density=
          Stg_ComponentFactory_GetDouble(cf,self->name,"bottomDensity",0.0);
}

/* Gets the actual values used by the StressBC (e.g. a float or a function). */
void _StressBC_GetValues(Stg_ComponentFactory* cf, void *stressBC,
                         char *direction, void *data)
{
          StressBC*          self             = (StressBC*)stressBC;
          char temp_str[20];
          char *type;

          switch(*direction)
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

          strcat(strcpy(temp_str,direction),"_type");
          type=Stg_ComponentFactory_GetString( cf, self->name, temp_str, "");
          strcat(strcpy(temp_str,direction),"_value");

          if(!strcasecmp(type,"double") || !strcasecmp(type,"float"))
            {
              self->_entryTbl[self->numEntries].type = StressBC_Double;
              self->_entryTbl[self->numEntries].DoubleValue =
                Stg_ComponentFactory_GetDouble( cf, self->name, temp_str, 0.0);
              (self->numEntries)++;
            }
          else if(!strcasecmp(type,"func"))
            {
              char *funcName =
                Stg_ComponentFactory_GetString( cf, self->name, temp_str, "");
              Index	cfIndex;

              cfIndex = ConditionFunction_Register_GetIndex
                ( self->conFunc_Register, funcName);
              self->_entryTbl[self->numEntries].type =
                StressBC_ConditionFunction;
              if ( cfIndex == (unsigned)-1 ) {	
                Stream*	errorStr = Journal_Register( Error_Type, self->type );
                
                Journal_Printf( errorStr, "Error- in %s: While parsing "
                                "definition of StressBC (applies to wall \"%s\"), the cond. func. applied to "
                                "direction \"%s\" - \"%s\" - wasn't found in the c.f. register.\n",
                                __func__, WallEnumToStr[self->_wall],
                                direction, funcName );
                Journal_Printf( errorStr, "(Available functions in the C.F. register are: ");	
                ConditionFunction_Register_PrintNameOfEachFunc
                  ( self->conFunc_Register, errorStr );
                Journal_Printf( errorStr, ")\n");	
                assert(0);
              }	
              self->_entryTbl[self->numEntries].CFIndex = cfIndex;
              (self->numEntries)++;
            }
          else if(!strcasecmp(type,"HydrostaticTerm"))
            {
              self->_entryTbl[self->numEntries].type = StressBC_HydrostaticTerm;
              self->_entryTbl[self->numEntries].hydrostaticTerm =
                Stg_ComponentFactory_ConstructByKey( cf, self->name, temp_str,
                                                     HydrostaticTerm, True,
                                                     data ) ;
              (self->numEntries)++;
            }
          else if(strlen(type)!=0)
            {
              Stream* errorStr = Journal_Register( Error_Type, self->type );
              Journal_Printf( errorStr, "Error- in %s: While parsing "
                              "definition of StressBC (applies to wall \"%s\"), the type of condition \"%s\"\nis not supported.  Supported types are \"double\" and \"function\".",
                                __func__, WallEnumToStr[self->_wall],
                                type );
              assert(0);
            }
}

void _StressBC_Build( void* forceTerm, void* data ) {
	StressBC*               self               = (StressBC*)forceTerm;
	Name                             name;

	_ForceTerm_Build( self, data );
}

void _StressBC_Initialise( void* forceTerm, void* data ) {
	StressBC*             self             = (StressBC*)forceTerm;
	Index                          i;

	_ForceTerm_Initialise( self, data );

}

void _StressBC_Execute( void* forceTerm, void* data ) {
	_ForceTerm_Execute( forceTerm, data );
}

void _StressBC_Destroy( void* forceTerm, void* data ) {
	_ForceTerm_Destroy( forceTerm, data );
}

double StressBC_compute_node_area(Wall wall, FeMesh *mesh, Index lElement_I, 
                                  Dimension_Index dim,
                                  const int block_map[9][8],
                                  const int &block,
                                  int *elementNodes);

void _StressBC_AssembleElement( void* forceTerm, ForceVector* forceVector, Element_LocalIndex lElement_I, double* elForceVec ) {
  StressBC*               self               = (StressBC*) forceTerm;
  Dimension_Index                  dim                = forceVector->dim;
  FeMesh*              mesh               = forceVector->feVariable->feMesh;
  Grid*		grid;
  Node_ElementLocalIndex           eNode_I;
  ElementType*                     elementType;
  Dof_Index                        nodeDofCount;
  double                           stress, area;
  IJK			ijk;
  int overcount;
  IArray *incidence;
  int elementNodeCount;
  int *elementNodes;
  double *coord;

  FeVariable* velocityField = NULL;
  velocityField = (FeVariable*)FieldVariable_Register_GetByName
    ( self->context->fieldVariable_Register, "VelocityField" );

  elementType       = FeMesh_GetElementType( mesh, lElement_I );
  nodeDofCount      = dim;
  
  grid=*(Grid**)ExtensionManager_Get(mesh->info, mesh, 
                                     ExtensionManager_GetHandle(mesh->info,
                                                                "vertexGrid"));
  
  incidence=IArray_New();
  Mesh_GetIncidence(mesh, Mesh_GetDimSize(mesh), lElement_I,
                    MT_VERTEX,incidence);
  elementNodes=IArray_GetPtr(incidence);
  elementNodeCount=IArray_GetSize(incidence);
                              
  int block_map[9][8];

  block_map[0][0]=0;
  block_map[0][1]=1;
  block_map[0][2]=2;
  block_map[0][3]=3;
  block_map[0][4]=4;
  block_map[0][5]=5;
  block_map[0][6]=6;
  block_map[0][7]=7;

  block_map[1][0]=0;
  block_map[1][1]=1;
  block_map[1][2]=3;
  block_map[1][3]=4;
  block_map[1][4]=9;
  block_map[1][5]=10;
  block_map[1][6]=12;
  block_map[1][7]=13;

  for(int i=0;i<8;i++)
    {
      block_map[2][i]=block_map[1][i]+1;
      block_map[3][i]=block_map[1][i]+3;
      block_map[4][i]=block_map[1][i]+4;
      block_map[5][i]=block_map[1][i]+9;
      block_map[6][i]=block_map[1][i]+10;
      block_map[7][i]=block_map[1][i]+12;
      block_map[8][i]=block_map[1][i]+13;
    }

  int block_start(0), nblocks(1), nodes_per_block(1<<dim);
  if(dim==2 && elementNodeCount==9)
    {
      block_start=1;
      nblocks=4;
    }
  else if(dim==3 && elementNodeCount==27)
    {
      block_start=1;
      nblocks=8;
    }
  for(int block=block_start;block<block_start+nblocks;++block)
    {
      area=StressBC_compute_node_area(self->_wall,mesh,lElement_I,dim,
                                      block_map,block,elementNodes);

      /* Apply the stress */
      for( eNode_I = 0 ; eNode_I < nodes_per_block; eNode_I++ ) {
        /* Make sure that we are on the boundary */
        int condition, entry_I;
        ConditionFunction* cf;
        RegularMeshUtils_Node_1DTo3D
          (mesh,
           Mesh_DomainToGlobal(mesh,MT_VERTEX,
                               elementNodes[block_map[block][eNode_I]]),
           ijk);
    
        switch(self->_wall)
          {
          case Wall_Left:
            condition=(ijk[0] == 0);
            break;
          case Wall_Right:
            condition=(ijk[0] == ( grid->sizes[0] - 1 ));
            break;
          case Wall_Bottom:
            condition=(ijk[1] == 0);
            break;
          case Wall_Top:
            condition=(ijk[1] == ( grid->sizes[1] - 1 ));
            break;
          case Wall_Front:
            condition=(ijk[2] == 0);
            break;
          case Wall_Back:
            condition=(ijk[2] == ( grid->sizes[2] - 1 ));
            break;
          }
    
        if(condition)
          {

            for(entry_I=0; entry_I<self->numEntries; ++entry_I)
              {
                switch(self->_entryTbl[entry_I].type)
                  {
                  case StressBC_Double:
                    stress=self->_entryTbl[entry_I].DoubleValue;
                    break;
                  case StressBC_ConditionFunction:
                    cf=self->conFunc_Register->
                      _cf[self->_entryTbl[entry_I].CFIndex];
                
                    /* We use a variable number of zero "0", because
                       we don't use the variable number and that one
                       is always going to exist. */
                    ConditionFunction_Apply(cf,elementNodes[block_map[block][eNode_I]],
                                            0,self->context,&stress);
                    break;
                  case StressBC_HydrostaticTerm:
                    if(self->_wall!=Wall_Top && self->_wall!=Wall_Bottom)
                      {
                        Stream* errorStr=Journal_Register( Error_Type, self->type );
                        Journal_Firewall(0,errorStr,"You can only apply a HydrostaticTerm StressBC to the top or bottom wall.\nYou applied it to the %s wall",WallVC_WallEnumToStr[self->_wall]);
                      }

                    coord=Mesh_GetVertex(mesh,elementNodes[block_map[block][eNode_I]]);
                    stress=
                      -HydrostaticTerm_Pressure(self->_entryTbl[entry_I].hydrostaticTerm,
                                                coord);
                    /* For the bottom, we need to add in the effects of
                       material outside the mesh. */

                    if(self->_wall==Wall_Bottom)
                      {
                        Coord bottom;
                        double dy;
                        bottom[0]=coord[0];
                        bottom[1]=0;
                        bottom[2]=coord[2];
                        dy=bottom[1] - coord[1];
                        stress+=HydrostaticTerm_Pressure(self->_entryTbl[entry_I].hydrostaticTerm,bottom);
                        stress-=self->bottom_density
                          * self->_entryTbl[entry_I].hydrostaticTerm->gravity
                          * dy;
                      }
                    
                    if(dim==2)
                      {
                        double dx, dy, *coord0, *coord1;
                        /* StGermain uses the ordering

                           0:x,y
                           1:x+,y
                           2:x+,y+
                           3:x,y+

                           So the top two vertices are 2 and 3 */

                        coord0=Mesh_GetVertex(mesh,elementNodes[block_map[block][3]]);
                        coord1=Mesh_GetVertex(mesh,elementNodes[block_map[block][2]]);

                        dx=coord1[0]-coord0[0];
                        dy=-(coord1[1]-coord0[1]);

                        elForceVec[block_map[block][eNode_I]*nodeDofCount + I_AXIS]+=
                          stress*dy/overcount;
                        elForceVec[block_map[block][eNode_I]*nodeDofCount + J_AXIS]+=
                          stress*dx/overcount;
                      }
                    else
                      {
                        double *coord0,*coord1,*coord2,*coord3,
                          vector0[3],vector1[3],normal[3];
                    
                        /* Decompose the quadrilateral into two triangles.
                           For each triangle, take the cross product and
                           apply the force in that direction. */

                        coord0=Mesh_GetVertex(mesh,elementNodes[block_map[block][2]]);
                        coord1=Mesh_GetVertex(mesh,elementNodes[block_map[block][3]]);
                        coord2=Mesh_GetVertex(mesh,elementNodes[block_map[block][7]]);
                        coord3=Mesh_GetVertex(mesh,elementNodes[block_map[block][6]]);
                    
                        StGermain_VectorSubtraction( vector0, coord1, coord0, dim) ;
                        StGermain_VectorSubtraction( vector1, coord2, coord1, dim) ;
                        StGermain_VectorCrossProduct( normal, vector0, vector1 );

                        elForceVec[block_map[block][eNode_I]*nodeDofCount + I_AXIS]+=
                          stress*normal[0]/(2*overcount);
                        elForceVec[block_map[block][eNode_I]*nodeDofCount + J_AXIS]+=
                          stress*normal[1]/(2*overcount);
                        elForceVec[block_map[block][eNode_I]*nodeDofCount + K_AXIS]+=
                          stress*normal[2]/(2*overcount);

                        StGermain_VectorSubtraction( vector0, coord2, coord1, dim) ;
                        StGermain_VectorSubtraction( vector1, coord3, coord2, dim) ;
                        StGermain_VectorCrossProduct( normal, vector0, vector1 );

                        elForceVec[block_map[block][eNode_I]*nodeDofCount + I_AXIS]+=
                          stress*normal[0]/(2*overcount);;
                        elForceVec[block_map[block][eNode_I]*nodeDofCount + J_AXIS]+=
                          stress*normal[1]/(2*overcount);;
                        elForceVec[block_map[block][eNode_I]*nodeDofCount + K_AXIS]+=
                          stress*normal[2]/(2*overcount);;
                      }
                    break;
                  }
                /* Actually apply the stress */
                if(self->_entryTbl[entry_I].type!=StressBC_HydrostaticTerm)
                  elForceVec[block_map[block][eNode_I]*nodeDofCount
                             + self->_entryTbl[entry_I].axis]+=
                    stress*area;
              }
          }
      }
    }
  NewClass_Delete(incidence);
}

double StressBC_compute_node_area(Wall wall, FeMesh *mesh, Index lElement_I, 
                                  Dimension_Index dim,
                                  const int block_map[9][8],
                                  const int &block,
                                  int *elementNodes)
{
  /* Compute the area of the face. */
  if(dim==2)
    {
      double *coord1, *coord2;
      int lower,upper,direction;
      switch(wall)
        {
        case Wall_Left:
          lower=0;
          upper=3;
          direction=1;
          break;
        case Wall_Right:
          lower=1;
          upper=2;
          direction=1;
          break;
        case Wall_Bottom:
          lower=0;
          upper=1;
          direction=0;
          break;
        case Wall_Top:
          lower=3;
          upper=2;
          direction=0;
          break;
        }

      coord1=Mesh_GetVertex(mesh,elementNodes[block_map[block][lower]]);
      coord2=Mesh_GetVertex(mesh,elementNodes[block_map[block][upper]]);
      /* Scale by 2 since each node only claims 1/2 of the area of a
         face. */
      return (coord2[direction]-coord1[direction])/2;
    }
  else
    {
      double *coord1, *coord2, *coord3, *coord4;
      int first, second, third, fourth;
                              
      /* StGermain uses the ordering
         0: x,y,z
         1: x+,y,z
         2: x,y+,z
         3: x+,y+,z
         4: x,y,z+
         5: x+,y,z+
         6: x,y+,z+
         7: x+,y+,z+
      */

      /* Get the indices for which wall we want to get the area of. */
      switch(wall)
        {
        case Wall_Left:
          first=0;
          second=2;
          third=6;
          fourth=4;
          break;
        case Wall_Right:
          first=1;
          second=3;
          third=7;
          fourth=5;
          break;
        case Wall_Bottom:
          first=0;
          second=1;
          third=5;
          fourth=4;
          break;
        case Wall_Top:
          first=2;
          second=3;
          third=7;
          fourth=6;
          break;
        case Wall_Front:
          first=0;
          second=1;
          third=4;
          fourth=3;
          break;
        case Wall_Back:
          first=4;
          second=5;
          third=7;
          fourth=6;
          break;
        }

      coord1=Mesh_GetVertex(mesh,elementNodes[block_map[block][first]]);
      coord2=Mesh_GetVertex(mesh,elementNodes[block_map[block][second]]);
      coord3=Mesh_GetVertex(mesh,elementNodes[block_map[block][third]]);
      coord4=Mesh_GetVertex(mesh,elementNodes[block_map[block][fourth]]);

      /* Scale by 4 since each node only claims 1/4 of the area of a
         face. */
      return StGermain_ConvexQuadrilateralArea(coord1,coord2,coord3,coord4,
                                               dim)/4;
    }
}


