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
                                     ForceVector* forceVector,
                                     Swarm* integrationSwarm,
                                     Stg_Shape* domainShape,
                                     FeMesh* geometryMesh,
                                     StressBC_Entry force,
                                     FiniteElementContext* context)
{
	DivergenceForce* self = (DivergenceForce*) _DivergenceForce_DefaultNew( name );

	DivergenceForce_InitAll( 
			self,
			forceVector,
			integrationSwarm,
                        domainShape,
                        geometryMesh,
                        force,
                        context);

	return self;
}

/* Creation implementation / Virtual constructor */
DivergenceForce* _DivergenceForce_New( 
		SizeT                                               sizeOfSelf,  
		Type                                                type,
		Stg_Class_DeleteFunction*                           _delete,
		Stg_Class_PrintFunction*                            _print,
		Stg_Class_CopyFunction*                             _copy, 
		Stg_Component_DefaultConstructorFunction*           _defaultConstructor,
		Stg_Component_ConstructFunction*                    _construct,
		Stg_Component_BuildFunction*                        _build,
		Stg_Component_InitialiseFunction*                   _initialise,
		Stg_Component_ExecuteFunction*                      _execute,
		Stg_Component_DestroyFunction*                      _destroy,
		ForceTerm_AssembleElementFunction*                  _assembleElement,		
		Name                                                name,
                FiniteElementContext*                               context)
{
	DivergenceForce* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(DivergenceForce) );
	self = (DivergenceForce*) _ForceTerm_New( 
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
		_assembleElement,
		name );
        self->context=context;
	
	return self;
}

void _DivergenceForce_Init(DivergenceForce* self,
                           Stg_Shape* domainShape, FeMesh *geometryMesh,
                           StressBC_Entry force,
                           FiniteElementContext* context)
{
	self->domainShape=domainShape;
        self->geometryMesh=geometryMesh;
        self->force=force;
        self->context=context;
}

void DivergenceForce_InitAll( 
		void*                                               forceTerm,
		ForceVector*                                        forceVector,
		Swarm*                                              integrationSwarm,
                Stg_Shape* domainShape,
                FeMesh* geometryMesh,
                StressBC_Entry force,
                FiniteElementContext* context)
{
	DivergenceForce* self = (DivergenceForce*) forceTerm;

	ForceTerm_InitAll( self, forceVector, integrationSwarm, NULL );
	_DivergenceForce_Init( self, domainShape, geometryMesh, force, context);
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
	return (void*)_DivergenceForce_New( 
		sizeof(DivergenceForce), 
		DivergenceForce_Type,
		_DivergenceForce_Delete,
		_DivergenceForce_Print,
		NULL,
		_DivergenceForce_DefaultNew,
		_DivergenceForce_Construct,
		_DivergenceForce_Build,
		_DivergenceForce_Initialise,
		_DivergenceForce_Execute,
		_DivergenceForce_Destroy,
		_DivergenceForce_AssembleElement,
		name,
                NULL);
}

void _DivergenceForce_Construct( void* forceTerm, Stg_ComponentFactory* cf, void* data ) {
	DivergenceForce*          self             = (DivergenceForce*)forceTerm;
	Dictionary*		dict;
        Stg_Shape* domainShape=NULL;
        FeMesh* geometryMesh=NULL;
        StressBC_Entry force;
        char *type;
        FiniteElementContext* context;

	/* Construct Parent */
	_ForceTerm_Construct( self, cf, data );

        context = (FiniteElementContext*)Stg_ComponentFactory_ConstructByName
          ( cf, "context", FiniteElementContext, True, data ) ;

	dict = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( cf->componentDict, self->name ) );
	domainShape =  Stg_ComponentFactory_ConstructByKey(  cf,  self->name,  "DomainShape", Stg_Shape,  True, data  ) ;
        type = Stg_ComponentFactory_GetString( cf, self->name, "force_type", "");

        if(!strcasecmp(type,"double") || !strcasecmp(type,"float"))
          {
            force.type = StressBC_Double;
            force.DoubleValue =
              Stg_ComponentFactory_GetDouble( cf, self->name, "force_value", 0.0);
          }
        else if(!strcasecmp(type,"func"))
          {
            char *funcName = Stg_ComponentFactory_GetString
              ( cf, self->name, "force_value", "");
            
            Index cfIndex;
            cfIndex = ConditionFunction_Register_GetIndex
              ( context->condFunc_Register, funcName);
            force.type = StressBC_ConditionFunction;
            if ( cfIndex == (unsigned)-1 ) {	
              Stream*	errorStr = Journal_Register( Error_Type, self->type );
              
              Journal_Printf( errorStr, "Error- in %s: While parsing "
                              "definition of DivergenceForce, the cond. func. "
                              " \"%s\" - wasn't found in the c.f. register.\n",
                              __func__, funcName );
              Journal_Printf( errorStr, "(Available functions in the C.F. register are: ");	
              ConditionFunction_Register_PrintNameOfEachFunc
                ( context->condFunc_Register, errorStr );
              Journal_Printf( errorStr, ")\n");	
              assert(0);
            }
            force.CFIndex = cfIndex;
          }
        else if(strlen(type)==0)
          {
            Stream* errorStr = Journal_Register( Error_Type, self->type );
            Journal_Printf( errorStr, "Error- in %s: While parsing "
                            "definition of DivergenceForce, force_type is not specified.\nSupported types are \"double\" and \"function\".\n",
                            __func__);
            assert(0);
          }
        else
          {
            Stream* errorStr = Journal_Register( Error_Type, self->type );
            Journal_Printf( errorStr, "Error- in %s: While parsing "
                            "definition of DivergenceForce, the type of condition \"%s\"\nis not supported.  Supported types are \"double\" and \"function\".\n",
                            __func__, type );
            assert(0);
          }
        
        geometryMesh=Stg_ComponentFactory_ConstructByKey(  cf,  self->name,  "GeometryMesh", FeMesh,  True, data  ) ;
        
	_DivergenceForce_Init( self, domainShape, geometryMesh, force, context);
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
  Node_ElementLocalIndex           eNode_I;
  Element_NodeIndex                elementNodeCount;
  Node_DomainIndex *elementNodes=NULL;

  double xi[3], force, factor;
  ElementType* geometryElementType;
  
  IArray *incidence;

  xi[0]=0;
  xi[1]=0;
  xi[2]=0;
  geometryElementType=FeMesh_GetElementType(self->geometryMesh,lElement_I);
  factor=ElementType_JacobianDeterminant( geometryElementType,
                                          self->geometryMesh,
                                          lElement_I,
                                          xi, forceVector->dim );
  incidence=IArray_New();
  Mesh_GetIncidence(mesh, Mesh_GetDimSize(mesh), lElement_I,
                    MT_VERTEX,incidence);
  elementNodeCount=IArray_GetSize(incidence);
  elementNodes=IArray_GetPtr(incidence);
  
  for( eNode_I = 0 ; eNode_I < elementNodeCount; eNode_I++ ) {
    if(Stg_Shape_IsCoordInside(self->domainShape,
                               Mesh_GetVertex(mesh,elementNodes[eNode_I])))
      {
        switch(self->force.type)
          {
          case StressBC_Double:
            force=self->force.DoubleValue;
            break;
          case StressBC_ConditionFunction:
            
            /* We use a variable number of zero "0", because
               we don't use the variable number and that one
               is always going to exist. */
            ConditionFunction_Apply
              (self->context->condFunc_Register->_cf[self->force.CFIndex],
               elementNodes[eNode_I],0,self->context,&force);
            break;
          }
        elForceVec[ eNode_I] += force*factor;
      }
  }
  NewClass_Delete(incidence);
}
