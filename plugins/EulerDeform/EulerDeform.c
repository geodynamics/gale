/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	David May, PhD Student Monash University, VPAC. (david.may@sci.maths.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
** $Id: LevelSetPlg.c 200 2005-07-08 08:24:41Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include "types.h"
#include "Context.h"
#include "EulerDeform.h"


Name		EULERDEFORM_PLUGIN_TAG = "EulerDeform";
const Type		Underworld_EulerDeform_Type = "EulerDeform";
ExtensionInfo_Index	EulerDeform_ContextHandle;


Index Underworld_EulerDeform_Register( PluginsManager* pluginsMgr ) {
	return PluginsManager_Submit( pluginsMgr, Underworld_EulerDeform_Type, (Name)"0", _Underworld_EulerDeform_DefaultNew );
}


void* _Underworld_EulerDeform_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(Codelet);
	Type                                                      type = Underworld_EulerDeform_Type;
	Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
	Stg_Class_PrintFunction*                                _print = _Codelet_Print;
	Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Underworld_EulerDeform_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _Underworld_EulerDeform_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _Underworld_EulerDeform_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _Codelet_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _Codelet_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _Underworld_EulerDeform_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return _Codelet_New(  CODELET_PASSARGS   );
}


void _Underworld_EulerDeform_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
	Codelet*					ed = (Codelet*)component;
	UnderworldContext*	uwCtx;
	EulerDeform_Context*	edCtx;

	assert( component );
	assert( cf );

	Journal_DPrintf( Underworld_Debug, "In: %s( void* )\n", __func__ );

	/* Retrieve context. */
	uwCtx = (UnderworldContext*)Stg_ComponentFactory_ConstructByName( cf, (Name)"context", UnderworldContext, True, data );
	ed->context = (AbstractContext* )uwCtx;

	/* Create new context. */
	EulerDeform_ContextHandle = ExtensionManager_Add( uwCtx->extensionMgr, (Name)Underworld_EulerDeform_Type, sizeof(EulerDeform_Context)  );
	edCtx = (EulerDeform_Context*)ExtensionManager_Get( uwCtx->extensionMgr, uwCtx, EulerDeform_ContextHandle );
	memset( edCtx, 0, sizeof(EulerDeform_Context) );
	edCtx->ctx = (AbstractContext*)uwCtx;

	/* Get the time integrator. */
	edCtx->timeIntegrator = Stg_ComponentFactory_ConstructByName( cf, (Name)"timeIntegrator", TimeIntegrator, True, data  );

	/* Grab the ArtDisplacementField from the dictionary */
	edCtx->artDField = Stg_ComponentFactory_ConstructByName( cf, (Name)"ArtDisplacementField", FeVariable, False, data  );
}


void _Underworld_EulerDeform_Build( void* component, void* data ) {
	Codelet*						ed	= (Codelet*)component;
	UnderworldContext*		uwCtx	= (UnderworldContext*)ed->context;
	EulerDeform_Context*		edCtx;
	Variable*					crdVar;
	TimeIntegrand*			crdAdvector;
	Stg_Component*				tiData[2];
	unsigned						sys_i;
	Dictionary_Entry_Value*	edDict;
	Dictionary_Entry_Value*	sysLst;

	assert( component );
	assert( uwCtx );

	edCtx = (EulerDeform_Context*)ExtensionManager_Get( uwCtx->extensionMgr, uwCtx, EulerDeform_ContextHandle );

	/* Get the dictionary. */
	edDict = Dictionary_Get( uwCtx->dictionary, (Dictionary_Entry_Key)"EulerDeform" );
	if( !edDict  ) {
		return;
	}

	/* Read system list. */
	sysLst = Dictionary_Entry_Value_GetMember( edDict, (Dictionary_Entry_Key)"systems" );
	if( sysLst ) {
		unsigned	sys_i;

		/* Allocate for systems. */
		edCtx->nSystems = Dictionary_Entry_Value_GetCount( sysLst  );
		edCtx->systems = Memory_Alloc_Array( EulerDeform_System, edCtx->nSystems, "EulerDeform->systems" );
		memset( edCtx->systems, 0, sizeof(EulerDeform_System) * edCtx->nSystems );

		for( sys_i = 0; sys_i < edCtx->nSystems; sys_i++ ) {
			EulerDeform_System*		sys = edCtx->systems + sys_i;
			Dictionary*					sysDict;
			Dictionary_Entry_Value*	varLst;
			char*							meshName;
			char*							remesherName;
			char*							velFieldName;
			char*							name;

			/* Get the dictionary for this system. */
			sysDict = Dictionary_Entry_Value_AsDictionary( Dictionary_Entry_Value_GetElement( sysLst, sys_i ) );
			assert( sysDict );

			/* Read contents. */
			meshName = Dictionary_GetString( sysDict, (Dictionary_Entry_Key)"mesh"  );
			remesherName = Dictionary_GetString( sysDict, (Dictionary_Entry_Key)"remesher"  );

			if( strcmp( remesherName, "" ) )
				sys->remesher = Stg_ComponentFactory_ConstructByName( uwCtx->CF, (Name)remesherName, Remesher, True, data  );
			name = Dictionary_GetString( sysDict, (Dictionary_Entry_Key)"displacementField" );

			if(strcmp(name, ""))
			    sys->dispField = Stg_ComponentFactory_ConstructByName( uwCtx->CF, (Name)name, FeVariable, True, data  );
			else
			    sys->dispField = NULL;

			velFieldName = Dictionary_GetString( sysDict, (Dictionary_Entry_Key)"VelocityField"  );
			sys->interval = Dictionary_GetInt_WithDefault( sysDict, (Dictionary_Entry_Key)"interval", -1  );
			sys->wrapTop = Dictionary_GetBool_WithDefault( sysDict, (Dictionary_Entry_Key)"wrapTop", False  );
			sys->wrapBottom = Dictionary_GetBool_WithDefault( sysDict, (Dictionary_Entry_Key)"wrapBottom", False  );
			sys->wrapLeft = Dictionary_GetBool_WithDefault( sysDict, (Dictionary_Entry_Key)"wrapLeft", False  );
			sys->mesh = Stg_ComponentFactory_ConstructByName( uwCtx->CF, (Name)meshName, Mesh, True, data  );
			/* This line is currently not working, have to manually set the velocity field name.
				This should be fixed once this plugin has been converted to a component. */
			/*sys->velField = Stg_ComponentFactory_ConstructByName( uwCtx->CF, (Name)velFieldName, FieldVariable, True, data  );*/
			sys->velField = Stg_ComponentFactory_ConstructByName( uwCtx->CF, (Name)"VelocityField", FieldVariable, True, data  );

			sys->staticTop = Dictionary_GetBool_WithDefault( sysDict, (Dictionary_Entry_Key)"staticTop", False  );
			sys->staticBottom = Dictionary_GetBool_WithDefault( sysDict, (Dictionary_Entry_Key)"staticBottom", False  );
			sys->staticLeft = Dictionary_GetBool_WithDefault( sysDict, (Dictionary_Entry_Key)"staticLeft", False  );
			sys->staticRight = Dictionary_GetBool_WithDefault( sysDict, (Dictionary_Entry_Key)"staticRight", False  );
			sys->staticFront = Dictionary_GetBool_WithDefault( sysDict, (Dictionary_Entry_Key)"staticFront", False  );
			sys->staticBack = Dictionary_GetBool_WithDefault( sysDict, (Dictionary_Entry_Key)"staticBack", False  );

			sys->staticLeftTop = Dictionary_GetBool_WithDefault( sysDict, "staticLeftTop", (sys->staticLeft && sys->staticTop)
                                                                             ? True : False);
			sys->staticRightTop = Dictionary_GetBool_WithDefault( sysDict, "staticRightTop", (sys->staticRight && sys->staticTop)
                                                                             ? True : False );
			sys->staticLeftTopFront = Dictionary_GetBool_WithDefault( sysDict, "staticLeftTopFront",
                                                                                  (sys->staticLeft && sys->staticTop && sys->staticFront)
                                                                             ? True : False );
			sys->staticRightTopFront = Dictionary_GetBool_WithDefault( sysDict, "staticRightTopFront",
                                                                                   (sys->staticRight && sys->staticTop && sys->staticFront)
                                                                             ? True : False );
			sys->staticLeftTopBack = Dictionary_GetBool_WithDefault( sysDict, "staticLeftTopBack",
                                                                                  (sys->staticLeft && sys->staticTop && sys->staticBack)
                                                                             ? True : False );
			sys->staticRightTopBack = Dictionary_GetBool_WithDefault( sysDict, "staticRightTopBack",
                                                                                   (sys->staticRight && sys->staticTop && sys->staticBack)
                                                                             ? True : False );

			sys->staticLeftBottom = Dictionary_GetBool_WithDefault( sysDict, "staticLeftBottom", (sys->staticLeft && sys->staticBottom)
                                                                             ? True : False );
			sys->staticRightBottom = Dictionary_GetBool_WithDefault( sysDict, "staticRightBottom", (sys->staticRight && sys->staticBottom)
                                                                             ? True : False );
			sys->staticLeftBottomFront = Dictionary_GetBool_WithDefault( sysDict, "staticLeftBottomFront",
                                                                                  (sys->staticLeft && sys->staticBottom && sys->staticFront)
                                                                             ? True : False );
			sys->staticRightBottomFront = Dictionary_GetBool_WithDefault( sysDict, "staticRightBottomFront",
                                                                                   (sys->staticRight && sys->staticBottom && sys->staticFront)
                                                                             ? True : False );
			sys->staticLeftBottomBack = Dictionary_GetBool_WithDefault( sysDict, "staticLeftBottomBack",
                                                                                  (sys->staticLeft && sys->staticBottom && sys->staticBack)
                                                                             ? True : False );
			sys->staticRightBottomBack = Dictionary_GetBool_WithDefault( sysDict, "staticRightBottomBack",
                                                                                   (sys->staticRight && sys->staticBottom && sys->staticBack)
                                                                             ? True : False );

			sys->staticLeftFront = Dictionary_GetBool_WithDefault( sysDict, "staticLeftFront", (sys->staticLeft && sys->staticFront)
                                                                             ? True : False );
			sys->staticRightFront = Dictionary_GetBool_WithDefault( sysDict, "staticRightFront", (sys->staticRight && sys->staticFront)
                                                                             ? True : False );
			sys->staticLeftBack = Dictionary_GetBool_WithDefault( sysDict, "staticLeftBack", (sys->staticLeft && sys->staticBack)
                                                                             ? True : False );
			sys->staticRightBack = Dictionary_GetBool_WithDefault( sysDict, "staticRightBack", (sys->staticRight && sys->staticBack)
                                                                             ? True : False );

			sys->staticTopFront = Dictionary_GetBool_WithDefault( sysDict, "staticTopFront", (sys->staticTop && sys->staticFront)
                                                                             ? True : False );
			sys->staticBottomFront = Dictionary_GetBool_WithDefault( sysDict, "staticBottomFront", (sys->staticBottom && sys->staticFront)
                                                                             ? True : False );
			sys->staticTopBack = Dictionary_GetBool_WithDefault( sysDict, "staticTopBack", (sys->staticTop && sys->staticBack)
                                                                             ? True : False );
			sys->staticBottomBack = Dictionary_GetBool_WithDefault( sysDict, "staticBottomBack", (sys->staticBottom && sys->staticBack)
                                                                             ? True : False );

			sys->floatLeftTop = Dictionary_GetBool_WithDefault( sysDict, "floatLeftTop", False );
			sys->floatRightTop = Dictionary_GetBool_WithDefault( sysDict, "floatRightTop", False );

			sys->staticSides = 
                          (sys->staticLeft
                           || sys->staticRight
                           || sys->staticTop
                           || sys->staticBottom
                           || sys->staticFront
                           || sys->staticBack
                           || sys->staticLeftTop
                           || sys->staticRightTop
                           || sys->staticLeftTopFront
                           || sys->staticRightTopFront
                           || sys->staticLeftTopBack
                           || sys->staticRightTopBack
                           || sys->staticLeftBottom
                           || sys->staticRightBottom
                           || sys->staticLeftBottomFront
                           || sys->staticRightBottomFront
                           || sys->staticLeftBottomBack
                           || sys->staticRightBottomBack
                           || sys->staticLeftFront
                           || sys->staticRightFront
                           || sys->staticLeftBack
                           || sys->staticRightBack
                           || sys->staticTopFront
                           || sys->staticBottomFront
                           || sys->staticTopBack
                           || sys->staticBottomBack)
                          ? True : False;


                        if(sys->staticRight && sys->wrapTop
                           && !sys->staticRightTop)
                          sys->x_right_coord =
                            Dictionary_GetDouble( uwCtx->dictionary, "maxX");

                        if(sys->staticLeft && sys->wrapTop
                           && !sys->staticLeftTop)
                          sys->x_left_coord =
                            Dictionary_GetDouble( uwCtx->dictionary, "minX");
                          
			/* Read the list of variables to interpolate. */
			varLst = Dictionary_Entry_Value_GetMember( Dictionary_Entry_Value_GetElement( sysLst, sys_i  ), "fields" );

			if( varLst ) {
				unsigned	var_i;

				sys->nFields = Dictionary_Entry_Value_GetCount( varLst );
				sys->fields = Memory_Alloc_Array( FieldVariable*, sys->nFields, "EulerDeform->systems[].fields" );
				sys->vars = Memory_Alloc_Array( Variable*, sys->nFields, "EulerDeform->systemsp[].vars" );

				for( var_i = 0; var_i <sys->nFields; var_i++ ) {
					Dictionary*	varDict;
					char*			varName;

					/* Get the dictionary for this field tuple. */
					varDict = Dictionary_Entry_Value_AsDictionary( Dictionary_Entry_Value_GetElement( varLst, var_i ) );
					assert( varDict );

					/* Get the field and its variable. */
					varName = Dictionary_GetString( varDict, (Dictionary_Entry_Key)"field"  );
					sys->fields[var_i] = Stg_ComponentFactory_ConstructByName( uwCtx->CF, (Name)varName, FieldVariable, True, data  ); 
					varName = Dictionary_GetString( varDict, (Dictionary_Entry_Key)"variable"  );
					sys->vars[var_i] = Stg_ComponentFactory_ConstructByName( uwCtx->CF, (Name)varName, Variable, True, data ); 
				}
			}
		}
	}

	for( sys_i = 0; sys_i < edCtx->nSystems; sys_i++  ) {
		EulerDeform_System*	sys = edCtx->systems + sys_i;

		/* Create a time integrand for the mesh's coordinates. */
		crdVar = EulerDeform_RegisterLocalNodeCoordsAsVariables( sys, uwCtx->variable_Register, NULL );
		Stg_Component_Build( crdVar, data, False );

		tiData[0] = (Stg_Component*)sys->velField;
		tiData[1] = (Stg_Component*)&sys->mesh->verts;
		crdAdvector = TimeIntegrand_New( "EulerDeform_Velocity", (DomainContext*)uwCtx, edCtx->timeIntegrator, crdVar, 2, tiData, True
			 /* Presume we need to allow fallback on edges of stretching mesh - PatrickSunter, 7 June 2006 */ );
		crdAdvector->_calculateTimeDeriv = EulerDeform_TimeDeriv;

		/* Add to live component register... */
		LiveComponentRegister_Add( uwCtx->CF->LCRegister, (Stg_Component*)crdAdvector );
		Stg_Component_Build( crdAdvector, data, False );
	}

	if( edCtx->nSystems > 0 ) {
		/* Insert the sync step. */
          TimeIntegrator_PrependSetupEP( edCtx->timeIntegrator, "EulerDeform_IntegrationSetup", (void*)EulerDeform_IntegrationSetup, "EulerDeform", edCtx );
	}

	/* Insert the remesh step. Note that this should look for the surface process
	   plugin's time integrator finish routine and ensure we enter the remesh step
	   after that one but before the particle updating routines. */
	TimeIntegrator_PrependFinishEP( edCtx->timeIntegrator, "EulerDeform_Execute", (void*)EulerDeform_Remesh, "EulerDeform", edCtx );
}


void _Underworld_EulerDeform_Destroy( void* component, void* data ) {
	Codelet*					ed	= (Codelet*)component;
	UnderworldContext*	uwCtx = (UnderworldContext*)ed->context;

	assert( component );
	assert( uwCtx );

	/* Clear the lot. */
	/* TODO */
}

/* This creates a set which makes sure not to include the corners if
   we have not made both sides of the corner static. */

IndexSet* EulerDeform_CreateStaticSet(EulerDeform_System* sys)
{
  Grid*		grid;
  unsigned	nNodes;
  unsigned	n_i;
  IndexSet	*set;
  IJK			ijk;

  grid = *(Grid**)ExtensionManager_Get ( sys->mesh->info, sys->mesh, ExtensionManager_GetHandle( sys->mesh->info, (Name)"vertexGrid" )  );

  nNodes = Mesh_GetDomainSize( sys->mesh, MT_VERTEX );
  set = IndexSet_New( nNodes );

  for( n_i = 0; n_i < nNodes; n_i++ ) {
    Bool add;
    add=False;
    RegularMeshUtils_Node_1DTo3D
      ( sys->mesh, Mesh_DomainToGlobal( sys->mesh, MT_VERTEX, n_i ), ijk );

    /* 2D */
    if(sys->mesh->topo->nDims == 2)
      {
        /* Left side and Corner */
        if(ijk[0]==0)
          {
            if(ijk[1]==0)
              {
                if(sys->staticLeftBottom)
                  add=True;
              }
            else if(ijk[1]==grid->sizes[1]-1)
              {
                if(sys->staticLeftTop)
                  add=True;
              }
            else if(sys->staticLeft)
              add=True;
          }
        /* Right side and corner */
        else if(ijk[0]==grid->sizes[0]-1)
          {
            if(ijk[1]==0)
              {
                if(sys->staticRightBottom)
                  add=True;
              }
            else if(ijk[1]==grid->sizes[1]-1)
              {
                if(sys->staticRightTop)
                  add=True;
              }
            else if(sys->staticRight)
              add=True;
          }
        /* Top and Bottom */
        else if((ijk[1]==0 && sys->staticBottom)
                || (ijk[1]==grid->sizes[1]-1 && sys->staticTop))
          add=True;
      }
    /* 3D */
    else if(sys->mesh->topo->nDims == 3)
      {
        /* Left side */
        if(ijk[0]==0)
          {
            /* Left Bottom */
            if(ijk[1]==0)
              {
                if(ijk[2]==0)
                  {
                    if(sys->staticLeftBottomBack)
                      add=True;
                  }
                else if(ijk[2]==grid->sizes[2]-1)
                  {
                    if(sys->staticLeftBottomFront)
                      add=True;
                  }
                else if(sys->staticLeftBottom)
                  add=True;
              }
            /* Left Top */
            else if(ijk[1]==grid->sizes[1]-1)
              {
                if(ijk[2]==0)
                  {
                    if(sys->staticLeftTopBack)
                      add=True;
                  }
                else if(ijk[2]==grid->sizes[2]-1)
                  {
                    if(sys->staticLeftTopFront)
                      add=True;
                  }
                else if(sys->staticLeftTop)
                  add=True;
              }
            /* Left Back */
            else if(ijk[2]==0)
              {
                if(sys->staticLeftBack)
                  add=True;
              }
            /* Left Front */
            else if(ijk[2]==grid->sizes[2]-1)
              {
                if(sys->staticLeftFront)
                  add=True;
              }
            /* Left */
            else if(sys->staticLeft)
              add=True;
          }
        /* Right side */
        else if(ijk[0]==grid->sizes[0]-1)
          {
            /* Right Bottom */
            if(ijk[1]==0)
              {
                if(ijk[2]==0)
                  {
                    if(sys->staticRightBottomBack)
                      add=True;
                  }
                else if(ijk[2]==grid->sizes[2]-1)
                  {
                    if(sys->staticRightBottomFront)
                      add=True;
                  }
                else if(sys->staticRightBottom)
                  add=True;
              }
            /* Right Top */
            else if(ijk[1]==grid->sizes[1]-1)
              {
                if(ijk[2]==0)
                  {
                    if(sys->staticRightTopBack)
                      add=True;
                  }
                else if(ijk[2]==grid->sizes[2]-1)
                  {
                    if(sys->staticRightTopFront)
                      add=True;
                  }
                else if(sys->staticRightTop)
                  add=True;
              }
            /* Right Back */
            else if(ijk[2]==0)
              {
                if(sys->staticRightBack)
                  add=True;
              }
            /* Right Front */
            else if(ijk[2]==grid->sizes[2]-1)
              {
                if(sys->staticRightFront)
                  add=True;
              }
            /* Right */
            else if(sys->staticRight)
              add=True;
          }
        /* Bottom */
        else if(ijk[1]==0)
          {
            if(ijk[2]==0)
              {
                if(sys->staticBottomBack)
                  add=True;
              }
            else if(ijk[2]==grid->sizes[2]-1)
              {
                if(sys->staticBottomFront)
                  add=True;
              }
            else if(sys->staticBottom)
              add=True;
          }
        /* Top */
        else if(ijk[1]==grid->sizes[1]-1)
          {
            if(ijk[2]==0)
              {
                if(sys->staticTopBack)
                  add=True;
              }
            else if(ijk[2]==grid->sizes[2]-1)
              {
                if(sys->staticTopFront)
                  add=True;
              }
            else if(sys->staticTop)
              add=True;
          }
        /*  Front and Back */
        else if((ijk[2]==0 && sys->staticBack)
                || (ijk[2]==grid->sizes[2]-1 && sys->staticFront))
          add=True;
      }

    if( add )
      IndexSet_Add( set, n_i );
  }
  return set;
}

Variable* EulerDeform_RegisterLocalNodeCoordsAsVariables( EulerDeform_System* sys, void* _variable_Register, Variable** variableList ) {
	FeMesh*					self = (FeMesh*)sys->mesh;
	Variable_Register*	variable_Register = (Variable_Register*) _variable_Register;
	Variable*				variable;
	char*						variableName;
	char*						variableNameX;
	char*						variableNameY;
	char*						variableNameZ;

	/* Allocate advection array. */
	sys->verts = AllocArray( double, Mesh_GetLocalSize( self, MT_VERTEX ) * Mesh_GetDimSize( self ) );
	
	/* Append Extension onto names */
	variableName  = Memory_Alloc_Array( char, strlen( self->name ) + strlen( "NodeCoords" ) + 1, "variableName" );
	sprintf( variableName , "%sNodeCoords", self->name );
	
	variableNameX = Memory_Alloc_Array( char, strlen( self->name ) + strlen( "NodeCoordX" ) + 1, "variableNameX" );
	sprintf( variableNameX, "%sNodeCoordX", self->name );

	variableNameY = Memory_Alloc_Array( char, strlen( self->name ) + strlen( "NodeCoordY" ) + 1, "variableNameY" );
	sprintf( variableNameY, "%sNodeCoordY", self->name );

	variableNameZ = Memory_Alloc_Array( char, strlen( self->name ) + strlen( "NodeCoordZ" ) + 1, "variableNameZ" );
	sprintf( variableNameZ, "%sNodeCoordZ", self->name );
	
	/* Construct */
	variable = Variable_NewVector( 
		variableName, 
		self->context,	
		Variable_DataType_Double, 
		Mesh_GetDimSize( self ), 
		(unsigned*)&((IGraph*)self->topo)->remotes[MT_VERTEX]->decomp->locals->size, 
		NULL,
		(void**)&sys->verts, 
		variable_Register, 
		variableNameX,
		variableNameY,
		variableNameZ );

	if ( variableList != NULL ) {
		variableList[ I_AXIS ] = Variable_Register_GetByName( variable_Register, variableNameX );
		variableList[ J_AXIS ] = Variable_Register_GetByName( variable_Register, variableNameY );
		variableList[ K_AXIS ] = Variable_Register_GetByName( variable_Register, variableNameZ );
	}

	/* Clean Up */
	Memory_Free( variableNameZ );
	Memory_Free( variableNameY );
	Memory_Free( variableNameX );
	Memory_Free( variableName );

	return variable;
}


void EulerDeform_IntegrationSetup( void* _timeIntegrator, void* context ) {
	EulerDeform_Context*	edCtx = (EulerDeform_Context*)context;
	unsigned					sys_i;

	FeVariable_SyncShadowValues( edCtx->systems[0].velField );

	/* 
	** We'll need to store side values that we require to be static here, for later
	** return to the mesh.
	*/

	for( sys_i = 0; sys_i < edCtx->nSystems; sys_i++ ) {
		EulerDeform_System* sys = edCtx->systems + sys_i;

		if( sys->staticSides ) {
			IndexSet	*tmpIndSet;
			unsigned	nInds, *inds;
			unsigned	nDims;
			unsigned	ind_i;

			/* Collect indices of all the sides. */

                        tmpIndSet = EulerDeform_CreateStaticSet(sys);
                        IndexSet_GetMembers( tmpIndSet, &nInds, &inds );

			/* Copy coords to temporary array. */
			nDims = Mesh_GetDimSize( sys->mesh );
			sys->sideCoords = AllocArray2D( double, nInds, nDims );
			for( ind_i = 0; ind_i < nInds; ind_i++ )
				memcpy( sys->sideCoords[ind_i], sys->mesh->verts[inds[ind_i]], nDims * sizeof(double) );
                        FreeObject( tmpIndSet );
                        FreeArray( inds );
		}
	}

	/* Update advection arrays. */
	for( sys_i = 0; sys_i < edCtx->nSystems; sys_i++ ) {
		EulerDeform_System*	sys = edCtx->systems + sys_i;
		unsigned					nDims;
		unsigned					nLocalNodes;
		unsigned					n_i;

		nDims = Mesh_GetDimSize( sys->mesh );
		nLocalNodes = Mesh_GetLocalSize( sys->mesh, MT_VERTEX );

		for( n_i = 0; n_i < nLocalNodes; n_i++ )
			memcpy( sys->verts + n_i * nDims, sys->mesh->verts[n_i], nDims * sizeof(double) );
	}
}


Bool EulerDeform_TimeDeriv( void* crdAdvector, Index arrayInd, double* timeDeriv ) {
	TimeIntegrand*		self = (TimeIntegrand*)crdAdvector;
	FeVariable*				velocityField = (FeVariable*)self->data[0];
	InterpolationResult	result = LOCAL;

	/* check if the node information is on the local proc */
	if (arrayInd >= Mesh_GetDomainSize(velocityField->feMesh,MT_VERTEX) )
		result = OTHER_PROC;

	FeVariable_GetValueAtNode( velocityField, arrayInd, timeDeriv );

	/* Check if periodic */
	if ( Stg_Class_IsInstance( velocityField->feMesh->generator, CartesianGenerator_Type ) ) {
		CartesianGenerator* cartesianGenerator = (CartesianGenerator*) velocityField->feMesh->generator;
		if ( cartesianGenerator->periodic[ I_AXIS ] )
			timeDeriv[I_AXIS] = 0.0;
		if ( cartesianGenerator->periodic[ J_AXIS ] )
			timeDeriv[J_AXIS] = 0.0;
		if ( cartesianGenerator->periodic[ K_AXIS ] )
			timeDeriv[K_AXIS] = 0.0;
	}
		
	if ( result == OTHER_PROC || result == OUTSIDE_GLOBAL || isinf(timeDeriv[0]) || isinf(timeDeriv[1]) || 
	     ( velocityField->dim == 3 && isinf(timeDeriv[2]) ) ) 
	{
#if 0
		Journal_Printf( Journal_Register( Error_Type, (Name)self->type  ),
				"Error in func '%s' for particle with index %u.\n\tPosition (%g, %g, %g)\n\tVelocity here is (%g, %g, %g)."
				"\n\tInterpolation result is %s.\n",
				__func__, array_I, coord[0], coord[1], coord[2], 

				InterpolationResultToStringMap[result]  );
#endif	
		return False;	
	}

	return True;
}

Bool _EulerDeform_LineInterp(double** crds,const double* pnt,
                             unsigned fromDim,unsigned toDim, 
			     double* val);
Bool _EulerDeform_QuadYInterp(double** crds,const double* pnt,double* val);

/* Make the top left or right corners float up to the height of the
   point just inside. */

void EulerDeform_FloatRightTop(EulerDeform_System* sys, Grid *grid,
                               double** crds)
{
  int i, max_z;
  IJK ijk, inside_ijk;
  unsigned ind, nLocalNodes, inside;

  nLocalNodes = Mesh_GetLocalSize( sys->mesh, MT_VERTEX );
  ijk[0]=grid->sizes[0]-1;
  ijk[1]=grid->sizes[1]-1;

  if(Mesh_GetDimSize(sys->mesh)==2)
    {
      max_z=1;
    }
  else
    {
      max_z=grid->sizes[2];
    }

  for(i=0;i<max_z;++i)
    {
      ijk[2]=i;
      ind=Grid_Project(grid,ijk);
      if( !(!Mesh_GlobalToDomain( sys->mesh, MT_VERTEX, ind, &ind )
            || ind >= nLocalNodes ))
        {
          inside_ijk[0]=ijk[0]-1;
          inside_ijk[1]=ijk[1];
          inside_ijk[2]=ijk[2];
          inside=Grid_Project(grid,inside_ijk);
          if(!(!Mesh_GlobalToDomain( sys->mesh, MT_VERTEX, inside, &inside )
               || inside >= nLocalNodes ))
            {
              crds[ind][1]=crds[inside][1];
            }
        }
    }
}

void EulerDeform_FloatLeftTop(EulerDeform_System* sys, Grid *grid,
                              double** crds)
{
  int i, max_z;
  IJK ijk, inside_ijk;
  unsigned ind, nLocalNodes, inside;

  nLocalNodes = Mesh_GetLocalSize( sys->mesh, MT_VERTEX );
  ijk[0]=0;
  ijk[1]=grid->sizes[1]-1;

  if(Mesh_GetDimSize(sys->mesh)==2)
    {
      max_z=1;
    }
  else
    {
      max_z=grid->sizes[2];
    }

  for(i=0;i<max_z;++i)
    {
      ijk[2]=i;
      ind=Grid_Project(grid,ijk);
      if( !(!Mesh_GlobalToDomain( sys->mesh, MT_VERTEX, ind, &ind )
            || ind >= nLocalNodes ))
        {
          inside_ijk[0]=ijk[0]+1;
          inside_ijk[1]=ijk[1];
          inside_ijk[2]=ijk[2];
          inside=Grid_Project(grid,inside_ijk);
          if(!(!Mesh_GlobalToDomain( sys->mesh, MT_VERTEX, inside, &inside )
               || inside >= nLocalNodes ))
            crds[ind][1]=crds[inside][1];
        }
    }
}

/* Remesh the left or right top corners in 2D or 3D */

void EulerDeform_Remesh_Corner(Mesh *mesh, int corner, int inside,
                               double side_coord) {
  IJK		ijk;
  Grid *grid;
  grid =
    *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
                                   ExtensionManager_GetHandle( mesh->info,
                                                               "vertexGrid" ) );
  ijk[0]=corner;
  ijk[1]=grid->sizes[1]-1;
  if(Mesh_GetDimSize(mesh)==2)
    {
      unsigned n_corner, n_interior, n, n_in;
      n_corner=RegularMeshUtils_Node_3DTo1D(mesh,ijk);
      ijk[0]=inside;
      n_interior=RegularMeshUtils_Node_3DTo1D(mesh,ijk);
      if(Mesh_GlobalToDomain(mesh,MT_VERTEX,
                             n_corner,&n)
         && Mesh_GlobalToDomain(mesh,MT_VERTEX,
                                n_interior,&n_in))
        {
          double *crds[2];

          crds[0]=mesh->verts[n];
          crds[1]=mesh->verts[n_in];
          if(!_EulerDeform_LineInterp(crds,&side_coord,0,1,
                                       &(mesh->verts[n][1])))
            {
              printf("The side is moving in the wrong direction.\n");
              printf("%g %g %g %g %g\n",side_coord,crds[0][0],crds[0][1],crds[1][0],crds[1][1]);
              abort();
            }
          mesh->verts[n][0]=side_coord;
        }
    }
  else /* 3D */
    {
      for(ijk[2]=0; ijk[2]<grid->sizes[2]; ++ijk[2])
        {
          unsigned n_corner, n_interior, n, n_in;
          ijk[0]=corner;
          n_corner=RegularMeshUtils_Node_3DTo1D(mesh,ijk);
          ijk[0]=inside;
          n_interior=RegularMeshUtils_Node_3DTo1D(mesh,ijk);
          if(Mesh_GlobalToDomain(mesh,MT_VERTEX,n_corner,&n)
             && Mesh_GlobalToDomain(mesh,MT_VERTEX,n_interior,&n_in))
            {
              double *crds[2];
              crds[0]=mesh->verts[n];
              crds[1]=mesh->verts[n_in];

              if(!_EulerDeform_LineInterp(crds,&side_coord,0,1,&(mesh->verts[n][1]))
                 || !_EulerDeform_LineInterp(crds,&side_coord,0,2,&(mesh->verts[n][2])))
                {
                  printf("The side is moving in the wrong direction.\n");
                  printf("%g %g %g %g %g\n",side_coord,crds[0][0],crds[0][1],crds[1][0],crds[1][1]);
                  abort();
                }
              mesh->verts[n][0]=side_coord;
            }
        }
    }
}

void EulerDeform_WrapSurface( EulerDeform_System* sys, double** oldCrds, int top );

void EulerDeform_Remesh( TimeIntegrand* crdAdvector, EulerDeform_Context* edCtx ) {
  Mesh_Algorithms	*tmpAlgs, *oldAlgs;
  unsigned	sys_i;

  assert( edCtx );

  /* We do the second system first, because that is the vertex
     centered system.  We want the cell centered and vertex
     centered systems to be compatible, so the cell centered
     system is based on the vertex centered one. */

  for( sys_i = 1; sys_i < edCtx->nSystems+1; sys_i++ ) {
    EulerDeform_System*	sys = edCtx->systems
      + sys_i%edCtx->nSystems;
    double**		oldCrds;
    double**		newCrds;
    unsigned		nDomainNodes;
    unsigned		nDims;
    unsigned		var_i, n_i, dof_i;
    Grid *grid;
    grid =
      *(Grid**)ExtensionManager_Get(sys->mesh->info, sys->mesh, 
                                    ExtensionManager_GetHandle( sys->mesh->info,
                                                                "vertexGrid" ));
    /* Update the displacement field. */
    if(sys->dispField) {
      double disp[3];
      int num_verts, num_dims;
      int ii, jj;

      num_dims = Mesh_GetDimSize(sys->mesh);
      num_verts = Mesh_GetLocalSize(sys->mesh, MT_VERTEX);
      for(ii = 0; ii < num_verts; ii++) {
        FeVariable_GetValueAtNode(sys->dispField, ii, disp);
        for(jj = 0; jj < num_dims; jj++)
          disp[jj] += sys->verts[ii*num_dims + jj] - sys->mesh->verts[ii][jj];
        FeVariable_SetValueAtNode(sys->dispField, ii, disp);
      }
    }

    nDims = Mesh_GetDimSize( sys->mesh );
    
    /* Update all local coordinates. */
    for( n_i = 0; n_i < Mesh_GetLocalSize( sys->mesh, MT_VERTEX ); n_i++ )
      memcpy( sys->mesh->verts[n_i], sys->verts + n_i * nDims, nDims * sizeof(double) );
    
    /* Revert side coordinates if required. */
    if( sys->staticSides ) {
      IndexSet	*tmpIndSet;
      unsigned	nInds, *inds;
      unsigned	ind_i;
      
      /* Collect indices of all the sides. */
      
      tmpIndSet = EulerDeform_CreateStaticSet(sys);
      IndexSet_GetMembers( tmpIndSet, &nInds, &inds );
      
      /* Copy back coords. */
      for( ind_i = 0; ind_i < nInds; ind_i++ )
        memcpy( sys->mesh->verts[inds[ind_i]], sys->sideCoords[ind_i], nDims * sizeof(double) );
      FreeObject( tmpIndSet );
      FreeArray( sys->sideCoords );
      
      if(sys->wrapTop)
        {
          if(sys->staticLeft && !sys->staticLeftTop && !sys->floatLeftTop)
            {
              EulerDeform_Remesh_Corner(sys->mesh,0,1,sys->x_left_coord);
            }
          if(sys->staticRight && !sys->staticRightTop && !sys->floatRightTop)
            {
              EulerDeform_Remesh_Corner(sys->mesh,grid->sizes[0]-1,
                                        grid->sizes[0]-2,
                                        sys->x_right_coord);
            }
        }

    }

    /* If we have regular mesh algorithms specified, set the
       algorithms temporarily to an irregular method. */
    if( !strcmp( sys->mesh->algorithms->type, "Mesh_RegularAlgorithms" ) && sys->remesher ) {
      tmpAlgs = Mesh_Algorithms_New( "", edCtx->ctx );
      oldAlgs = sys->mesh->algorithms;
      sys->mesh->algorithms = NULL;
      Mesh_SetAlgorithms( sys->mesh, tmpAlgs );
    }
    else
      tmpAlgs = NULL;

    /* Every system should synchronise the mesh coordinates. */
    Mesh_Sync( sys->mesh );
    Mesh_DeformationUpdate( sys->mesh );

    /* Only if remesher specified. */
    if( !sys->remesher ) {
      continue;
    }

    /* If a remesh interval is requested, check now. */
    if( sys->interval > 0 && edCtx->ctx->timeStep % sys->interval > 0 ) {
      Journal_Printf( Underworld_Info,
                      "*** EulerDeform: Not remeshing this timestep.\n" );
      continue;
    }
    Journal_Printf( Underworld_Info, "*** EulerDeform: Remeshing.\n" );

    /* Store old coordinates. */
    nDomainNodes = FeMesh_GetNodeDomainSize( sys->mesh );
    oldCrds = AllocArray2D( double, nDomainNodes, nDims );
    for( n_i = 0; n_i < nDomainNodes; n_i++ )
      memcpy( oldCrds[n_i], sys->mesh->verts[n_i], nDims * sizeof(double) );

    /* Remesh the system. */
    Stg_Component_Execute( sys->remesher, NULL, True );
    Mesh_Sync( sys->mesh );

    /* Shrink wrap the top/bottom surface. */
    if( sys->wrapTop )
      EulerDeform_WrapSurface( sys, oldCrds, 1 );
    if( sys->wrapBottom )
      EulerDeform_WrapSurface( sys, oldCrds, 0 );

    /* Swap old coordinates back in temporarily. */
    newCrds = sys->mesh->verts;
    sys->mesh->verts = oldCrds;

    /* Interpolate the variables. */
    for( var_i = 0; var_i < sys->nFields; var_i++ )
      EulerDeform_InterpVar( sys->fields[var_i],
                             NULL/*sys->vars[var_i]*/, sys->mesh, newCrds );

    /* Float the top left and right corners if needed.  We do this
       after interpolating, because these points almost certainly are
       outside of the domain, and so can not be interpolated to. */
    if(sys->floatLeftTop)
      EulerDeform_FloatLeftTop(sys,grid,newCrds);
    if(sys->floatRightTop)
      EulerDeform_FloatRightTop(sys,grid,newCrds);

    /* Create an artificial displacement field from the nodal
     * displacements between newCrds and oldCrds.  This displacement
     * is currently used to correct the advDiffEqn's nodal velocity
     * input */
    if( edCtx->artDField ) {
      double artDis[3]; /* temporary displacement vector */

      for( n_i = 0 ; n_i < nDomainNodes; n_i++ ) {
        for( dof_i = 0 ; dof_i < nDims ; dof_i++ ) {
          artDis[dof_i] = newCrds[n_i][dof_i] - oldCrds[n_i][dof_i];
        }
        FeVariable_SetValueAtNode( edCtx->artDField, n_i, artDis );
      }
        
    }

    /* Swap back coordinates and free memory. */
    sys->mesh->verts = newCrds;
    FreeArray( oldCrds );
    
    /* Swap back old algorithms. */
    if( tmpAlgs ) {
      Mesh_SetAlgorithms( sys->mesh, oldAlgs );
    }

    /* Re-sync with new coordinates. */
    Mesh_Sync( sys->mesh );
    Mesh_DeformationUpdate( sys->mesh );
    for( var_i = 0; var_i < sys->nFields; var_i++ )
      FeVariable_SyncShadowValues( sys->fields[var_i] );
  }
}


void EulerDeform_InterpVar( FieldVariable* field, Variable* var, Mesh* mesh, double** newCrds ) {
	double*		newVals;
	unsigned	curValInd = 0;
	unsigned	nLocalNodes;
	unsigned	n_i, c_i;

	assert( field );
	/*assert( var );*/
	assert( newCrds );

	/* Allocate for new values. */
	nLocalNodes = Mesh_GetLocalSize( mesh, MT_VERTEX );
	newVals = Memory_Alloc_Array( double, field->fieldComponentCount * nLocalNodes, "EulerDeform_InterpVar::newVals" );

	/* Interpolate using new node coordinates. */
	for( n_i = 0; n_i < nLocalNodes; n_i++ ) {
		InterpolationResult	res;

		/* Interpolate the value. */
		res = FieldVariable_InterpolateValueAt( field, newCrds[n_i], newVals + n_i * field->fieldComponentCount );
		if( res == OTHER_PROC || res == OUTSIDE_GLOBAL ) {
			FeVariable_GetValueAtNode( (FeVariable*)field, n_i, 
						   newVals + n_i * field->fieldComponentCount );
		}
	}

	/* Transfer the new values back to the variable. */
	for( n_i = 0; n_i < nLocalNodes; n_i++ ) {
	   for( c_i = 0; c_i < field->fieldComponentCount; c_i++ )
	      DofLayout_SetValueDouble( ((FeVariable*)field)->dofLayout, n_i, c_i, newVals[curValInd++] );
	}

	/* Free the values array. */
	FreeArray( newVals );
}

void EulerDeform_InternalLoop( EulerDeform_System* sys, Grid* grm, double** oldCrds, unsigned* ijk, unsigned curDim, int top );

void EulerDeform_WrapSurface( EulerDeform_System* sys, double** oldCrds, int top ) {
	IJK	ijk;
	Grid*	grm;
	Mesh*	mesh;

	assert( sys );
	assert( oldCrds );

	/* Loop over top internal surface. */
	mesh = sys->mesh;
	grm = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
					     ExtensionManager_GetHandle( mesh->info, (Name)"vertexGrid" )  );
	EulerDeform_InternalLoop( sys, grm, oldCrds, ijk, 0, top );
}

#if 0
void EulerDeform_WrapLeftSurface( EulerDeform_System* sys, double** oldCrds ) {
	IJK	ijk;
	GRM	grm;

	assert( sys );
	assert( oldCrds );

	/* Loop over top internal surface. */
	RegMesh_Generalise( sys->mesh, &grm );
	EulerDeform_LeftInternalLoop( sys, &grm, oldCrds, ijk, 0 );
}
#endif


void _EulerDeform_TriBarycenter( double** tri, const double* pnt, double* dst ) {
	double	a = tri[0][0] - tri[2][0];
	double	b = tri[1][0] - tri[2][0];
	double	c = tri[2][0] - pnt[0];
	double	d = tri[0][1] - tri[2][1];
	double	e = tri[1][1] - tri[2][1];
	double	f = tri[2][1] - pnt[1];
	double	g = tri[0][2] - tri[2][2];
	double	h = tri[1][2] - tri[2][2];
	double	i = tri[2][2] - pnt[2];

	dst[0] = (b * (f + i) - c * (e + h)) / (a * (e + h) - b * (d + g));
	dst[1] = (a * (f + i) - c * (d + g)) / (b * (d + g) - a * (e + h));
	dst[2] = 1.0 - dst[0] - dst[1];
}


Bool _EulerDeform_QuadYInterp( double** crds, const double* pnt, double* val ) {
	double*		modCrds[4];
	double		modCrds0[2], modCrds1[2], modCrds2[2], modCrds3[2];
	unsigned*	inds[2];
	unsigned	inds0[3], inds1[3];
	unsigned	inc[4];
	double		modPnt[3];
	unsigned	inside;
	double		bc[3];

	modCrds[0] = modCrds0;
	modCrds[1] = modCrds1;
	modCrds[2] = modCrds2;
	modCrds[3] = modCrds3;
	modCrds[0][0] = crds[0][0]; modCrds[0][1] = crds[0][2];
	modCrds[1][0] = crds[1][0]; modCrds[1][1] = crds[1][2];
	modCrds[2][0] = crds[2][0]; modCrds[2][1] = crds[2][2];
	modCrds[3][0] = crds[3][0]; modCrds[3][1] = crds[3][2];
	modPnt[0] = pnt[0]; modPnt[1] = pnt[2];

	inds[0] = inds0;
	inds[1] = inds1;
	inds[0][0] = 0; inds[0][1] = 1; inds[0][2] = 2;
	inds[1][0] = 1; inds[1][1] = 3; inds[1][2] = 2;
	inc[0] = 0; inc[1] = 1; inc[2] = 2; inc[3] = 3;

	if( Simplex_Search2D( modCrds, inc, 2, inds, 
			      modPnt, bc, &inside ) )
	{
		*val = bc[0] * crds[inds[inside][0]][1] + bc[1] * crds[inds[inside][1]][1] + 
			bc[2] * crds[inds[inside][2]][1];
		return True;
	}
	else
		return False;
}


Bool _EulerDeform_FindBarycenter1D( const double* crds, const double pnt, 
				    double* bcs )
{
	assert( crds );
	assert( bcs );

	bcs[1] = (pnt - crds[0])/(crds[1] - crds[0]);
	bcs[0] = 1.0 - bcs[1];

	return (bcs[0] >= 0.0 && bcs[0] <= 1.0 && bcs[1] >= 0.0 && bcs[1] <= 1.0) ? True : False;
}


Bool _EulerDeform_LineInterp(double** crds, const double* pnt, unsigned fromDim, unsigned toDim, 
			      double* val ) {
	double	bcCrds[2];
	double	bcs[2];

	assert( crds );
	assert( val );

	bcCrds[0] = crds[0][fromDim];
	bcCrds[1] = crds[1][fromDim];
	if( _EulerDeform_FindBarycenter1D( bcCrds, pnt[fromDim], bcs ) ) {
		*val = bcs[0]*crds[0][toDim] + bcs[1]*crds[1][toDim];
		return True;
	}

	return False;
}


#if 0
Bool _EulerDeform_QuadZInterp( double** crds, const double* pnt, double* val ) {
	double		modCrds[4][3];
	double		modPnt[3];
	unsigned	inds[3];
	double		bc[3];

	modCrds[0][0] = crds[0][0]; modCrds[0][1] = crds[0][1]; modCrds[0][2] = 0.0;
	modCrds[1][0] = crds[1][0]; modCrds[1][1] = crds[1][1]; modCrds[1][2] = 0.0;
	modCrds[2][0] = crds[2][0]; modCrds[2][1] = crds[2][1]; modCrds[2][2] = 0.0;
	modCrds[3][0] = crds[3][0]; modCrds[3][1] = crds[3][1]; modCrds[3][2] = 0.0;
	modPnt[0] = pnt[0]; modPnt[1] = pnt[1]; modPnt[2] = 0.0;

	if( _HexaEL_FindTriBarycenter( (const double**)modCrds, modPnt, bc, inds, INCLUSIVE_UPPER_BOUNDARY, NULL, 0 ) ) {
		*val = bc[0]*crds[inds[0]][1] + bc[1]*crds[inds[1]][1] + bc[2]*crds[inds[2]][1];
		return True;
	}
	else
		return False;
}
#endif


void EulerDeform_InternalLoop( EulerDeform_System* sys, Grid* grm, double** oldCrds, unsigned* ijk, unsigned curDim, int top ) {
	unsigned	nDims;
	XYZ		newCrd, oldCrd;
	double*		crds[4];
	double		crds0[3], crds1[3], crds2[3], crds3[3];
	unsigned	centerInd;
	Mesh*		mesh;
	unsigned	nLocalNodes;
	unsigned	ind;
        double          fudge_factor;
        
        /* fudge_factor nudges the coordinate inside the mesh a little
           to avoid failed interpolations */
        fudge_factor=1.0e-10;

	if( curDim < grm->nDims ) {
		if( curDim == 1 ) {
                        if(top) {
                          ijk[1] = grm->sizes[curDim] - 1;
                        }
                        else {
                          ijk[1]=0;
                        }
			EulerDeform_InternalLoop( sys, grm, oldCrds, ijk, curDim + 1, top );
		}
		else {
			for( ijk[curDim] = 0; ijk[curDim] < grm->sizes[curDim]; ijk[curDim]++ ) {
                          EulerDeform_InternalLoop( sys, grm, oldCrds, ijk, curDim + 1, top );
			}
		}
	}
	else {
		if( grm->nDims == 2 ) {
			mesh = sys->mesh;
			nDims = Mesh_GetDimSize( mesh );

			crds[0] = crds0;
			crds[1] = crds1;
			nLocalNodes = Mesh_GetLocalSize( mesh, MT_VERTEX );

			/* Skip corners. */
			if( ijk[0] == 0 || ijk[0] == grm->sizes[0] - 1 ) {
				return;
			}

			/* Get old and new coordinate. */
			centerInd = Grid_Project( grm, ijk );
			if( !Mesh_GlobalToDomain( mesh, MT_VERTEX, centerInd, &centerInd ) || centerInd >= nLocalNodes )
				return;

			newCrd[0] = mesh->verts[centerInd][0];
			newCrd[1] = mesh->verts[centerInd][1];
			oldCrd[0] = oldCrds[centerInd][0];
			oldCrd[1] = oldCrds[centerInd][1];

			/* Are we left or right? */
			if( newCrd[0] < oldCrd[0] ) {
				ijk[0]--; ind = Grid_Project( grm, ijk ); ijk[0]++;
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
				memcpy( crds[0], oldCrds[ind], nDims * sizeof(double) );
				memcpy( crds[1], oldCrd, nDims * sizeof(double) );
			}
			else if( newCrd[0] > oldCrd[0] ) {
				ijk[0]++; ind = Grid_Project( grm, ijk ); ijk[0]--;
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
				memcpy( crds[1], oldCrds[ind], nDims * sizeof(double) );
				memcpy( crds[0], oldCrd, nDims * sizeof(double) );
			}

                        if(newCrd[0]==oldCrd[0]) {
                          mesh->verts[centerInd][1]=oldCrd[1];
                        }
                        else {

                          /* Interpolate. */
#ifndef NDEBUG
			assert( _EulerDeform_LineInterp(crds, newCrd, 0, 1, &mesh->verts[centerInd][1] ) );
#else
			_EulerDeform_LineInterp(crds, newCrd, 0, 1, &mesh->verts[centerInd][1] );
#endif
                        if((mesh->verts[centerInd][1]>0) ^ top)
                          {
                            mesh->verts[centerInd][1] *= 1+fudge_factor;
                          }
                        else
                          {
                            mesh->verts[centerInd][1] *= 1-fudge_factor;
                          }
                        }
		}
		else if( grm->nDims == 3 ) {
			mesh = sys->mesh;
			nDims = Mesh_GetDimSize( mesh );

			crds[0] = crds0; crds[1] = crds1; crds[2] = crds2; crds[3] = crds3;
			nLocalNodes = Mesh_GetLocalSize( mesh, MT_VERTEX );

			/* Skip corners. */
			if( (ijk[0] == 0 || ijk[0] == grm->sizes[0] - 1) && 
			    (ijk[2] == 0 || ijk[2] == grm->sizes[2] - 1))
			{
				return;
			}

			/* Get old and new coordinate. */
			centerInd = Grid_Project( grm, ijk );
			if( !Mesh_GlobalToDomain( mesh, MT_VERTEX, centerInd, &centerInd ) || centerInd >= nLocalNodes )
				return;

			newCrd[0] = mesh->verts[centerInd][0];
			newCrd[1] = mesh->verts[centerInd][1];
			newCrd[2] = mesh->verts[centerInd][2];
			oldCrd[0] = oldCrds[centerInd][0];
			oldCrd[1] = oldCrds[centerInd][1];
			oldCrd[2] = oldCrds[centerInd][2];

			/* Handle internal nodes. */
			if( ijk[0] > 0 && ijk[2] > 0 ) {
				ijk[0]--; ijk[2]--; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
				memcpy( crds[0], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]++; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
				memcpy( crds[1], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]--; ijk[2]++; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
				memcpy( crds[2], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]++; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
				memcpy( crds[3], oldCrds[ind], nDims * sizeof(double) );

				if( _EulerDeform_QuadYInterp( crds, newCrd, &mesh->verts[centerInd][1] ) ) {
                                  if((mesh->verts[centerInd][1]>0) ^ top)
                                    {
                                      mesh->verts[centerInd][1] *= 1+fudge_factor;
                                    }
                                  else
                                    {
                                      mesh->verts[centerInd][1] *= 1-fudge_factor;
                                    }
                                  return;
				}
			}

			if( ijk[0] > 0 && ijk[2] < grm->sizes[2] - 1 ) {
				ijk[0]--; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
				memcpy( crds[0], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]++; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
				memcpy( crds[1], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]--; ijk[2]++; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
				memcpy( crds[2], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]++; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
				memcpy( crds[3], oldCrds[ind], nDims * sizeof(double) );

				ijk[2]--;
				if( _EulerDeform_QuadYInterp( crds, newCrd, &mesh->verts[centerInd][1] ) ) {
                                  if((mesh->verts[centerInd][1]>0) ^ top)
                                    {
                                      mesh->verts[centerInd][1] *= 1+fudge_factor;
                                    }
                                  else
                                    {
                                      mesh->verts[centerInd][1] *= 1-fudge_factor;
                                    }
                                  return;
				}
			}

			if( ijk[0] < grm->sizes[0] - 1 && ijk[2] > 0 ) {
				ijk[2]--; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
				memcpy( crds[0], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]++; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
				memcpy( crds[1], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]--; ijk[2]++; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
				memcpy( crds[2], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]++; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
				memcpy( crds[3], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]--;
				if( _EulerDeform_QuadYInterp( crds, newCrd, &mesh->verts[centerInd][1] ) ) {
                                  if((mesh->verts[centerInd][1]>0) ^ top)
                                    {
                                      mesh->verts[centerInd][1] *= 1+fudge_factor;
                                    }
                                  else
                                    {
                                      mesh->verts[centerInd][1] *= 1-fudge_factor;
                                    }
                                  return;
				}
			}

			if( ijk[0] < grm->sizes[0] - 1 && ijk[2] < grm->sizes[2] - 1 ) {
				ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
				memcpy( crds[0], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]++; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
				memcpy( crds[1], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]--; ijk[2]++; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
				memcpy( crds[2], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]++; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
				memcpy( crds[3], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]--; ijk[2]--;
				if( _EulerDeform_QuadYInterp( crds, newCrd, &mesh->verts[centerInd][1] ) ) {
                                  if((mesh->verts[centerInd][1]>0) ^ top)
                                    {
                                      mesh->verts[centerInd][1] *= 1+fudge_factor;
                                    }
                                  else
                                    {
                                      mesh->verts[centerInd][1] *= 1-fudge_factor;
                                    }
                                  return;
				}
			}

			assert( 0 );
		}
		else {
			assert( 0 );
		}
	}
}

#if 0
void EulerDeform_LeftInternalLoop( EulerDeform_System* sys, GRM* grm, Coord* oldCrds, unsigned* ijk, unsigned curDim ) {
	if( curDim < grm->nDims ) {
		if( curDim == 0 ) {
			ijk[0] = 0;
			EulerDeform_LeftInternalLoop( sys, grm, oldCrds, ijk, curDim + 1 );
		}
		else {
			for( ijk[curDim] = 1; ijk[curDim] < grm->nNodes[curDim] - 1; ijk[curDim]++ ) {
				EulerDeform_LeftInternalLoop( sys, grm, oldCrds, ijk, curDim + 1 );
			}
		}
	}
	else {
		if( grm->nDims == 2 ) {
			XYZ		newCrd, oldCrd;
			unsigned	centerInd;
			Mesh*		mesh = sys->mesh;

			/* Get old and new coordinate. */
			GRM_Project( grm, ijk, &centerInd );
			newCrd[0] = mesh->nodeCoord[centerInd][0];
			newCrd[1] = mesh->nodeCoord[centerInd][1];
			oldCrd[0] = oldCrds[centerInd][0];
			oldCrd[1] = oldCrds[centerInd][1];

			/* Are we above or below? */
			if( newCrd[1] < oldCrd[1] ) {
				XYZ		leftCrd;
				unsigned	leftInd;
				double		a0, a1;

				/* Get left old coord. */
				ijk[1]--;
				GRM_Project( grm, ijk, &leftInd );
				ijk[1]++;
				leftCrd[0] = oldCrds[leftInd][0];
				leftCrd[1] = oldCrds[leftInd][1];

				/* Calc barycenter. */
				a1 = (newCrd[1] - leftCrd[1]) / (oldCrd[1] - leftCrd[1]);
				a0 = 1.0 - a1;
				mesh->nodeCoord[centerInd][0] = (a0 * leftCrd[0] + a1 * oldCrd[0])*(1.0+1.0e-10);
			}
			else {
				XYZ		rightCrd;
				unsigned	rightInd;
				double		a0, a1;

				/* Get right old coord. */
				ijk[1]++;
				GRM_Project( grm, ijk, &rightInd );
				ijk[1]--;
				rightCrd[0] = oldCrds[rightInd][0];
				rightCrd[1] = oldCrds[rightInd][1];

				/* Calc barycenter. */
				a1 = (newCrd[1] - oldCrd[1]) / (rightCrd[1] - oldCrd[1]);
				a0 = 1.0 - a1;
				mesh->nodeCoord[centerInd][0] = (a0 * oldCrd[0] + a1 * rightCrd[0])*(1.0+1.0e-10);
			}
		}
		else if( grm->nDims == 3 ) {
			assert( 0 );
		}
		else {
			assert( 0 );
		}
	}
}
#endif


