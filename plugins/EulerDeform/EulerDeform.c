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
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include "types.h"
#include "Context.h"
#include "EulerDeform.h"


const char*		EULERDEFORM_PLUGIN_TAG = "EulerDeform";
const Type		Underworld_EulerDeform_Type = "EulerDeform";
ExtensionInfo_Index	EulerDeform_ContextHandle;


Index _Underworld_EulerDeform_Register( PluginsManager* pluginsMgr ) {
	return PluginsManager_Submit( pluginsMgr, 
				      Underworld_EulerDeform_Type, 
				      "0", 
				      _Underworld_EulerDeform_DefaultNew );
}


void* _Underworld_EulerDeform_DefaultNew( Name name ) {
	return _Codelet_New( sizeof(Codelet), 
			     Underworld_EulerDeform_Type, 
			     _Codelet_Delete, 
			     _Codelet_Print, 
			     _Codelet_Copy, 
			     _Underworld_EulerDeform_DefaultNew, 
			     _Underworld_EulerDeform_Construct, 
			     _Underworld_EulerDeform_Build, 
			     _Codelet_Initialise, 
			     _Codelet_Execute, 
			     _Underworld_EulerDeform_Destroy, 
			     name );
}


void _Underworld_EulerDeform_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	UnderworldContext*	uwCtx;
	EulerDeform_Context*	edCtx;
	Dictionary_Entry_Value*	edDict;
	Dictionary_Entry_Value*	sysLst;

	assert( component );
	assert( cf );

	Journal_DPrintf( Underworld_Debug, "In: %s( void* )\n", __func__ );

	/* Retrieve context. */
	uwCtx = (UnderworldContext*)Stg_ComponentFactory_ConstructByName( 
		cf,
		"context", 
		UnderworldContext, 
		True,
		data );

	/* Create new context. */
	EulerDeform_ContextHandle = ExtensionManager_Add( uwCtx->extensionMgr, 
							  Underworld_EulerDeform_Type,
							  sizeof(EulerDeform_Context) );
	edCtx = ExtensionManager_Get( uwCtx->extensionMgr, 
				      uwCtx, 
				      EulerDeform_ContextHandle );
	memset( edCtx, 0, sizeof(EulerDeform_Context) );

	/* Get the dictionary. */
	edDict = Dictionary_Get( uwCtx->dictionary, "EulerDeform" );
	if( !edDict ) {
		return;
	}

	/* Read system list. */
	sysLst = Dictionary_Entry_Value_GetMember( edDict, "systems" );
	if( sysLst ) {
		unsigned	sys_i;

		/* Allocate for systems. */
		edCtx->nSystems = Dictionary_Entry_Value_GetCount( sysLst );
		edCtx->systems = Memory_Alloc_Array( EulerDeform_System, edCtx->nSystems, "EulerDeform->systems" );
		memset( edCtx->systems, 0, sizeof(EulerDeform_System) * edCtx->nSystems );

		for( sys_i = 0; sys_i < edCtx->nSystems; sys_i++ ) {
			EulerDeform_System*	sys = edCtx->systems + sys_i;
			Dictionary*		sysDict;
			Dictionary_Entry_Value*	varLst;
			char*			meshName;
			char*			remesherName;
			char*			velFieldName;

			/* Get the dictionary for this system. */
			sysDict = Dictionary_Entry_Value_AsDictionary( Dictionary_Entry_Value_GetElement( sysLst, sys_i ) );
			assert( sysDict );

			/* Read contents. */
			meshName = Dictionary_GetString( sysDict, "mesh" );
			remesherName = Dictionary_GetString( sysDict, "remesher" );
			if( strcmp( remesherName, "" ) )
				sys->remesher = Stg_ComponentFactory_ConstructByName( cf, remesherName, Remesher, True, data );
			velFieldName = Dictionary_GetString( sysDict, "velocityField" );
			sys->wrapTop = Dictionary_GetBool_WithDefault( sysDict, "wrapTop", False );
			sys->wrapBottom = Dictionary_GetBool_WithDefault( sysDict, "wrapBottom", False );
			sys->wrapLeft = Dictionary_GetBool_WithDefault( sysDict, "wrapLeft", False );
			sys->mesh = Stg_ComponentFactory_ConstructByName( cf, meshName, Mesh, True, data );
			sys->velField = Stg_ComponentFactory_ConstructByName( cf, velFieldName, FieldVariable, True, data );

			sys->staticSides = Dictionary_GetBool_WithDefault( sysDict, "staticSides", False );

			/* Read the list of variables to interpolate. */
			varLst = Dictionary_Entry_Value_GetMember( Dictionary_Entry_Value_GetElement( sysLst, sys_i ), "fields" );
			if( varLst ) {
				unsigned	var_i;

				sys->nFields = Dictionary_Entry_Value_GetCount( varLst );
				sys->fields = Memory_Alloc_Array( FieldVariable*, sys->nFields, "EulerDeform->systems[].fields" );
				sys->vars = Memory_Alloc_Array( Variable*, sys->nFields, "EulerDeform->systemsp[].vars" );
				for( var_i = 0; var_i <sys->nFields; var_i++ ) {
					Dictionary*	varDict;
					char*		varName;

					/* Get the dictionary for this field tuple. */
					varDict = Dictionary_Entry_Value_AsDictionary( 
						Dictionary_Entry_Value_GetElement( varLst, var_i ) );
					assert( varDict );

					/* Get the field and its variable. */
					varName = Dictionary_GetString( varDict, "field" );
					sys->fields[var_i] = Stg_ComponentFactory_ConstructByName( 
						cf, 
						varName, 
						FieldVariable, 
						True, 
						data ); 
					varName = Dictionary_GetString( varDict, "variable" );
					sys->vars[var_i] = Stg_ComponentFactory_ConstructByName( 
						cf, 
						varName, 
						Variable, 
						True, 
						data ); 
				}
			}
		}
	}
}


void _Underworld_EulerDeform_Build( void* component, void* data ) {
	UnderworldContext*	uwCtx = (UnderworldContext*)data;
	EulerDeform_Context*	edCtx;
	Variable*		crdVar;
	TimeIntegratee*		crdAdvector;
	Stg_Component*		tiData[2];
	unsigned		sys_i;

	assert( component );
	assert( uwCtx );

	edCtx = ExtensionManager_Get( uwCtx->extensionMgr, 
				      uwCtx, 
				      EulerDeform_ContextHandle );

	for( sys_i = 0; sys_i < edCtx->nSystems; sys_i++ ) {
		EulerDeform_System*	sys = edCtx->systems + sys_i;

		/* Create a time integratee for the mesh's coordinates. */
		crdVar = EulerDeform_RegisterLocalNodeCoordsAsVariables( sys, uwCtx->variable_Register, NULL );
		Build( crdVar, data, False );

		tiData[0] = (Stg_Component*)sys->velField;
		tiData[1] = (Stg_Component*)&sys->mesh->verts;
		crdAdvector = TimeIntegratee_New( "EulerDeform_Velocity", 
						  uwCtx->timeIntegrator, 
						  crdVar, 
						  2, 
						  tiData,
						  True /* Presume we need to allow fallback on edges of
							  stretching mesh - PatrickSunter, 7 June 2006 */ );
		crdAdvector->_calculateTimeDeriv = EulerDeform_TimeDeriv;

		/* Add to live component register... */
		LiveComponentRegister_Add( uwCtx->CF->LCRegister, (Stg_Component*)crdAdvector );
		Build( crdAdvector, data, False );
	}

	if( edCtx->nSystems > 0 ) {
		/* Insert the sync step. */
		TimeIntegrator_PrependSetupEP( uwCtx->timeIntegrator, 
					       "EulerDeform_IntegrationSetup", 
					       EulerDeform_IntegrationSetup, 
					       "EulerDeform", 
					       edCtx );
	}

	/* Insert the remesh step. Note that this should look for the surface process
	   plugin's time integrator finish routine and ensure we enter the remesh step
	   after that one but before the particle updating routines. */
	TimeIntegrator_PrependFinishEP( uwCtx->timeIntegrator, 
					"EulerDeform_Execute", 
					EulerDeform_Remesh, 
					"EulerDeform", 
					edCtx );
}


void _Underworld_EulerDeform_Destroy( void* component, void* data ) {
	UnderworldContext*	uwCtx = (UnderworldContext*)data;

	assert( component );
	assert( uwCtx );

	/* Clear the lot. */
	/* TODO */
}


Variable* EulerDeform_RegisterLocalNodeCoordsAsVariables( EulerDeform_System* sys, void* _variable_Register, 
							  Variable** variableList )
{
	FeMesh* self              = (FeMesh*)sys->mesh;
	Variable_Register*  variable_Register = (Variable_Register*) _variable_Register;
	Variable*           variable;
	Name                variableName;
	Name                variableNameX;
	Name                variableNameY;
	Name                variableNameZ;

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
		Variable_DataType_Double, 
		Mesh_GetDimSize( self ), 
		&self->topo->domains[MT_VERTEX]->decomp->nLocals, 
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
	TimeIntegrator*		timeIntegrator = (TimeIntegrator*)_timeIntegrator;
	EulerDeform_Context*	edCtx = (EulerDeform_Context*)context;
	unsigned		sys_i;

	FeVariable_SyncShadowValues( edCtx->systems[0].velField );

	/* 
	** We'll need to store side values that we require to be static here, for later
	** return to the mesh.
	*/

	for( sys_i = 0; sys_i < edCtx->nSystems; sys_i++ ) {
		EulerDeform_System*	sys = edCtx->systems + sys_i;

		if( sys->staticSides ) {
			IndexSet	*tmpIndSet;
			RangeSet	*sideSet, *tmpSet;
			unsigned	nInds, *inds;
			unsigned	nDims;
			unsigned	ind_i;

			/* Collect indices of all sides except top surface. */
			sideSet = RangeSet_New();
			tmpSet = RangeSet_New();

			tmpIndSet = RegularMeshUtils_CreateGlobalLeftSet( sys->mesh );
			IndexSet_GetMembers( tmpIndSet, &nInds, &inds );
			RangeSet_SetIndices( sideSet, nInds, inds );
			FreeArray( inds );
			FreeObject( tmpIndSet );

			tmpIndSet = RegularMeshUtils_CreateGlobalRightSet( sys->mesh );
			IndexSet_GetMembers( tmpIndSet, &nInds, &inds );
			RangeSet_SetIndices( tmpSet, nInds, inds );
			RangeSet_Union( sideSet, tmpSet );
			FreeArray( inds );
			FreeObject( tmpIndSet );

			tmpIndSet = RegularMeshUtils_CreateGlobalBottomSet( sys->mesh );
			IndexSet_GetMembers( tmpIndSet, &nInds, &inds );
			RangeSet_SetIndices( tmpSet, nInds, inds );
			RangeSet_Union( sideSet, tmpSet );
			FreeArray( inds );
			FreeObject( tmpIndSet );

			if( sys->mesh->topo->nDims == 3 ) {
				tmpIndSet = RegularMeshUtils_CreateGlobalBottomSet( sys->mesh );
				IndexSet_GetMembers( tmpIndSet, &nInds, &inds );
				RangeSet_SetIndices( tmpSet, nInds, inds );
				RangeSet_Union( sideSet, tmpSet );
				FreeArray( inds );
				FreeObject( tmpIndSet );

				tmpIndSet = RegularMeshUtils_CreateGlobalBottomSet( sys->mesh );
				IndexSet_GetMembers( tmpIndSet, &nInds, &inds );
				RangeSet_SetIndices( tmpSet, nInds, inds );
				RangeSet_Union( sideSet, tmpSet );
				FreeArray( inds );
				FreeObject( tmpIndSet );
			}

			RangeSet_GetIndices( sideSet, &nInds, &inds );
			FreeObject( sideSet );
			FreeObject( tmpSet );

			/* Copy coords to temporary array. */
			nDims = Mesh_GetDimSize( sys->mesh );
			sys->sideCoords = AllocArray2D( double, nInds, nDims );
			for( ind_i = 0; ind_i < nInds; ind_i++ )
				memcpy( sys->sideCoords[ind_i], sys->mesh->verts[inds[ind_i]], nDims * sizeof(double) );
		}
	}

	/* Update advection arrays. */
	for( sys_i = 0; sys_i < edCtx->nSystems; sys_i++ ) {
		EulerDeform_System*	sys = edCtx->systems + sys_i;
		unsigned		nDims;
		unsigned		nLocalNodes;
		unsigned		n_i;

		nDims = Mesh_GetDimSize( sys->mesh );
		nLocalNodes = Mesh_GetLocalSize( sys->mesh, MT_VERTEX );
		for( n_i = 0; n_i < nLocalNodes; n_i++ )
			memcpy( sys->verts + n_i * nDims, sys->mesh->verts[n_i], nDims * sizeof(double) );
	}
}


Bool EulerDeform_TimeDeriv( void* crdAdvector, Index arrayInd, double* timeDeriv ) {
	TimeIntegratee*	self = (TimeIntegratee*)crdAdvector;
	FieldVariable*	velocityField = (FieldVariable*)self->data[0];
	double**	crds = *(double***)self->data[1];
	InterpolationResult  result;

	result = FieldVariable_InterpolateValueAt( velocityField, crds[arrayInd], timeDeriv );

	if ( result == OTHER_PROC || result == OUTSIDE_GLOBAL || isinf(timeDeriv[0]) || isinf(timeDeriv[1]) || 
	     ( velocityField->dim == 3 && isinf(timeDeriv[2]) ) ) 
	{
#if 0
		Journal_Printf( Journal_Register( Error_Type, self->type ),
				"Error in func '%s' for particle with index %u.\n\tPosition (%g, %g, %g)\n\tVelocity here is (%g, %g, %g)."
				"\n\tInterpolation result is %s.\n",
				__func__, array_I, coord[0], coord[1], coord[2], 

				InterpolationResultToStringMap[result]  );
#endif	
		return False;	
	}

	return True;
}


void EulerDeform_Remesh( TimeIntegratee* crdAdvector, EulerDeform_Context* edCtx ) {
	unsigned	sys_i;

	assert( edCtx );

	for( sys_i = 0; sys_i < edCtx->nSystems; sys_i++ ) {
		EulerDeform_System*	sys = edCtx->systems + sys_i;
		double**		oldCrds;
		double**		newCrds;
		unsigned		nDomainNodes;
		unsigned		nDims;
		unsigned		var_i, n_i;

		nDims = Mesh_GetDimSize( sys->mesh );

		/* Update all local coordinates. */
		for( n_i = 0; n_i < Mesh_GetLocalSize( sys->mesh, MT_VERTEX ); n_i++ )
			memcpy( sys->mesh->verts[n_i], sys->verts + n_i * nDims, nDims * sizeof(double) );

		/* Revert side coordinates if required. */
		if( sys->staticSides ) {
			IndexSet	*tmpIndSet;
			RangeSet	*sideSet, *tmpSet;
			unsigned	nInds, *inds;
			unsigned	ind_i;

			/* Collect indices of all sides except top surface. */
			sideSet = RangeSet_New();
			tmpSet = RangeSet_New();

			tmpIndSet = RegularMeshUtils_CreateGlobalLeftSet( sys->mesh );
			IndexSet_GetMembers( tmpIndSet, &nInds, &inds );
			RangeSet_SetIndices( sideSet, nInds, inds );
			FreeArray( inds );
			FreeObject( tmpIndSet );

			tmpIndSet = RegularMeshUtils_CreateGlobalRightSet( sys->mesh );
			IndexSet_GetMembers( tmpIndSet, &nInds, &inds );
			RangeSet_SetIndices( tmpSet, nInds, inds );
			RangeSet_Union( sideSet, tmpSet );
			FreeArray( inds );
			FreeObject( tmpIndSet );

			tmpIndSet = RegularMeshUtils_CreateGlobalBottomSet( sys->mesh );
			IndexSet_GetMembers( tmpIndSet, &nInds, &inds );
			RangeSet_SetIndices( tmpSet, nInds, inds );
			RangeSet_Union( sideSet, tmpSet );
			FreeArray( inds );
			FreeObject( tmpIndSet );

			if( sys->mesh->topo->nDims == 3 ) {
				tmpIndSet = RegularMeshUtils_CreateGlobalBottomSet( sys->mesh );
				IndexSet_GetMembers( tmpIndSet, &nInds, &inds );
				RangeSet_SetIndices( tmpSet, nInds, inds );
				RangeSet_Union( sideSet, tmpSet );
				FreeArray( inds );
				FreeObject( tmpIndSet );

				tmpIndSet = RegularMeshUtils_CreateGlobalBottomSet( sys->mesh );
				IndexSet_GetMembers( tmpIndSet, &nInds, &inds );
				RangeSet_SetIndices( tmpSet, nInds, inds );
				RangeSet_Union( sideSet, tmpSet );
				FreeArray( inds );
				FreeObject( tmpIndSet );
			}

			RangeSet_GetIndices( sideSet, &nInds, &inds );
			FreeObject( sideSet );
			FreeObject( tmpSet );

			/* Copy back coords. */
			for( ind_i = 0; ind_i < nInds; ind_i++ )
				memcpy( sys->mesh->verts[inds[ind_i]], sys->sideCoords[ind_i], nDims * sizeof(double) );
			FreeArray( sys->sideCoords );
		}

		/* Every system should synchronise the mesh coordinates. */
		Mesh_Sync( sys->mesh );
		Mesh_DeformationUpdate( sys->mesh );

		/* Only if remesher specified. */
		if( !sys->remesher ) {
			continue;
		}

		/* Store old coordinates. */
		nDomainNodes = FeMesh_GetNodeDomainSize( sys->remesher->mesh );
		oldCrds = AllocArray2D( double, nDomainNodes, nDims );
		for( n_i = 0; n_i < nDomainNodes; n_i++ )
			memcpy( oldCrds[n_i], sys->remesher->mesh->verts[n_i], nDims * sizeof(double) );

		/* Remesh the system. */
		Execute( sys->remesher, NULL, True );
		Mesh_Sync( sys->mesh );

		/* Shrink wrap the top/bottom surface. */
		if( sys->wrapTop )
			EulerDeform_WrapTopSurface( sys, oldCrds );

		/* Swap old coordinates back in temporarily. */
		newCrds = sys->remesher->mesh->verts;
		sys->remesher->mesh->verts = oldCrds;

		/* Interpolate the variables. */
		for( var_i = 0; var_i < sys->nFields; var_i++ )
			EulerDeform_InterpVar( sys->fields[var_i], sys->vars[var_i], sys->remesher->mesh, newCrds );

		/* Swap back coordinates and free memory. */
		sys->remesher->mesh->verts = newCrds;
		FreeArray( oldCrds );

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
	unsigned	n_i;

	assert( field );
	assert( var );
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
			assert( 0 );
		}
	}

	/* Transfer the new values back to the variable. */
	if( field->fieldComponentCount > 1 ) {
		unsigned	c_i;

		for( n_i = 0; n_i < nLocalNodes; n_i++ ) {
			for( c_i = 0; c_i < field->fieldComponentCount; c_i++ )
				Variable_SetValueAtDouble( var, n_i, c_i, newVals[curValInd++] );
		}
	}
	else {
		for( n_i = 0; n_i < nLocalNodes; n_i++ )
			Variable_SetValueDouble( var, n_i, newVals[curValInd++] );
	}

	/* Free the values array. */
	FreeArray( newVals );
}


void EulerDeform_WrapTopSurface( EulerDeform_System* sys, double** oldCrds ) {
	IJK	ijk;
	Grid*	grm;
	Mesh*	mesh;

	assert( sys );
	assert( oldCrds );

	/* Loop over top internal surface. */
	mesh = sys->remesher->mesh;
	grm = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
					     ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );
	EulerDeform_TopInternalLoop( sys, grm, oldCrds, ijk, 0 );
}

#if 0
void EulerDeform_WrapBottomSurface( EulerDeform_System* sys, double** oldCrds ) {
	IJK	ijk;
	GRM	grm;

	assert( sys );
	assert( oldCrds );

	/* Loop over top internal surface. */
	RegMesh_Generalise( sys->remesher->mesh, &grm );
	EulerDeform_BottomInternalLoop( sys, &grm, oldCrds, ijk, 0 );
}


void EulerDeform_WrapLeftSurface( EulerDeform_System* sys, double** oldCrds ) {
	IJK	ijk;
	GRM	grm;

	assert( sys );
	assert( oldCrds );

	/* Loop over top internal surface. */
	RegMesh_Generalise( sys->remesher->mesh, &grm );
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
	assert( pnt );
	assert( bcs );

	bcs[1] = (pnt - crds[0])/(crds[1] - crds[0]);
	bcs[0] = 1.0 - bcs[1];

	return (bcs[0] >= 0.0 && bcs[0] <= 1.0 && bcs[1] >= 0.0 && bcs[1] <= 1.0);
}


Bool _EulerDeform_LineInterp( const double** crds, const double* pnt, unsigned fromDim, unsigned toDim, 
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


void EulerDeform_TopInternalLoop( EulerDeform_System* sys, Grid* grm, double** oldCrds, unsigned* ijk, unsigned curDim ) {
	unsigned	nDims;
	XYZ		newCrd, oldCrd;
	double*		crds[4];
	double		crds0[3], crds1[3], crds2[3], crds3[3];
	unsigned	centerInd;
	Mesh*		mesh;
	unsigned	nLocalNodes;
	unsigned	ind;

	if( curDim < grm->nDims ) {
		if( curDim == 1 ) {
			ijk[1] = grm->sizes[curDim] - 1;
			EulerDeform_TopInternalLoop( sys, grm, oldCrds, ijk, curDim + 1 );
		}
		else {
			for( ijk[curDim] = 0; ijk[curDim] < grm->sizes[curDim]; ijk[curDim]++ ) {
				EulerDeform_TopInternalLoop( sys, grm, oldCrds, ijk, curDim + 1 );
			}
		}
	}
	else {
		if( grm->nDims == 2 ) {
			mesh = sys->remesher->mesh;
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
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ) );
				memcpy( crds[0], oldCrds[ind], nDims * sizeof(double) );
				memcpy( crds[1], oldCrd, nDims * sizeof(double) );
			}
			else {
				ijk[0]++; ind = Grid_Project( grm, ijk ); ijk[0]--;
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ) );
				memcpy( crds[1], oldCrds[ind], nDims * sizeof(double) );
				memcpy( crds[0], oldCrd, nDims * sizeof(double) );
			}

			/* Interpolate. */
#ifndef NDEBUG
			assert( _EulerDeform_LineInterp( (const double**)crds, newCrd, 0, 1, &mesh->verts[centerInd][1] ) );
#else
			_EulerDeform_LineInterp( (const double**)crds, newCrd, 0, 1, &mesh->verts[centerInd][1] );
#endif
			mesh->verts[centerInd][1] -= 1e-15;
		}
		else if( grm->nDims == 3 ) {
			mesh = sys->remesher->mesh;
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
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ) );
				memcpy( crds[0], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]++; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ) );
				memcpy( crds[1], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]--; ijk[2]++; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ) );
				memcpy( crds[2], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]++; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ) );
				memcpy( crds[3], oldCrds[ind], nDims * sizeof(double) );

				if( _EulerDeform_QuadYInterp( crds, newCrd, &mesh->verts[centerInd][1] ) ) {
					mesh->verts[centerInd][1] -= 1e-13;
					return;
				}
			}

			if( ijk[0] > 0 && ijk[2] < grm->sizes[2] - 1 ) {
				ijk[0]--; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ) );
				memcpy( crds[0], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]++; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ) );
				memcpy( crds[1], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]--; ijk[2]++; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ) );
				memcpy( crds[2], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]++; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ) );
				memcpy( crds[3], oldCrds[ind], nDims * sizeof(double) );

				ijk[2]--;
				if( _EulerDeform_QuadYInterp( crds, newCrd, &mesh->verts[centerInd][1] ) ) {
					mesh->verts[centerInd][1] -= 1e-13;
					return;
				}
			}

			if( ijk[0] < grm->sizes[0] - 1 && ijk[2] > 0 ) {
				ijk[2]--; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ) );
				memcpy( crds[0], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]++; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ) );
				memcpy( crds[1], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]--; ijk[2]++; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ) );
				memcpy( crds[2], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]++; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ) );
				memcpy( crds[3], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]--;
				if( _EulerDeform_QuadYInterp( crds, newCrd, &mesh->verts[centerInd][1] ) ) {
					mesh->verts[centerInd][1] -= 1e-13;
					return;
				}
			}

			if( ijk[0] < grm->sizes[0] - 1 && ijk[2] < grm->sizes[2] - 1 ) {
				ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ) );
				memcpy( crds[0], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]++; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ) );
				memcpy( crds[1], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]--; ijk[2]++; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ) );
				memcpy( crds[2], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]++; ind = Grid_Project( grm, ijk );
				insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ) );
				memcpy( crds[3], oldCrds[ind], nDims * sizeof(double) );

				ijk[0]--; ijk[2]--;
				if( _EulerDeform_QuadYInterp( crds, newCrd, &mesh->verts[centerInd][1] ) ) {
					mesh->verts[centerInd][1] -= 1e-13;
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
void EulerDeform_BottomInternalLoop( EulerDeform_System* sys, GRM* grm, Coord* oldCrds, unsigned* ijk, unsigned curDim ) {
	if( curDim < grm->nDims ) {
		if( curDim == 1 ) {
			ijk[1] = 0;
			EulerDeform_BottomInternalLoop( sys, grm, oldCrds, ijk, curDim + 1 );
		}
		else {
			for( ijk[curDim] = 0; ijk[curDim] < grm->nNodes[curDim]; ijk[curDim]++ ) {
				EulerDeform_BottomInternalLoop( sys, grm, oldCrds, ijk, curDim + 1 );
			}
		}
	}
	else {
		if( grm->nDims == 2 ) {
			XYZ		newCrd, oldCrd;
			unsigned	centerInd;
			Mesh*		mesh = sys->remesher->mesh;

			/* Skip corners. */
			if( (ijk[0] == 0 || ijk[0] == grm->nNodes[0] - 1) && 
			    (ijk[1] == 0 || ijk[1] == grm->nNodes[1] - 1))
			{
				return;
			}

			/* Get old and new coordinate. */
			GRM_Project( grm, ijk, &centerInd );
			newCrd[0] = mesh->nodeCoord[centerInd][0];
			newCrd[1] = mesh->nodeCoord[centerInd][1];
			oldCrd[0] = oldCrds[centerInd][0];
			oldCrd[1] = oldCrds[centerInd][1];

			/* Are we left or right? */
			if( newCrd[0] < oldCrd[0] ) {
				XYZ		leftCrd;
				unsigned	leftInd;
				double		a0, a1;

				/* Get left old coord. */
				ijk[0]--;
				GRM_Project( grm, ijk, &leftInd );
				ijk[0]++;
				leftCrd[0] = oldCrds[leftInd][0];
				leftCrd[1] = oldCrds[leftInd][1];

				/* Calc barycenter. */
				a1 = (newCrd[0] - leftCrd[0]) / (oldCrd[0] - leftCrd[0]);
				a0 = 1.0 - a1;
				mesh->nodeCoord[centerInd][1] = a0 * leftCrd[1] + a1 * oldCrd[1] + 1e-15;
			}
			else {
				XYZ		rightCrd;
				unsigned	rightInd;
				double		a0, a1;

				/* Get right old coord. */
				ijk[0]++;
				GRM_Project( grm, ijk, &rightInd );
				ijk[0]--;
				rightCrd[0] = oldCrds[rightInd][0];
				rightCrd[1] = oldCrds[rightInd][1];

				/* Calc barycenter. */
				a1 = (newCrd[0] - oldCrd[0]) / (rightCrd[0] - oldCrd[0]);
				a0 = 1.0 - a1;
				mesh->nodeCoord[centerInd][1] = a0 * oldCrd[1] + a1 * rightCrd[1] + 1e-15;
			}
		}
		else if( grm->nDims == 3 ) {
			XYZ		newCrd, oldCrd;
			Coord		crds[4];
			unsigned	centerInd;
			Mesh*		mesh = sys->remesher->mesh;
			unsigned	ind;

			/* Skip corners. */
			if( (ijk[0] == 0 || ijk[0] == grm->nNodes[0] - 1) && 
			    (ijk[1] == 0 || ijk[1] == grm->nNodes[1] - 1) && 
			    (ijk[2] == 0 || ijk[2] == grm->nNodes[2] - 1))
			{
				return;
			}

			/* Get old and new coordinate. */
			GRM_Project( grm, ijk, &centerInd );
			newCrd[0] = mesh->nodeCoord[centerInd][0];
			newCrd[1] = mesh->nodeCoord[centerInd][1];
			newCrd[2] = mesh->nodeCoord[centerInd][2];
			oldCrd[0] = oldCrds[centerInd][0];
			oldCrd[1] = oldCrds[centerInd][1];
			oldCrd[2] = oldCrds[centerInd][2];

			/* Figure out what qudrant we're in. */
			if( ijk[0] > 0 && ijk[2] > 0 ) {
				ijk[0]--; ijk[2]--;	GRM_Project( grm, ijk, &ind ); memcpy( crds[0], oldCrds[ind], sizeof(Coord) );
				ijk[0]++;		GRM_Project( grm, ijk, &ind ); memcpy( crds[1], oldCrds[ind], sizeof(Coord) );
				ijk[0]--; ijk[2]++;	GRM_Project( grm, ijk, &ind ); memcpy( crds[2], oldCrds[ind], sizeof(Coord) );
				ijk[0]++;		GRM_Project( grm, ijk, &ind ); memcpy( crds[3], oldCrds[ind], sizeof(Coord) );
				if( _EulerDeform_QuadYInterp( crds, newCrd, &mesh->nodeCoord[centerInd][1] ) ) {
					mesh->nodeCoord[centerInd][1] += 1e-15;
					return;
				}
			}

			if( ijk[0] > 0 && ijk[2] < grm->nNodes[2] - 1 ) {
				ijk[0]--;		GRM_Project( grm, ijk, &ind ); memcpy( crds[0], oldCrds[ind], sizeof(Coord) );
				ijk[0]++;		GRM_Project( grm, ijk, &ind ); memcpy( crds[1], oldCrds[ind], sizeof(Coord) );
				ijk[0]--; ijk[2]++;	GRM_Project( grm, ijk, &ind ); memcpy( crds[2], oldCrds[ind], sizeof(Coord) );
				ijk[0]++;		GRM_Project( grm, ijk, &ind ); memcpy( crds[3], oldCrds[ind], sizeof(Coord) );
				ijk[2]--;
				if( _EulerDeform_QuadYInterp( crds, newCrd, &mesh->nodeCoord[centerInd][1] ) ) {
					mesh->nodeCoord[centerInd][1] += 1e-15;
					return;
				}
			}

			if( ijk[0] < grm->nNodes[0] - 1 && ijk[2] > 0 ) {
				ijk[2]--;		GRM_Project( grm, ijk, &ind ); memcpy( crds[0], oldCrds[ind], sizeof(Coord) );
				ijk[0]++;		GRM_Project( grm, ijk, &ind ); memcpy( crds[1], oldCrds[ind], sizeof(Coord) );
				ijk[0]--; ijk[2]++;	GRM_Project( grm, ijk, &ind ); memcpy( crds[2], oldCrds[ind], sizeof(Coord) );
				ijk[0]++;		GRM_Project( grm, ijk, &ind ); memcpy( crds[3], oldCrds[ind], sizeof(Coord) );
				ijk[0]--;
				if( _EulerDeform_QuadYInterp( crds, newCrd, &mesh->nodeCoord[centerInd][1] ) ) {
					mesh->nodeCoord[centerInd][1] += 1e-15;
					return;
				}
			}

			if( ijk[0] < grm->nNodes[0] - 1 && ijk[2] < grm->nNodes[2] - 1 ) {
				GRM_Project( grm, ijk, &ind ); memcpy( crds[0], oldCrds[ind], sizeof(Coord) );
				ijk[0]++;		GRM_Project( grm, ijk, &ind ); memcpy( crds[1], oldCrds[ind], sizeof(Coord) );
				ijk[0]--; ijk[2]++;	GRM_Project( grm, ijk, &ind ); memcpy( crds[2], oldCrds[ind], sizeof(Coord) );
				ijk[0]++;		GRM_Project( grm, ijk, &ind ); memcpy( crds[3], oldCrds[ind], sizeof(Coord) );
				ijk[0]--; ijk[2]--;
				if( _EulerDeform_QuadYInterp( crds, newCrd, &mesh->nodeCoord[centerInd][1] ) ) {
					mesh->nodeCoord[centerInd][1] += 1e-15;
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
			Mesh*		mesh = sys->remesher->mesh;

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
				mesh->nodeCoord[centerInd][0] = a0 * leftCrd[0] + a1 * oldCrd[0] + 1e-15;
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
				mesh->nodeCoord[centerInd][0] = a0 * oldCrd[0] + a1 * rightCrd[0] + 1e-15;
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
