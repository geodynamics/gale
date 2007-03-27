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
** $Id:  $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>

#include "types.h"
#include "Polynomial.h"
#include "Quadrature.h"
#include "LevelSet.h"
#include "Context.h"
#include "LevelSetPlugin.h"


const Type		StgFEM_LevelSetPlugin_Type = "LevelSetPlugin";
ExtensionInfo_Index	LevelSetPlugin_ContextHandle;

const char*	LEVELSET_PLUGIN_TAG = "LevelSet";


Index _StgFEM_LevelSetPlugin_Register( PluginsManager* pluginsMgr ) {
	/*
	** Register the level set component.
	*/

	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), 
				   LevelSet_Type, 
				   "0", 
				   (Stg_Component_DefaultConstructorFunction*)LevelSet_DefaultNew );
	RegisterParent( LevelSet_Type, FieldVariable_Type );


	/* Complete plugin bit. */
	return PluginsManager_Submit( pluginsMgr, 
				      StgFEM_LevelSetPlugin_Type, 
				      "0", 
				      _StgFEM_LevelSetPlugin_DefaultNew );
}


void* _StgFEM_LevelSetPlugin_DefaultNew( Name name ) {
	return _Codelet_New( sizeof(Codelet), 
			     StgFEM_LevelSetPlugin_Type, 
			     _Codelet_Delete, 
			     _Codelet_Print, 
			     _Codelet_Copy, 
			     _StgFEM_LevelSetPlugin_DefaultNew, 
			     _StgFEM_LevelSetPlugin_Construct, 
			     _StgFEM_LevelSetPlugin_Build, 
			     _Codelet_Initialise, 
			     _Codelet_Execute, 
			     _StgFEM_LevelSetPlugin_Destroy, 
			     name );
}


void _StgFEM_LevelSetPlugin_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	FiniteElementContext*	feCtx;
	LevelSetPlugin_Context*	lsCtx;

	assert( component );
	assert( cf );

	/* Retrieve context. */
	feCtx = (FiniteElementContext*)Stg_ComponentFactory_ConstructByName( 
		cf,
		"context", 
		FiniteElementContext, 
		True,
		data );

	Journal_DPrintf( StgFEM_Debug, "In: %s( void* )\n", __func__ );


	/*
	** Extend the context for storing level sets.
	*/

	LevelSetPlugin_ContextHandle = ExtensionManager_Add( feCtx->extensionMgr,
							     StgFEM_LevelSetPlugin_Type,
							     sizeof(LevelSetPlugin_Context) );
	lsCtx = ExtensionManager_Get( feCtx->extensionMgr,
				      feCtx, 
				      LevelSetPlugin_ContextHandle );
	lsCtx->nLSets = 0;
	lsCtx->lSets = NULL;


	/*
	** Append to the solve phase, note that this must execute after all SLEs.
	*/

	EntryPoint_Append( Context_GetEntryPoint( feCtx, AbstractContext_EP_Solve ),
			   "default", 
			   _StgFEM_LevelSetPlugin_Solve, 
			   StgFEM_LevelSetPlugin_Type );
}


void _StgFEM_LevelSetPlugin_Build( void* component, void* data ) {
	FiniteElementContext*	feCtx;

	assert( component );
	assert( data );

	/* Get the contexts. */
	feCtx = (FiniteElementContext*)data;
}


void _StgFEM_LevelSetPlugin_Destroy( void* component, void* data ) {
	FiniteElementContext*	feCtx;
	LevelSetPlugin_Context*	lsCtx;

	assert( component );
	assert( data );

	/* Get the contexts. */
	feCtx = (FiniteElementContext*)data;
	lsCtx = ExtensionManager_Get( feCtx->extensionMgr,
				      feCtx, 
				      LevelSetPlugin_ContextHandle );
	FreeArray( lsCtx->lSets );
}


void _StgFEM_LevelSetPlugin_Solve( void* ctx ) {
	FiniteElementContext*	feCtx = (FiniteElementContext*)ctx;
	LevelSetPlugin_Context*	lsCtx;
	unsigned		ls_i;

	assert( feCtx );

	lsCtx = ExtensionManager_Get( feCtx->extensionMgr,
				      feCtx, 
				      LevelSetPlugin_ContextHandle );

	for( ls_i = 0; ls_i < lsCtx->nLSets; ls_i++ ) {
		Execute( lsCtx->lSets[ls_i], feCtx, True );
	}
}
