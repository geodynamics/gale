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
** $Id: FeVariable.c 796 2007-03-27 07:22:24Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include "units.h"
#include "types.h"
#include "shortcuts.h"
#include "ElementType.h"
#include "ElementType_Register.h"
#include "Element.h"
#include "Mesh.h"
#include "FeEquationNumber.h"
#include "FeVariable.h"
#include "LinkedDofInfo.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

const Type   FeVariable_Type = "FeVariable";
const Name   defaultFeVariableFeEquationNumberName = "defaultFeVariableFeEqName";
const char*  StgFEM_Native_ImportExportType = "StgFEM_Native";

/* MPI Tags */
static const int DOF_VALUES_TAG = 10;

/* Global objects */
Stg_ObjectList*    FeVariable_FileFormatImportExportList = NULL;

void* FeVariable_DefaultNew( Name name )
{
	return _FeVariable_New(
		sizeof(FeVariable),
		FeVariable_Type,
		_FeVariable_Delete,
		_FeVariable_Print,
		_FeVariable_Copy,
		(Stg_Component_DefaultConstructorFunction*)FeVariable_DefaultNew,
		_FeVariable_Construct,
		_FeVariable_Build, 
		_FeVariable_Initialise,
		_FeVariable_Execute,
		_FeVariable_Destroy,
		name,
		False,
		_FeVariable_InterpolateValueAt,
		_FeVariable_GetMinGlobalFieldMagnitude,
		_FeVariable_GetMaxGlobalFieldMagnitude,
		_FeVariable_GetMinAndMaxLocalCoords,
		_FeVariable_GetMinAndMaxGlobalCoords,
		_FeVariable_InterpolateNodeValuesToElLocalCoord,
		_FeVariable_GetValueAtNode,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		0,
		0,
		True,
		NULL,
		NULL,
		MPI_COMM_WORLD,
		NULL );
}


FeVariable* FeVariable_New(
		Name                                            name,
		void*                                           feMesh,
		void*                                           geometryMesh,
		DofLayout*                                      dofLayout, 
		void*                                          	bcs,
		void*                                          	ics,
		void*                                           linkedDofInfo,
		Dimension_Index                                 dim,
		Bool                                            isCheckpointedAndReloaded,
		const char* const                               importFormatType,
		const char* const                               exportFormatType,
		FieldVariable_Register*                         fV_Register )		
{
	return FeVariable_New_Full(
		name,
		feMesh,
		geometryMesh,
		dofLayout,
		bcs,
		ics,
		linkedDofInfo,
		NULL,
		dofLayout->_totalVarCount,
		dim,
	    isCheckpointedAndReloaded,
		importFormatType,
		exportFormatType,
		((Mesh*)feMesh)->layout->decomp->communicator,
		fV_Register );
}


FeVariable* FeVariable_New_FromTemplate(
		Name                                            name,
		void*                                          	_templateFeVariable,
		DofLayout*                                      dofLayout, 
		void*                                          	ics,
		const char* const                               importFormatType,
		const char* const                               exportFormatType,
		FieldVariable_Register*                         fV_Register )
{
	FeVariable*	templateFeVariable = _templateFeVariable;
	FeVariable*	newFeVariable = NULL;
	
	newFeVariable = FeVariable_New_Full(
		name, 
		templateFeVariable->feMesh,
		templateFeVariable->geometryMesh,
		dofLayout,
		templateFeVariable->bcs,
		ics,
		templateFeVariable->linkedDofInfo,
		templateFeVariable,
		templateFeVariable->fieldComponentCount,
		templateFeVariable->dim,
		templateFeVariable->isCheckpointedAndReloaded,
		importFormatType,
		exportFormatType,
		templateFeVariable->communicator,
		fV_Register );

	newFeVariable->templateFeVariable = templateFeVariable;
	return newFeVariable;
}


FeVariable* FeVariable_New_Full(
		Name                                            name,
		void*                                           feMesh,
		void*                                           geometryMesh,
		DofLayout*                                      dofLayout, 
		void*                                          	bcs,
		void*                                          	ics,
		void*                                           linkedDofInfo,
		void*                                          	templateFeVariable,
		Index                                           fieldComponentCount,
		Dimension_Index                                 dim,
		Bool                                            isCheckpointedAndReloaded,
		const char* const                               importFormatType,
		const char* const                               exportFormatType,
		MPI_Comm                                        communicator,
		FieldVariable_Register*                         fV_Register )		
{
	return _FeVariable_New(
		sizeof(FeVariable),
		FeVariable_Type,
		_FeVariable_Delete,
		_FeVariable_Print,
		_FeVariable_Copy,
		(Stg_Component_DefaultConstructorFunction*)FeVariable_DefaultNew,
		_FeVariable_Construct,
		_FeVariable_Build, 
		_FeVariable_Initialise,
		_FeVariable_Execute,
		_FeVariable_Destroy,
		name,
		True,
		_FeVariable_InterpolateValueAt,
		_FeVariable_GetMinGlobalFieldMagnitude,
		_FeVariable_GetMaxGlobalFieldMagnitude,
		_FeVariable_GetMinAndMaxLocalCoords,
		_FeVariable_GetMinAndMaxGlobalCoords,
		_FeVariable_InterpolateNodeValuesToElLocalCoord,
		_FeVariable_GetValueAtNode,
		feMesh,
		geometryMesh,
		dofLayout,
		bcs,
		ics,
		linkedDofInfo,
		templateFeVariable,
		fieldComponentCount,
		dim,
		isCheckpointedAndReloaded,
		importFormatType,
		exportFormatType,
		communicator,
		fV_Register );
}


FeVariable* _FeVariable_New( 
		SizeT                                           _sizeOfSelf,
		Type                                            type,
		Stg_Class_DeleteFunction*                       _delete,
		Stg_Class_PrintFunction*                        _print,
		Stg_Class_CopyFunction*                         _copy, 
		Stg_Component_DefaultConstructorFunction*       _defaultConstructor,
		Stg_Component_ConstructFunction*                _construct,
		Stg_Component_BuildFunction*                    _build,
		Stg_Component_InitialiseFunction*               _initialise,
		Stg_Component_ExecuteFunction*                  _execute,
		Stg_Component_DestroyFunction*                  _destroy,
		Name                                            name,
		Bool                                            initFlag,
		FieldVariable_InterpolateValueAtFunction*       _interpolateValueAt,
		FieldVariable_GetValueFunction*                 _getMinGlobalFieldMagnitude,
		FieldVariable_GetValueFunction*                 _getMaxGlobalFieldMagnitude,
		FieldVariable_GetCoordFunction*                 _getMinAndMaxLocalCoords,
		FieldVariable_GetCoordFunction*                 _getMinAndMaxGlobalCoords,		
		FeVariable_InterpolateWithinElementFunction*    _interpolateWithinElement,	
		FeVariable_GetValueAtNodeFunction*              _getValueAtNode,	
		void*                                           feMesh,
		void*                                           geometryMesh,
		DofLayout*                                      dofLayout, 
		void*                                          	bcs,
		void*                                          	ics,
		void*                                           linkedDofInfo,
		void*                                          	templateFeVariable,
		Index                                           fieldComponentCount,
		Dimension_Index                                 dim,
		Bool                                            isCheckpointedAndReloaded,
		const char* const                               importFormatType,
		const char* const                               exportFormatType,
		MPI_Comm                                        communicator,
		FieldVariable_Register*                         fV_Register )
{
	FeVariable*		   self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(FeVariable) );
	
	self = (FeVariable*)
		_FieldVariable_New(
				_sizeOfSelf,
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
				name,
				initFlag,
				_interpolateValueAt,
				_getMinGlobalFieldMagnitude,
				_getMaxGlobalFieldMagnitude,
				_getMinAndMaxLocalCoords,
				_getMinAndMaxGlobalCoords,
				fieldComponentCount,
				dim,
				isCheckpointedAndReloaded,
				communicator,
				fV_Register );
	
	/* General info */
	
	/* Virtual functions */
	self->_interpolateWithinElement = _interpolateWithinElement;
	self->_getValueAtNode           = _getValueAtNode;
	
	/* FeVariable info */
	if( initFlag ){
		_FeVariable_Init( self, feMesh, geometryMesh, dofLayout, bcs, ics, linkedDofInfo,
			templateFeVariable, importFormatType, exportFormatType );
	}
	
	return self;
}


void _FeVariable_Init( 
		FeVariable*                                     self,
		void*                                           feMesh,
		void*                                           geometryMesh,
		DofLayout*                                      dofLayout, 
		void*                                           bcs,
		void*                                           ics,
		void*                                           linkedDofInfo,
		void*                                           templateFeVariable,
		const char* const                               importFormatType,
		const char* const                               exportFormatType )
{
	/* General and Virtual info should already be set */
	
	/* FeVariable info */
	self->isConstructed = True;
	self->debug = Stream_RegisterChild( StgFEM_Discretisation_Debug, self->type );
	self->feMesh = Stg_CheckType( feMesh, FiniteElement_Mesh );
	/* Set pointer for geometry mesh - if none is provided then it'll use the feMesh */
	self->geometryMesh = ( geometryMesh ? 
			Stg_CheckType( geometryMesh, FiniteElement_Mesh ) : 
			Stg_CheckType( feMesh, FiniteElement_Mesh ) );
	self->dofLayout = dofLayout;
	if ( bcs )
		self->bcs = Stg_CheckType( bcs, VariableCondition );
	if ( ics )
		self->ics = Stg_CheckType( ics, VariableCondition );
	if ( linkedDofInfo )
		self->linkedDofInfo = Stg_CheckType( linkedDofInfo, LinkedDofInfo );
	self->shadowValuesSynchronised = False;

	if ( templateFeVariable )
		self->templateFeVariable = Stg_CheckType( templateFeVariable, FeVariable );
	if ( self->templateFeVariable ) {
		self->eqNum = self->templateFeVariable->eqNum;
	}
	else {
		self->eqNum = FeEquationNumber_New( defaultFeVariableFeEquationNumberName, self->feMesh,
			self->dofLayout, self->bcs, linkedDofInfo );
	}

	self->importFormatType = StG_Strdup( importFormatType );
	/* check the given file format is actually among the registered list. If not, print them and exit. */
	if( NULL == Stg_ObjectList_Get( FeVariable_FileFormatImportExportList, (char*)importFormatType ) ) {
		Stream*    errorStream = Journal_Register( Error_Type, self->type );
		
		Journal_Printf( errorStream, "Error - in %s() - for FeVariable \"%s\": you specified this "
			"FeVariable's import type as %s, which is not in the register of known "
			"FeVariable import/export types\n.", __func__, self->name, importFormatType );
		Journal_Printf( errorStream, "Currently registered import/export types are:\n" );
		Stream_Indent( errorStream );
		Stg_ObjectList_PrintAllEntryNames( FeVariable_FileFormatImportExportList, errorStream );
		Stream_UnIndent( errorStream );
		
		Journal_Firewall( 0, errorStream, "Exiting.\n" );
	}	

	self->exportFormatType = StG_Strdup( exportFormatType );
	/* check the given file format is actually among the registered list. If not, print them and exit. */
	if( NULL == Stg_ObjectList_Get( FeVariable_FileFormatImportExportList, (char*)exportFormatType ) ) {
		Stream*    errorStream = Journal_Register( Error_Type, self->type );
		
		Journal_Printf( errorStream, "Error - in %s() - for FeVariable \"%s\": you specified this "
			"FeVariable's export type as %s, which is not in the register of known "
			"FeVariable import/export types\n.", __func__, self->name, exportFormatType );
		Journal_Printf( errorStream, "Currently registered import/export types are:\n" );
		Stream_Indent( errorStream );
		Stg_ObjectList_PrintAllEntryNames( FeVariable_FileFormatImportExportList, errorStream );
		Stream_UnIndent( errorStream );
		
		Journal_Firewall( 0, errorStream, "Exiting.\n" );
	}	
}


void _FeVariable_Delete( void* variable ) {
	FeVariable* self = (FeVariable*)variable;
	
	Journal_DPrintf( self->debug, "In %s- for \"%s\":\n", __func__, self->name );
	Stream_IndentBranch( StgFEM_Debug );
	if( self->eqNum && ( NULL == self->templateFeVariable ) ) {
		Stg_Class_Delete( self->eqNum );
	}
	Memory_Free( self->importFormatType );
	Memory_Free( self->exportFormatType );
	/* feMesh bc and doflayout are purposely not deleted */

	/* Stg_Class_Delete parent*/
	_Stg_Component_Delete( self );
	Stream_UnIndentBranch( StgFEM_Debug );
}


/* --- Virtual Function Implementations --- */

void _FeVariable_Print( void* variable, Stream* stream ) {
	FeVariable* self = (FeVariable*)variable;
	
	/* General info */
	Journal_Printf( stream, "FeVariable (ptr): %p\n", self );
	
	/* Print parent */
	_Stg_Component_Print( self, stream );
	
	/* Virtual info */
	
	/* FeVariable info */
	Print( self->feMesh, stream );
	if ( self->dofLayout ) {
		Print( self->dofLayout, stream );
	}
	else {
		Journal_Printf( stream, "\tdofLayout: (null)... not provided (may be Operator type)\n" );
	}
	if ( self->bcs ) {
		Print( self->bcs, stream );
	}
	else {
		Journal_Printf( stream, "\tbcs: (null)... not provided (may be Operator type)\n" );
	}
	if ( self->ics ) {
		Print( self->ics, stream );
	}
	else {
		Journal_Printf( stream, "\tics: (null)... not provided (may be Operator type)\n" );
	}

	if ( self->linkedDofInfo ) {
		Print( self->linkedDofInfo, stream );
	}
	else {
		Journal_Printf( stream, "\tlinkedDofInfo: (null)... not provided\n" );
	}
	
	if( self->eqNum ) {
		Print( self->eqNum, stream );
	}
	else {
		Journal_Printf( stream, "\teqNum: (null)... not built yet\n" );
	}
}


void* _FeVariable_Copy( void* feVariable, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	FeVariable*	self = (FeVariable*)feVariable;
	FeVariable*	newFeVariable;
	PtrMap*		map = ptrMap;
	Bool		ownMap = False;
	
	if( !map ) {
		map = PtrMap_New( 10 );
		ownMap = True;
	}
	
	newFeVariable = _FieldVariable_Copy( self, dest, deep, nameExt, map );

	newFeVariable->templateFeVariable = self->templateFeVariable;
	
	if( deep ) {
		newFeVariable->debug = self->debug; 
		newFeVariable->feMesh = (FiniteElement_Mesh*)Stg_Class_Copy( self->feMesh, NULL, deep, nameExt, map );
		newFeVariable->dofLayout = (DofLayout*)Stg_Class_Copy( self->dofLayout, NULL, deep, nameExt, map );
		newFeVariable->bcs = (VariableCondition*)Stg_Class_Copy( self->bcs, NULL, deep, nameExt, map );
		newFeVariable->ics = self->ics ? (VariableCondition*)Stg_Class_Copy( self->ics, NULL, deep, nameExt, map ) : NULL;
		if ( self->linkedDofInfo == NULL ) {
			newFeVariable->linkedDofInfo = NULL;
		}
		else {
			newFeVariable->linkedDofInfo = (LinkedDofInfo*)Stg_Class_Copy( self->linkedDofInfo, NULL,
				deep, nameExt, map );
		}

		if ( self->templateFeVariable ) {
			newFeVariable->eqNum = self->eqNum;
		}
		else {
			newFeVariable->eqNum = (FeEquationNumber*)Stg_Class_Copy( self->eqNum, NULL, deep, nameExt, map );
		}
	}
	else {
		newFeVariable->debug = self->debug;
		newFeVariable->feMesh = self->feMesh;
		newFeVariable->geometryMesh = self->geometryMesh;
		newFeVariable->dofLayout = self->dofLayout;
		newFeVariable->bcs = self->bcs;
		newFeVariable->ics = self->ics;
		newFeVariable->linkedDofInfo = self->linkedDofInfo;
		newFeVariable->eqNum = self->eqNum;
	}
	
	if( ownMap ) {
		Stg_Class_Delete( map );
	}
	
	return (void*)newFeVariable;
}


void _FeVariable_Build( void* variable, void* data ) {
	FeVariable* self = (FeVariable*)variable;
	DiscretisationContext*   context = (DiscretisationContext*)data;
	
	if ( False == self->isBuilt ) {
		self->isBuilt = True;
	
		Journal_DPrintf( self->debug, "In %s- for %s:\n", __func__, self->name );
		Stream_IndentBranch( StgFEM_Debug );

		/* build the BCs */
		Build( self->feMesh, data, False );
		Build( self->dofLayout, data, False );
		if ( self->bcs ){
			Build( self->bcs, data, False );
		}
		/* only bother building the ics specified via XML/construct if we are not in restart mode
		  - otherwise, we will use the checkpointed values anyway */
		if ( self->ics && !(context && (True == context->loadFromCheckPoint) ) ) {
			Build( self->ics, data, False );
		}
		if ( self->linkedDofInfo ) {
			Build( self->linkedDofInfo, data, False );
		}
		
		/* build the e.q. number array */
		FeEquationNumber_Build( self->eqNum );
		
		Stream_UnIndentBranch( StgFEM_Debug );
	}
}

void _FeVariable_Construct( void* variable, Stg_ComponentFactory* cf, void* data ) 
{
	FeVariable*         self          = (FeVariable*)variable;
	FiniteElement_Mesh* feMesh        = NULL;
	FiniteElement_Mesh* geometryMesh  = NULL;
	DofLayout*          dofLayout     = NULL;
	VariableCondition*  bc            = NULL;
	VariableCondition*  ic            = NULL;
	LinkedDofInfo*      linkedDofInfo = NULL;
	char*               importFormatType = NULL;
	char*               exportFormatType = NULL;

	_FieldVariable_Construct( self, cf, data );

	feMesh        = Stg_ComponentFactory_ConstructByKey( cf, self->name, "FEMesh",        FiniteElement_Mesh, True, data );
	geometryMesh  = Stg_ComponentFactory_ConstructByKey( cf, self->name, "GeometryMesh",  FiniteElement_Mesh, False, data );
	dofLayout     = Stg_ComponentFactory_ConstructByKey( cf, self->name, DofLayout_Type,  DofLayout,          True, data );
	importFormatType =  Stg_ComponentFactory_GetString( cf, self->name, "importFormatType",
		StgFEM_Native_ImportExportType );
	exportFormatType =  Stg_ComponentFactory_GetString( cf, self->name, "exportFormatType",
		StgFEM_Native_ImportExportType );

	ic            = Stg_ComponentFactory_ConstructByKey( cf, self->name, "IC",            VariableCondition,  False, data );
	bc            = Stg_ComponentFactory_ConstructByKey( cf, self->name, "BC",            VariableCondition,  False, data );
	linkedDofInfo = Stg_ComponentFactory_ConstructByKey( cf, self->name, "LinkedDofInfo", LinkedDofInfo,      False, data );

	self->fieldComponentCount = dofLayout->_totalVarCount;

	_FeVariable_Init( self, feMesh, geometryMesh, dofLayout, bc, ic, linkedDofInfo, NULL,
		importFormatType, exportFormatType );
}

void _FeVariable_Initialise( void* variable, void* data ) {
	FeVariable*              self = (FeVariable*)variable;
	DiscretisationContext*   context = (DiscretisationContext*)data;
	
	Journal_DPrintf( self->debug, "In %s- for %s:\n", __func__, self->name );
	Stream_IndentBranch( StgFEM_Debug );
	
	/* do basic mesh initialisation */
	Initialise( self->feMesh, data, False );
	Initialise( self->dofLayout, data, False );

	if ( self->linkedDofInfo ) {
		Initialise( self->linkedDofInfo, data, False );
	}
	
	FeEquationNumber_Initialise( self->eqNum );
	/*Seting up whether to load from checkpointing */
	if ( self->ics || ((context && (True == context->loadFromCheckPoint) )&& (self->isCheckpointedAndReloaded)) )  {
		Journal_DPrintf( self->debug, "applying the I.C.s for this Variable:\n" ); 
		Stream_Indent( self->debug );
	
		if ( self->ics && !(context && (True == context->loadFromCheckPoint) ) ) {
			Journal_DPrintf( self->debug, "regular (non-restart) mode -> applying ICs specified in XML/constructor\n" );
			Initialise( self->ics, data, False );
			VariableCondition_Apply( self->ics, data );
		}
		else {
			char*                  inputPathString = NULL;
			Index                  inputStrLen = 0;

			Journal_DPrintf( self->debug, "restart from checkpoint mode -> loading checkpointed "
				"nodal values as initial conditions, ignoring ics specified via XML/constructor\n" );

			inputStrLen = strlen(context->outputPath) + 1 + 1;
			if ( strlen(context->checkPointPrefixString) > 0 ) {
				inputStrLen += strlen(context->checkPointPrefixString) + 1;
			}
			inputPathString = Memory_Alloc_Array( char, inputStrLen, "inputPathString" );

			if ( strlen(context->checkPointPrefixString) > 0 ) {
				sprintf( inputPathString, "%s/%s.", context->outputPath, context->checkPointPrefixString );
			}
			else {
				sprintf( inputPathString, "%s/", context->outputPath );
			}

			FeVariable_ReadFromFile( self, inputPathString, context->restartTimestep );
		
			/* TODO: maybe we want a mechanism in future to over-ride the checkpointed ICs in certain regions too
			so the user can introduce new phenomena into the model */
		}	
	}
	Stream_UnIndent( self->debug );

	if ( self->bcs ) {
		Initialise( self->bcs, data, False );
		Journal_DPrintf( self->debug, "applying the B.C.s for this Variable.\n" ); 
		VariableCondition_Apply( self->bcs, data );
	}

	Stream_UnIndentBranch( StgFEM_Debug );
}


void FeVariable_ApplyBCs( void* variable, void* data ) {
	FeVariable* self = (FeVariable*)variable;

	if ( self->bcs ) {
		Journal_DPrintf( self->debug, "In %s- for %s:\n", __func__, self->name );
		Journal_DPrintf( self->debug, "applying the B.C.s for this Variable.\n" ); 
		VariableCondition_Apply( self->bcs, data );
	}
}


void _FeVariable_Execute( void* variable, void* data ) {
}

void _FeVariable_Destroy( void* variable, void* data ) {
}

void FeVariable_PrintLocalDiscreteValues( void* variable, Stream* stream ) {
	FeVariable* self = (FeVariable*)variable;

	Journal_Printf( stream, "In %s: for FeVariable \"%s\":\n", __func__, self->name );

	_FeVariable_PrintLocalOrDomainValues( variable, self->feMesh->nodeLocalCount, stream );
}


unsigned _FeVariable_ClosestNode( FeVariable* self, Coord crd ) {
	Bool		done;
	Mesh*		mesh = (Mesh*)self->feMesh;
	Coord*		nodeCrds = mesh->nodeCoord;
	unsigned	curNode;
	unsigned	nDims = self->dim;

	/* Begin somewhere in the middle. */
	curNode = mesh->nodeLocalCount / 2;

	/* Loop until we've found closest local node. */
	do {
		unsigned	nNbrs = mesh->nodeNeighbourCountTbl[curNode];
		unsigned*	nbrs = mesh->nodeNeighbourTbl[curNode];
		double		dist;
		double		tmp;
		unsigned	nbr_i, d_i;

		/* Assume we'll be done after this loop. */
		done = True;

		/* Calc distance squared to current node. */
		tmp = nodeCrds[curNode][0] - crd[0];
		dist = tmp * tmp;
		for( d_i = 1; d_i < nDims; d_i++ ) {
			tmp = nodeCrds[curNode][d_i] - crd[d_i];
			dist += tmp * tmp;
		}

		/* Compare to neighbours. */
		for( nbr_i = 0; nbr_i < nNbrs; nbr_i++ ) {
			double	nbrDist;

			/* Just in case... */
			if( nbrs[nbr_i] >= mesh->nodeLocalCount )
				continue;

			tmp = nodeCrds[nbrs[nbr_i]][0] - crd[0];
			nbrDist = tmp * tmp;
			for( d_i = 1; d_i < nDims; d_i++ ) {
				tmp = nodeCrds[nbrs[nbr_i]][d_i] - crd[d_i];
				nbrDist += tmp * tmp;
			}

			if( nbrDist < dist ) {
				curNode = nbrs[nbr_i];
				dist = nbrDist;
				done = False;
			}
		}
	}
	while( !done );

	return curNode;
}


InterpolationResult _FeVariable_InterpolateValueAt( void* variable, Coord globalCoord, double* value ) {
	FeVariable*		self = (FeVariable*)variable;
	Element_DomainIndex	elementCoordIn = (unsigned)-1;
	Coord			elLocalCoord={0,0,0};
	InterpolationResult	retValue;


	retValue = FeVariable_GetElementLocalCoordAtGlobalCoord( self, globalCoord, elLocalCoord, &elementCoordIn );
	
	if ( retValue == LOCAL ) {
		/* Now interpolate the value at that coordinate, using shape functions */
		self->_interpolateWithinElement( self, elementCoordIn, elLocalCoord, value );
	}
	else if ( retValue == SHADOW ) {
		if ( False == self->shadowValuesSynchronised ) {
			Stream* warningStr = Journal_Register( Error_Type, self->type );
			Journal_Printf( warningStr, "Warning - in %s: user asking to interpolate a value at "
				"coord (%g,%g,%g), which is in shadow space, but "
				"FeVariable_SyncShadowValues() hasn't been called yet.\n", 
				__func__, globalCoord[0], globalCoord[1], globalCoord[2] );
			return retValue;		
		}
		/* Now interpolate the value at that coordinate, using shape functions */
		self->_interpolateWithinElement( self, elementCoordIn, elLocalCoord, value );
	}
	
	return retValue;
}


double _FeVariable_GetMinGlobalFieldMagnitude( void* feVariable ) {
	FeVariable*			self = (FeVariable*)feVariable;
	FiniteElement_Mesh*		mesh = self->feMesh;
	int				node_lI=0;
	int				nodeLocalCount = mesh->nodeLocalCount;
	double				min = 0;
	double				globalMin = 0;
	double				currValue;

	min = FeVariable_GetScalarAtNode( self, 0 );
	
	/* Find upper and lower bounds on this processor */
	for ( node_lI = 0 ; node_lI < nodeLocalCount ; node_lI++ ) {
		currValue = FeVariable_GetScalarAtNode( self, node_lI );
		if ( currValue < min ) {
			min = currValue;
		}	
	}

	/* Find upper and lower bounds on all processors */
	MPI_Allreduce( &min, &globalMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
	return globalMin;
}


double _FeVariable_GetMaxGlobalFieldMagnitude( void* feVariable ) {
	FeVariable*			self = (FeVariable*)feVariable;
	FiniteElement_Mesh*		mesh = self->feMesh;
	int				node_lI=0;
	int				nodeLocalCount = mesh->nodeLocalCount;
	double				max = 0;
	double				globalMax = 0;
	double				currValue;

	max = FeVariable_GetScalarAtNode( self, 0 );
	
	/* Find upper and lower bounds on this processor */
	for ( node_lI = 0 ; node_lI < nodeLocalCount ; node_lI++ ) {
		currValue = FeVariable_GetScalarAtNode( self, node_lI );
		if ( currValue > max ) {
			max = currValue;
		}	
	}

	/* Find upper and lower bounds on all processors */
	MPI_Allreduce( &max, &globalMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
	return globalMax;
}

void _FeVariable_GetMinAndMaxLocalCoords( void* feVariable, Coord min, Coord max ) {
	FeVariable*		self = (FeVariable*)feVariable;
	FiniteElement_Mesh*	mesh = self->feMesh;
	ElementLayout*		eLayout = mesh->layout->elementLayout;
	double*			currCoord = NULL;
	Dimension_Index		dim_I = 0;
	Node_LocalIndex		node_lI = 0;

	if ( True == eLayout->getStaticMinAndMaxLocalCoords( eLayout, min, max ) ) {
		return;
	}
	else {
		/* Make sure that the mesh is initialised */
		assert( mesh->isInitialised );

		/* No shortcut, so just loop through every node coord, and record min & max */
		for (dim_I = 0; dim_I < self->dim ; dim_I++ ) {
			min[dim_I] = mesh->nodeCoord[0][dim_I];
			max[dim_I] = mesh->nodeCoord[0][dim_I];
		}

		for ( node_lI = 1; node_lI < mesh->nodeLocalCount; node_lI++ ) {
			currCoord = mesh->nodeCoord[node_lI];

			for (dim_I = 0; dim_I < self->dim ; dim_I++ ) {
				if ( currCoord[dim_I] < min[dim_I] ) {
					min[dim_I] = currCoord[dim_I];
				}
				else if ( currCoord[dim_I] > max[dim_I] ) {
					max[dim_I] = currCoord[dim_I];
				}
			}
		}
	}
}


void _FeVariable_GetMinAndMaxGlobalCoords( void* feVariable, Coord min, Coord max ) {
	FeVariable*		self = (FeVariable*)feVariable;
	FiniteElement_Mesh*	mesh = self->feMesh;
	ElementLayout*		eLayout = mesh->layout->elementLayout;

	if ( True == eLayout->getStaticMinAndMaxGlobalCoords( eLayout, min, max ) ) {
		return;
	}
	else {
		/* Use the default approach of each processor calculating min, then reduce */
		_FieldVariable_GetMinAndMaxGlobalCoords( feVariable, min, max );
	}
}

double FeVariable_GetScalarAtNode( void* feVariable, Node_LocalIndex lNode_I ) {
	FeVariable*	self = (FeVariable*)feVariable;
	Dof_Index	dofCountThisNode = 0;
	Dof_Index	nodeLocalDof_I = 0;
	double      value[ MAX_FIELD_COMPONENTS ];
	
	if ( self->dofLayout )
		dofCountThisNode = self->dofLayout->dofCounts[lNode_I];
	else 
		dofCountThisNode = self->fieldComponentCount;

	FeVariable_GetValueAtNode( self, lNode_I, value );

	if ( dofCountThisNode > 1) {
		double		magnitude=0;
		for ( nodeLocalDof_I=0; nodeLocalDof_I < dofCountThisNode; nodeLocalDof_I++ ) {
			magnitude += value[ nodeLocalDof_I ] * value[ nodeLocalDof_I ];
		}
		return sqrt( magnitude );
	}
	else 
		return value[0];

}	


void _FeVariable_GetValueAtNode( void* feVariable, Node_DomainIndex dNode_I, double* value ) {
	FeVariable*	self = (FeVariable*)feVariable;
	Variable*	currVariable = NULL;
	Dof_Index	dofCountThisNode = 0;
	Dof_Index	nodeLocalDof_I = 0;

	dofCountThisNode = self->dofLayout->dofCounts[dNode_I];
	
	for ( nodeLocalDof_I=0; nodeLocalDof_I < dofCountThisNode; nodeLocalDof_I++ ) {
		currVariable = DofLayout_GetVariable( self->dofLayout, dNode_I, nodeLocalDof_I );
		value[ nodeLocalDof_I ] = Variable_GetValueDouble( currVariable, dNode_I );
	}
}

void FeVariable_ZeroField( void* feVariable ) {
	FeVariable* self = (FeVariable*) feVariable;
	double*     values =  Memory_Alloc_Array( double, self->fieldComponentCount, "tempValues" );
	Index       lNode_I;

	memset( values, 0, self->fieldComponentCount * sizeof(double) );

	for( lNode_I = 0 ; lNode_I < self->feMesh->nodeLocalCount ; lNode_I++ ) {
		FeVariable_SetValueAtNode( self, lNode_I, values );
	}

	Memory_Free( values );
}
	
/* --- Public Functions --- */

InterpolationResult FeVariable_GetElementLocalCoordAtGlobalCoord( void* feVariable, Coord globalCoord, Coord elLocalCoord,
		Element_DomainIndex* elementCoordInPtr )
{
	FeVariable*		self = (FeVariable*)feVariable;
	MeshLayout*		mLayout = self->feMesh->layout;
	ElementLayout*		eLayout = mLayout->elementLayout;
	InterpolationResult	retValue;

	(*elementCoordInPtr) = (unsigned)-1;
	
	/* locate which mesh element given coord is in : use inclusive upper boundaries to save
		the need to use shadow space if possible */
	(*elementCoordInPtr) = Mesh_ElementWithPoint( self->feMesh, globalCoord, INCLUSIVE_UPPER_BOUNDARY );
#if 0
	if( eLayout->type == ParallelPipedHexaEL_Type ) {
		(*elementCoordInPtr) = eLayout->elementWithPoint( eLayout, mLayout->decomp, globalCoord,
							    INCLUSIVE_UPPER_BOUNDARY, 0, NULL );
	}
	else {
		unsigned	cNode;

		/* Find closest node to point. */
		cNode = _FeVariable_ClosestNode( self, globalCoord );

		/* Find with hint of incident elements. */
		(*elementCoordInPtr) = eLayout->elementWithPoint( eLayout, mLayout->decomp, globalCoord,
							    INCLUSIVE_UPPER_BOUNDARY, 
							    self->feMesh->nodeElementCountTbl[cNode], self->feMesh->nodeElementTbl[cNode] );

		/* If still no cigar, brute force. */
		if ( (*elementCoordInPtr) >= self->feMesh->elementDomainCount ) {
			(*elementCoordInPtr) = eLayout->elementWithPoint( eLayout, mLayout->decomp, globalCoord,
								    INCLUSIVE_UPPER_BOUNDARY, 0, NULL );
		}
	}
#endif

	if ( (*elementCoordInPtr) >= self->feMesh->elementDomainCount ) {
		Bool			outsideGlobal = False;
		Coord			min, max;
		Dimension_Index		dim_I=0;
		Bool			checkResult;
		Stream*			errorStr = Journal_Register( Error_Type, self->type );

		checkResult = ElementLayout_GetStaticMinAndMaxGlobalCoords( mLayout->elementLayout, min, max );
		Journal_Firewall( True == checkResult, errorStr, "Error - in %s: current mesh doesn't"
			"know how to calculate global max and min values.\n", __func__ );
		
		for ( dim_I = 0; dim_I < self->dim; dim_I++ ) {
			if ( ( globalCoord[dim_I] < min[dim_I] ) || (globalCoord[dim_I] > max[dim_I] ) ) {
				outsideGlobal = True;
			}
		}

		if ( outsideGlobal == True ) {
			return OUTSIDE_GLOBAL;
		}
		else {
			return OTHER_PROC;
		}	
	}	
	else /* We found the coord is within a local or shadow element */ {
		Node_LocalIndex		currElementNodeCount=0;
		Coord**			globalNodeCoordPtrs=NULL;
		ElementType*		elementType = NULL;
	
		if ( (*elementCoordInPtr) < self->feMesh->elementLocalCount ) {
			retValue = LOCAL;
		}
		else {
			retValue = SHADOW;
		}

		/* convert global coordinate to local co-ordinates of element the coord is in */
		currElementNodeCount = self->feMesh->elementNodeCountTbl[(*elementCoordInPtr)];
		globalNodeCoordPtrs = Memory_Alloc_Array( Coord*, currElementNodeCount, "globalNodeCoordPtrs" );
		Mesh_GetNodeCoordPtrsOfElement( self->feMesh, (*elementCoordInPtr), globalNodeCoordPtrs );

		elementType = FeMesh_ElementTypeAt( self->feMesh, (*elementCoordInPtr) );
		ElementType_ConvertGlobalCoordToElLocal( elementType, eLayout,
			(const Coord**) globalNodeCoordPtrs, globalCoord, elLocalCoord );

		Memory_Free( globalNodeCoordPtrs );
	}	

	return retValue;
}


void FeVariable_SetValueAtNode( void* feVariable, Node_DomainIndex dNode_I, double* componentValues ) {
	FeVariable*	self = (FeVariable*)feVariable;
	Dof_Index	dofCountThisNode = 0;
	Dof_Index	nodeLocalDof_I = 0;

	dofCountThisNode = self->dofLayout->dofCounts[dNode_I];
	
	for ( nodeLocalDof_I=0; nodeLocalDof_I < dofCountThisNode; nodeLocalDof_I++ ) {
		DofLayout_SetValueDouble( (self)->dofLayout, dNode_I, nodeLocalDof_I, componentValues[nodeLocalDof_I] );
	}
}


void FeVariable_PrintLocalDiscreteValues_2dBox( void* variable, Stream* stream ) {
	FeVariable*			self = (FeVariable*)variable;
	FiniteElement_Mesh*		mesh = self->feMesh;
	ParallelPipedHexaEL*		elementLayout = (ParallelPipedHexaEL*)mesh->layout->elementLayout;
	Node_LocalIndex			node_lI=0;
	Index				x_I, y_I;
	Index				ii;
	Dof_Index			dof_I=0;
	Dof_Index			currNodeNumDofs=0;
	BlockGeometry*			bGeometry = (BlockGeometry*)elementLayout->geometry;
	Index				nx = 0;
	Index				ny = 0;
	double				dx = 0;
	double				dy = 0;
	DofLayout*			dofLayout = self->dofLayout;
	Stream*				eStream = Journal_Register( Error_Type, self->type );
	HexaMD*				hexaMD = (HexaMD*)mesh->layout->decomp;
	Index				minLocalNodeX;
	Index				minLocalNodeY;
	Index				maxLocalNodeX;
	Index				maxLocalNodeY;


	if ( elementLayout->type != ParallelPipedHexaEL_Type && elementLayout->dim == 2 ) {
		Journal_Printf( eStream, "Warning: %s called on variable \"%s\", but this isn't stored on a "
			"regular 2D (%s) mesh - so just returning.\n", __func__, self->name, ParallelPipedHexaEL_Type );
		return;
	}	
	
	nx = bGeometry->size[I_AXIS];
	ny = bGeometry->size[J_AXIS];
	dx = elementLayout->elementLengthEachDim[I_AXIS];
	dy = elementLayout->elementLengthEachDim[J_AXIS];

	minLocalNodeX = hexaMD->_nodeOffsets[hexaMD->rank][I_AXIS];
	minLocalNodeY = hexaMD->_nodeOffsets[hexaMD->rank][J_AXIS];
	maxLocalNodeX = minLocalNodeX + hexaMD->nodeLocal3DCounts[hexaMD->rank][I_AXIS];
	maxLocalNodeY = minLocalNodeY + hexaMD->nodeLocal3DCounts[hexaMD->rank][J_AXIS];

	Journal_Printf( stream, "display of Values in 2D box X:{%5.2f-%5.2f}, Y:{%5.2f-%5.2f}\n",
		bGeometry->min[I_AXIS], bGeometry->max[I_AXIS],
		bGeometry->min[J_AXIS], bGeometry->max[J_AXIS] );
	Journal_Printf( stream, "\twith %d elements in X (dx=%5.2f) and %d elements in Y (dy=%5.2f)\n\n",
		nx-1, dx, ny-1, dy );

	/*Header*/
	for (ii=0;ii<10;ii++) Journal_Printf( stream, " " );
	for ( x_I=0; x_I < nx; x_I++ ) {
		Journal_Printf( stream, "|  xNode=%3d   ", x_I );
	}
	Journal_Printf( stream, "|\n", x_I );

	for ( y_I= ny-1; y_I != (unsigned)-1; y_I-- ) {
		/*Blocks */
		for (ii=0;ii<10;ii++) Journal_Printf( stream, " " );
		for ( x_I=0; x_I < nx; x_I++ ) {
			if (y_I == ny-1) {
				Journal_Printf( stream, "-" );
			}
			else if (x_I==0) {
				Journal_Printf( stream, "|" );
			}	
			else {
				Journal_Printf( stream, "*" );
			}	
			for (ii=0;ii<14;ii++) Journal_Printf( stream, "-" );
		}
		if (y_I == ny-1) {
			Journal_Printf( stream, "-\n" );
		}	
		else {
			Journal_Printf( stream, "|\n" );
		}	


		/* Now a row of y values */
		Journal_Printf( stream, "yNode=%3d |", y_I );
		for ( x_I=0; x_I < nx; x_I++ ) {
		
			if ( ( y_I >= minLocalNodeY ) && ( y_I < maxLocalNodeY )
				&& ( x_I >= minLocalNodeX ) && ( x_I < maxLocalNodeX ) ) {

				node_lI = RegularMeshUtils_Node_Global3DToLocal1D( hexaMD, x_I, y_I, 0 );
				currNodeNumDofs = dofLayout->dofCounts[node_lI];

				if ( currNodeNumDofs == 1 ) {
					Journal_Printf( stream, "   " );
				}
				Journal_Printf( stream, "(" );
				for ( dof_I=0; dof_I < currNodeNumDofs - 1 ; dof_I++ ) {
					Journal_Printf( stream, "%5.2f,", DofLayout_GetValueDouble( dofLayout, node_lI, dof_I ) );
				}
				Journal_Printf( stream, "%5.2f )", DofLayout_GetValueDouble( dofLayout, node_lI, dof_I ) );
				
				if ( currNodeNumDofs == 1 ) {
					Journal_Printf( stream, "   " );
				}
				Journal_Printf( stream, "|" );
			}
			else {
				for (ii=0;ii<14;ii++) Journal_Printf( stream, "X" );
				Journal_Printf( stream, "|" );
			}
		}
		Journal_Printf( stream, "\n" );
	}
	
	/*Blocks */
	for (ii=0;ii<10;ii++) Journal_Printf( stream, " " );
	for ( x_I=0; x_I < nx; x_I++ ) {
		Journal_Printf( stream, "-" );
		for (ii=0;ii<14;ii++) Journal_Printf( stream, "-" );
	}
	Journal_Printf( stream, "-\n", x_I );
}


Bool FeVariable_InterpolateDerivativesAt( void* variable, Coord globalCoord, double* value ) {
	FeVariable*	        self                 = (FeVariable*)variable;
	MeshLayout*	        mLayout              = self->feMesh->layout;
	ElementLayout*	    eLayout              = mLayout->elementLayout;
	Element_DomainIndex	elementCoordIn       = (unsigned)-1;
	Node_LocalIndex	    currElementNodeCount = 0;
	Coord**			    globalNodeCoordPtrs  = NULL;
	ElementType*		elementType          = NULL;
	Coord               elLocalCoord         = {0,0,0};

	/* Need a special rule for points on this processor's boundary: instead of the normal
	   rule, "round" the point to lie inside the local space, rather than shadow */
	
	/* locate which mesh element given coord is in : use inclusive upper boundaries to save
		the need to use shadow space if possible */
	elementCoordIn = Mesh_ElementWithPoint( self->feMesh, globalCoord, INCLUSIVE_UPPER_BOUNDARY );

	if ( elementCoordIn >= self->feMesh->elementDomainCount ) {
		/* If coord isn't inside domain elements list, bail out */
		return False;
	}	
	else /* We found the coord is within a local or shadow element */ {
		if ( elementCoordIn >= self->feMesh->elementLocalCount ) {
			if ( False == self->shadowValuesSynchronised ) {
				Stream* warningStr = Journal_Register( Error_Type, self->type );
				Journal_Printf( warningStr, "Warning - in %s: user asking to interpolate derivatives "
					"to coord (%g,%g,%g), which is in shadow space, but "
					"FeVariable_SyncShadowValues() hasn't been called yet.\n", 
					__func__, globalCoord[0], globalCoord[1], globalCoord[2] );
				return False;	
			}
		}
		/* convert global coordinate to local co-ordinates of element the coord is in */
		currElementNodeCount = self->feMesh->elementNodeCountTbl[elementCoordIn];
		globalNodeCoordPtrs = Memory_Alloc_Array( Coord*, currElementNodeCount, "globalNodeCoordPtrs" );
		Mesh_GetNodeCoordPtrsOfElement( self->feMesh, elementCoordIn, globalNodeCoordPtrs );

		elementType = FeMesh_ElementTypeAt( self->feMesh, elementCoordIn );
		ElementType_ConvertGlobalCoordToElLocal( elementType, eLayout,
			(const Coord**) globalNodeCoordPtrs, globalCoord, elLocalCoord );

		/* Now interpolate the value at that coordinate, using shape functions */
		FeVariable_InterpolateDerivativesToElLocalCoord( self, elementCoordIn, elLocalCoord, value );
		Memory_Free( globalNodeCoordPtrs );
	}	
	
	return True;
}

void FeVariable_InterpolateDerivativesToElLocalCoord( void* _feVariable, Element_DomainIndex lElement_I, Coord elLocalCoord, double* value ) {
	FeVariable*    self             = (FeVariable*) _feVariable;
	ElementType*            elementType      = FeMesh_ElementTypeAt( self->feMesh, lElement_I );
	Node_Index              elementNodeCount = elementType->nodeCount;
	double**                GNx; 
	double                  detJac;
	Dimension_Index         dim         = self->dim;

	GNx = Memory_Alloc_2DArray( double, dim, elementNodeCount, "Global Shape Function Derivatives" );

	/* Evaluate Global Shape Functions */
	ElementType_ShapeFunctionsGlobalDerivs( 
			elementType,
			self->feMesh, lElement_I,
			elLocalCoord, dim, &detJac, GNx );

	/* Do Interpolation */
	FeVariable_InterpolateDerivatives_WithGNx( self, lElement_I, GNx, value );
	
	Memory_Free(GNx);
}

void FeVariable_InterpolateDerivatives_WithGNx( void* _feVariable, Element_LocalIndex lElement_I, double** GNx, double* value ) {
	FeVariable*             self        = (FeVariable*) _feVariable;
	ElementType*            elementType = FeMesh_ElementTypeAt( self->feMesh, lElement_I );
	Node_ElementLocalIndex  elLocalNode_I;
	Node_ElementLocalIndex  elNodeCount = elementType->nodeCount;
	Node_LocalIndex         lNode_I;
	Dof_Index               dof_I;
	Dof_Index               dofCount;
	Variable*               dofVariable;
	double                  nodeValue;
	Dimension_Index         dim         = self->dim;

	/* Gets number of degrees of freedom - assuming it is the same throughout the mesh */
	dofCount = self->dofLayout->dofCounts[0];

	/* Initialise */
	memset( value, 0, sizeof( double ) * dofCount * dim );

	for ( dof_I = 0 ; dof_I < dofCount ; dof_I++ ) {
		/* Interpolate derivative from nodes */
		for ( elLocalNode_I = 0 ; elLocalNode_I < elNodeCount ; elLocalNode_I++) {
			lNode_I      = self->feMesh->elementNodeTbl[ lElement_I ][ elLocalNode_I ];
			dofVariable  = DofLayout_GetVariable( self->dofLayout, lNode_I, dof_I );
			nodeValue    = Variable_GetValueDouble( dofVariable, lNode_I );
			
			value[dof_I*dim + 0] += GNx[0][elLocalNode_I] * nodeValue;
			value[dof_I*dim + 1] += GNx[1][elLocalNode_I] * nodeValue;
			if( dim == 3 ) 
				value[dof_I*dim + 2] += GNx[2][elLocalNode_I] * nodeValue;	
		}
	}
}
	

void FeVariable_GetMinimumSeparation( void* feVariable, double* minSeparationPtr, double minSeparationEachDim[3] )
{
	FeVariable*            self = (FeVariable*) feVariable;
	ParallelPipedHexaEL*   eLayout = (ParallelPipedHexaEL*) self->feMesh->layout->elementLayout;
	FiniteElement_Mesh*    mesh = self->feMesh;

	if ( Stg_Class_IsInstance( eLayout, ParallelPipedHexaEL_Type ) ) {
		(*minSeparationPtr) = eLayout->elementLengthEachDim[I_AXIS];
		minSeparationEachDim[I_AXIS] = eLayout->elementLengthEachDim[I_AXIS];
		
		if ( eLayout->elementLengthEachDim[J_AXIS] < *minSeparationPtr ) {
			*minSeparationPtr = eLayout->elementLengthEachDim[J_AXIS];
		}
		minSeparationEachDim[J_AXIS] = eLayout->elementLengthEachDim[J_AXIS];

		if ( self->dim == 3 ) {
			if ( eLayout->elementLengthEachDim[K_AXIS] < *minSeparationPtr ) {
				*minSeparationPtr = eLayout->elementLengthEachDim[K_AXIS];
			}
			minSeparationEachDim[K_AXIS] = eLayout->elementLengthEachDim[K_AXIS];
		}
		else {
			minSeparationEachDim[K_AXIS] = 0;
		}
	}
	else if ( Stg_Class_IsInstance( eLayout, HexaEL_Type ) ) {
		double              currElementSeparation[3] = {0,0,0};
		double              currPairSeparation[3] = {0,0,0};
		Element_LocalIndex  element_lI = 0;

		minSeparationEachDim[I_AXIS] = HUGE_VAL;
		minSeparationEachDim[J_AXIS] = HUGE_VAL;
		minSeparationEachDim[K_AXIS] = HUGE_VAL;
		*minSeparationPtr = HUGE_VAL;

		for ( element_lI = 0; element_lI < mesh->elementLocalCount; element_lI++ ) {
			/* Axis 0 (I) */
			currPairSeparation[0] = mesh->nodeCoord[mesh->elementNodeTbl[element_lI][1]][0]
				- mesh->nodeCoord[mesh->elementNodeTbl[element_lI][0]][0];
			currElementSeparation[0] = currPairSeparation[0];

			currPairSeparation[0] = mesh->nodeCoord[mesh->elementNodeTbl[element_lI][2]][0]
				- mesh->nodeCoord[mesh->elementNodeTbl[element_lI][3]][0];
			if ( currPairSeparation[0] < currElementSeparation[0] ) {
				currElementSeparation[0] = currPairSeparation[0];
			}
			
			if ( self->dim == 3 ) {
				currPairSeparation[0] = mesh->nodeCoord[mesh->elementNodeTbl[element_lI][5]][0]
					- mesh->nodeCoord[mesh->elementNodeTbl[element_lI][4]][0];
				if ( currPairSeparation[0] < currElementSeparation[0] ) {
					currElementSeparation[0] = currPairSeparation[0];
				}
				currPairSeparation[0] = mesh->nodeCoord[mesh->elementNodeTbl[element_lI][6]][0]
					- mesh->nodeCoord[mesh->elementNodeTbl[element_lI][7]][0];
				if ( currPairSeparation[0] < currElementSeparation[0] ) {
					currElementSeparation[0] = currPairSeparation[0];
				}
			}
			if ( currElementSeparation[0] < minSeparationEachDim[0] ) {
				minSeparationEachDim[0] = currElementSeparation[0];
			}
			/* Axis 1 (J) */
			currPairSeparation[1] = mesh->nodeCoord[mesh->elementNodeTbl[element_lI][3]][1]
				- mesh->nodeCoord[mesh->elementNodeTbl[element_lI][0]][1];
			currElementSeparation[1] = currPairSeparation[1];

			currPairSeparation[1] = mesh->nodeCoord[mesh->elementNodeTbl[element_lI][2]][1]
				- mesh->nodeCoord[mesh->elementNodeTbl[element_lI][1]][1];
			if ( currPairSeparation[1] < currElementSeparation[1] ) {
				currElementSeparation[1] = currPairSeparation[1];
			}
			
			if ( self->dim == 3 ) {
				currPairSeparation[1] = mesh->nodeCoord[mesh->elementNodeTbl[element_lI][7]][1]
					- mesh->nodeCoord[mesh->elementNodeTbl[element_lI][4]][1];
				if ( currPairSeparation[1] < currElementSeparation[1] ) {
					currElementSeparation[1] = currPairSeparation[1];
				}
				currPairSeparation[1] = mesh->nodeCoord[mesh->elementNodeTbl[element_lI][6]][1]
					- mesh->nodeCoord[mesh->elementNodeTbl[element_lI][5]][1];
				if ( currPairSeparation[1] < currElementSeparation[1] ) {
					currElementSeparation[1] = currPairSeparation[1];
				}
			}
			if ( currElementSeparation[1] < minSeparationEachDim[1] ) {
				minSeparationEachDim[1] = currElementSeparation[1];
			}
			/* Axis 2 (K) */
			if ( self->dim == 3 ) {
				currPairSeparation[2] = mesh->nodeCoord[mesh->elementNodeTbl[element_lI][4]][2]
					- mesh->nodeCoord[mesh->elementNodeTbl[element_lI][0]][2];
				currElementSeparation[2] = currPairSeparation[2];

				currPairSeparation[2] = mesh->nodeCoord[mesh->elementNodeTbl[element_lI][5]][2]
					- mesh->nodeCoord[mesh->elementNodeTbl[element_lI][1]][2];
				if ( currPairSeparation[2] < currElementSeparation[2] ) {
					currElementSeparation[2] = currPairSeparation[2];
				}
			
				currPairSeparation[2] = mesh->nodeCoord[mesh->elementNodeTbl[element_lI][7]][2]
					- mesh->nodeCoord[mesh->elementNodeTbl[element_lI][3]][2];
				if ( currPairSeparation[2] < currElementSeparation[2] ) {
					currElementSeparation[2] = currPairSeparation[2];
				}
				currPairSeparation[2] = mesh->nodeCoord[mesh->elementNodeTbl[element_lI][6]][2]
					- mesh->nodeCoord[mesh->elementNodeTbl[element_lI][2]][2];
				if ( currPairSeparation[2] < currElementSeparation[2] ) {
					currElementSeparation[2] = currPairSeparation[2];
				}
				if ( currElementSeparation[2] < minSeparationEachDim[2] ) {
					minSeparationEachDim[2] = currElementSeparation[2];
				}
			}
			else {
				minSeparationEachDim[2] = 0;
			}
		}	
		*minSeparationPtr = minSeparationEachDim[0];
		*minSeparationPtr = minSeparationEachDim[1] < *minSeparationPtr ? minSeparationEachDim[1] : *minSeparationPtr;
		if ( self->dim == 3 ) {
			*minSeparationPtr = minSeparationEachDim[2] < *minSeparationPtr ? minSeparationEachDim[2] : *minSeparationPtr;
		}
	}
	else  {
		Stream*    errorStr = Journal_Register( Error_Type, self->type );
		Journal_Firewall( 0, errorStr, "Error: in %s - Don't know how to find minSeparation for element type %s.\n",
			__func__, eLayout->type );
	}
}	


void FeVariable_SyncShadowValues( void* feVariable ) {
	FeVariable*			self = (FeVariable*)feVariable;
	Neighbour_Index			nbr_I = 0;
	Node_Index			node_stI = 0;
	Node_DomainIndex		node_dI = 0;
	Node_LocalIndex			node_lI = 0;
	Dof_Index			nodalDof_I = 0;
	Processor_Index			nbrRank = 0;
	Index*				incomingDofTotals = NULL;			
	double**			incomingDofValues = NULL;
	MPI_Request**			incomingDofValRequests = NULL;
	Node_Index			myShadowNodesOnThisNbrCount = 0;
	Dof_Index			incomingDof_I;
	Index*				outgoingDofTotals = NULL;			
	double**			outgoingDofValues = NULL;
	MPI_Request**			outgoingDofValRequests = NULL;
	Node_Index			nbrShadowNodesOnMeCount = 0;
	Dof_Index			outgoingDof_I;
	FiniteElement_Mesh*		mesh = self->feMesh;
	MPI_Status			status;
	int				incomingDofValueSetsYetToReceive;
	Bool*				incomingDofValueSetsReceived;

	Journal_DPrintf( self->debug, "In %s- for feVariable \"%s\":\n", __func__, self->name );
	Stream_Indent( self->debug );

	if ( ( 1 == mesh->layout->decomp->procsInUse ) || ( 0 == mesh->layout->decomp->shadowDepth ) ) {
		Journal_DPrintf( self->debug, "No shadow nodes: nothing to do - returning.\n" );
		Stream_UnIndent( self->debug );
		return;
	}

	self->shadowValuesSynchronised = True;

	/* allocate memory for incoming info */
	incomingDofTotals = Memory_Alloc_Array( Index, mesh->procNbrInfo->procNbrCnt, "incomingDofTotals" );
	incomingDofValues = Memory_Alloc_Array( double*, mesh->procNbrInfo->procNbrCnt, "incomingDofValues" );
	incomingDofValRequests = Memory_Alloc_Array( MPI_Request*, mesh->procNbrInfo->procNbrCnt, "incomingDofValRequests" );
	incomingDofValueSetsReceived = Memory_Alloc_Array( Bool, mesh->procNbrInfo->procNbrCnt, "incomingDofValueSets" );
	incomingDofValueSetsYetToReceive = 0;
	for ( nbr_I=0; nbr_I < mesh->procNbrInfo->procNbrCnt; nbr_I++ ) {
		myShadowNodesOnThisNbrCount = mesh->nodeShadowInfo->procShadowCnt[nbr_I];
		incomingDofTotals[nbr_I] = 0;
		incomingDofValues[nbr_I] = NULL;
		incomingDofValRequests[nbr_I] = NULL;
		incomingDofValueSetsReceived[nbr_I] = False;

		if ( myShadowNodesOnThisNbrCount > 0 ) {
			for( node_stI = 0; node_stI < myShadowNodesOnThisNbrCount; node_stI++ ) {
				node_dI = mesh->nodeShadowInfo->procShadowTbl[nbr_I][node_stI];
				incomingDofTotals[nbr_I] += self->dofLayout->dofCounts[node_dI];
			}	
			incomingDofValues[nbr_I] = Memory_Alloc_Array( double, incomingDofTotals[nbr_I],
				"incomingDofValues[]" );
			incomingDofValRequests[nbr_I] = Memory_Alloc( MPI_Request, "incomingDofValRequest" );	
			incomingDofValueSetsYetToReceive++;
		}
	}
	/* allocate memory for outgoing info */
	outgoingDofTotals = Memory_Alloc_Array( Index, mesh->procNbrInfo->procNbrCnt, "outgoingDofTotals" );
	outgoingDofValues = Memory_Alloc_Array( double*, mesh->procNbrInfo->procNbrCnt, "outgoingDofValues" );
	outgoingDofValRequests = Memory_Alloc_Array( MPI_Request*, mesh->procNbrInfo->procNbrCnt, "outgoingDofValRequests" );
	for ( nbr_I=0; nbr_I < mesh->procNbrInfo->procNbrCnt; nbr_I++ ) {
		nbrShadowNodesOnMeCount = mesh->nodeShadowInfo->procShadowedCnt[nbr_I];
		outgoingDofTotals[nbr_I] = 0;
		outgoingDofValues[nbr_I] = NULL;
		outgoingDofValRequests[nbr_I] = NULL;

		if ( nbrShadowNodesOnMeCount > 0 ) {
			for( node_stI = 0; node_stI < nbrShadowNodesOnMeCount; node_stI++ ) {
				node_lI = mesh->nodeShadowInfo->procShadowedTbl[nbr_I][node_stI];
				outgoingDofTotals[nbr_I] += self->dofLayout->dofCounts[node_lI];
			}	
			outgoingDofValues[nbr_I] = Memory_Alloc_Array( double, outgoingDofTotals[nbr_I],
				"outgoingDofValues[]" );
			outgoingDofValRequests[nbr_I] = Memory_Alloc( MPI_Request, "outgoingDofValRequest" );	
		}
	}	

	Journal_DPrintfL( self->debug, 2, "Starting non-blocking recv's of incoming values\n" );
	Stream_Indent( self->debug );
	for ( nbr_I=0; nbr_I < mesh->procNbrInfo->procNbrCnt; nbr_I++ ) {
		myShadowNodesOnThisNbrCount = mesh->nodeShadowInfo->procShadowCnt[nbr_I];
		if ( myShadowNodesOnThisNbrCount > 0 ) {
			nbrRank = mesh->procNbrInfo->procNbrTbl[nbr_I];
			Journal_DPrintfL( self->debug, 2, "Start recv from proc %u - %u values\n", nbrRank, incomingDofTotals[nbr_I] );
			MPI_Irecv( incomingDofValues[nbr_I], incomingDofTotals[nbr_I], MPI_DOUBLE, nbrRank,
				DOF_VALUES_TAG, self->communicator, incomingDofValRequests[nbr_I] );
		}
	}
	Stream_UnIndent( self->debug );

	Journal_DPrintfL( self->debug, 2, "Non-blocking send out required shadow values to neighbours\n" );
	Stream_Indent( self->debug );
	for ( nbr_I=0; nbr_I < mesh->procNbrInfo->procNbrCnt; nbr_I++ ) {
		if ( mesh->nodeShadowInfo->procShadowedCnt[nbr_I] > 0 ) {
			nbrRank = mesh->procNbrInfo->procNbrTbl[nbr_I];

			outgoingDof_I = 0;
			for( node_stI = 0; node_stI < mesh->nodeShadowInfo->procShadowedCnt[nbr_I]; node_stI++ ) {
				node_lI = mesh->nodeShadowInfo->procShadowedTbl[nbr_I][node_stI];
				for ( nodalDof_I=0; nodalDof_I < self->dofLayout->dofCounts[node_lI]; nodalDof_I++ ) {
					outgoingDofValues[nbr_I][outgoingDof_I] =
						DofLayout_GetValueDouble( self->dofLayout, node_lI, nodalDof_I );
					outgoingDof_I++;
				}
			}
			Journal_DPrintfL( self->debug, 2, "Start send to proc %u - %u values\n", nbrRank, outgoingDofTotals[nbr_I] );
			MPI_Isend( outgoingDofValues[nbr_I], outgoingDofTotals[nbr_I], MPI_DOUBLE, nbrRank,
				DOF_VALUES_TAG, self->communicator, outgoingDofValRequests[nbr_I] );
		}
	}
	Stream_UnIndent( self->debug );

	Journal_DPrintfL( self->debug, 2, "Receiving and updating shadow values I need from neighbours:\n" );
	Stream_Indent( self->debug );
	while ( incomingDofValueSetsYetToReceive > 0 ) {
		int	testFlag = 0;

		for ( nbr_I=0; nbr_I < mesh->procNbrInfo->procNbrCnt; nbr_I++ ) {

			if ( ( mesh->nodeShadowInfo->procShadowCnt[nbr_I] > 0 )
				&& ( False == incomingDofValueSetsReceived[nbr_I] ) )
			{
				MPI_Test( incomingDofValRequests[nbr_I], &testFlag, &status );
				if ( False == testFlag ) {
					continue;
				}
				else {
					Journal_DPrintfL( self->debug, 2, "Recv'd a batch of values from proc %u: updating...\n",
						nbrRank );
					/* update the appropriate values from recv'd set */
					incomingDof_I = 0;
					for( node_stI = 0; node_stI < mesh->nodeShadowInfo->procShadowCnt[nbr_I]; node_stI++ ) {
						node_dI = mesh->nodeShadowInfo->procShadowTbl[nbr_I][node_stI];
						for ( nodalDof_I=0; nodalDof_I < self->dofLayout->dofCounts[node_dI]; nodalDof_I++ ) {
							DofLayout_SetValueDouble( self->dofLayout, node_dI, nodalDof_I,
								incomingDofValues[nbr_I][incomingDof_I] );
							incomingDof_I++;	
						}
					}
					incomingDofValueSetsReceived[nbr_I] = True;
					incomingDofValueSetsYetToReceive--;
				}
			}	
		}	
	}
	Journal_DPrintfL( self->debug, 2, "Done.\n" );
	Stream_UnIndent( self->debug );

	Journal_DPrintfL( self->debug, 2, "Making sure outgoing sends have completed...\n" );
	for ( nbr_I=0; nbr_I < mesh->procNbrInfo->procNbrCnt; nbr_I++ ) {
		if ( mesh->nodeShadowInfo->procShadowCnt[nbr_I] > 0 ) {
			MPI_Wait( outgoingDofValRequests[nbr_I], &status );
		}
	}
	Journal_DPrintfL( self->debug, 2, "Done.\n" );

	/* clean up temporary memory */
	for ( nbr_I=0; nbr_I < mesh->procNbrInfo->procNbrCnt; nbr_I++ ) {
		myShadowNodesOnThisNbrCount = mesh->nodeShadowInfo->procShadowCnt[nbr_I];
		if ( myShadowNodesOnThisNbrCount > 0 ) {
			Memory_Free( incomingDofValues[nbr_I] );
			Memory_Free( incomingDofValRequests[nbr_I] );
		}	
	}
	Memory_Free( incomingDofTotals );
	Memory_Free( incomingDofValues );
	Memory_Free( incomingDofValRequests );
	for ( nbr_I=0; nbr_I < mesh->procNbrInfo->procNbrCnt; nbr_I++ ) {
		nbrShadowNodesOnMeCount = mesh->nodeShadowInfo->procShadowedCnt[nbr_I];
		if ( nbrShadowNodesOnMeCount > 0 ) {
			Memory_Free( outgoingDofValues[nbr_I] );
			Memory_Free( outgoingDofValRequests[nbr_I] );
		}	
	}
	Memory_Free( outgoingDofTotals );
	Memory_Free( outgoingDofValues );
	Memory_Free( outgoingDofValRequests );
	
	Stream_UnIndent( self->debug );
}


void FeVariable_PrintDomainDiscreteValues( void* variable, Stream* stream ) {
	FeVariable* self = (FeVariable*)variable;

	Journal_Printf( stream, "In %s: for FeVariable \"%s\":\n", __func__, self->name );

	_FeVariable_PrintLocalOrDomainValues( variable, self->feMesh->nodeDomainCount, stream );
}

void FeVariable_PrintCoordsAndValues( void* _feVariable, Stream* stream ) {
	FeVariable*         self            = (FeVariable*) _feVariable;
	FiniteElement_Mesh* mesh            = self->feMesh;
	Node_LocalIndex     node_I          = 0;
	Node_LocalIndex     nodeLocalCount  = mesh->nodeLocalCount;
	Dof_Index           currNodeNumDofs;
	Dof_Index           nodeLocalDof_I;
	Variable*           currVariable;
	double*             nodeCoord;
	
	/* Print Header of stream */
	Journal_Printf( stream, "# FeVariable - %s\n", self->name );
	Journal_Printf( stream, "#    x coord   |    y coord   |    z coord" );
	currNodeNumDofs = self->dofLayout->dofCounts[ 0 ];
	for ( nodeLocalDof_I = 0; nodeLocalDof_I < currNodeNumDofs; nodeLocalDof_I++ ) {
		currVariable = DofLayout_GetVariable( self->dofLayout, node_I, nodeLocalDof_I );
		Journal_Printf( stream, "  |  %s", currVariable->name );
	}
	Journal_Printf(stream, "\n");
	
	/* Loop over local nodes */
	for( node_I=0; node_I < nodeLocalCount ; node_I++ ) {
		currNodeNumDofs = self->dofLayout->dofCounts[ node_I ];

		/* Get Coordinate of Node */
		nodeCoord = Mesh_CoordAt( mesh, node_I );
		Journal_Printf( stream, "%12.6g   %12.6g   %12.6g   ", 
				nodeCoord[ I_AXIS ], nodeCoord[ J_AXIS ], nodeCoord[ K_AXIS ] );
		
		/* Print each dof */
		for ( nodeLocalDof_I = 0; nodeLocalDof_I < currNodeNumDofs; nodeLocalDof_I++ ) {
			currVariable = DofLayout_GetVariable( self->dofLayout, node_I, nodeLocalDof_I );
			Journal_Printf( stream, "%12.6g   ", Variable_GetValueDouble( currVariable, node_I ) );
		}
		Journal_Printf( stream, "\n" );
	}	
}


/* --- Private Functions --- */

void _FeVariable_InterpolateNodeValuesToElLocalCoord( void* feVariable, Element_DomainIndex element_lI, Coord elLocalCoord, double* value ) {
	FeVariable*         self       = (FeVariable*) feVariable;
	ElementType*		elementType=NULL;
	Dof_Index		nodeLocalDof_I=0;
	Dof_Index		dofCountThisNode=0;
	Node_ElementLocalIndex	elLocalNode_I=0;
	Node_ElementLocalIndex	nodeCountThisElement=0;
	double*			shapeFuncsEvaluatedAtCoord=NULL;
	Node_LocalIndex		lNode_I=0;
	Variable*		currVariable=NULL;
	double			dofValueAtCurrNode=0;

	nodeCountThisElement = self->feMesh->elementNodeCountTbl[element_lI];
	/* Gets number of degrees of freedom - assuming it is the same throughout the mesh */
	dofCountThisNode = self->dofLayout->dofCounts[lNode_I];

	/* evaluate shape function values of current element at elLocalCoords */
	elementType = FeMesh_ElementTypeAt( self->feMesh, element_lI );
	shapeFuncsEvaluatedAtCoord = Memory_Alloc_Array_Unnamed( double, elementType->nodeCount );
	ElementType_EvaluateShapeFunctionsAt( elementType, elLocalCoord, shapeFuncsEvaluatedAtCoord );

	for ( nodeLocalDof_I=0; nodeLocalDof_I < dofCountThisNode; nodeLocalDof_I++ ) {
		value[nodeLocalDof_I] = 0;
	}

	/* Now for each node, add that node's contribution at point */
	for ( elLocalNode_I=0; elLocalNode_I < nodeCountThisElement; elLocalNode_I++ ) {
		lNode_I = self->feMesh->elementNodeTbl[element_lI][elLocalNode_I];
		
		for ( nodeLocalDof_I=0; nodeLocalDof_I < dofCountThisNode; nodeLocalDof_I++ ) {
			currVariable = DofLayout_GetVariable( self->dofLayout, lNode_I, nodeLocalDof_I );
			dofValueAtCurrNode = Variable_GetValueDouble( currVariable, lNode_I );
			value[nodeLocalDof_I] += dofValueAtCurrNode * shapeFuncsEvaluatedAtCoord[elLocalNode_I];
		}	
	}
	Memory_Free( shapeFuncsEvaluatedAtCoord );
}


void _FeVariable_PrintLocalOrDomainValues( void* variable, Index localOrDomainCount, Stream* stream ) {
	FeVariable* self = (FeVariable*)variable;
	Node_Index 		node_I = 0;
	Node_GlobalIndex 	gNode_I = 0;
	Dof_Index		currNodeNumDofs;
	Dof_Index		nodeLocalDof_I;
	Dof_EquationNumber	currEqNum;
	Variable*		currVariable;
	
	for( node_I=0; node_I < localOrDomainCount; node_I++ ) {
		gNode_I = self->feMesh->nodeD2G[node_I];
		Journal_Printf( stream, "node %d (global index %d):\n", node_I, gNode_I );
		
		currNodeNumDofs = self->fieldComponentCount;
		
		
		/* Print each dof */
		for ( nodeLocalDof_I = 0; nodeLocalDof_I < currNodeNumDofs; nodeLocalDof_I++ ) {
			currVariable = DofLayout_GetVariable( self->dofLayout, node_I, nodeLocalDof_I );
			Journal_Printf( stream, "\tdof %d \"%s\": %6g - ", nodeLocalDof_I, currVariable->name,
				Variable_GetValueDouble( currVariable, node_I ) );
			currEqNum = self->eqNum->destinationArray[node_I][nodeLocalDof_I];
			if ( currEqNum == -1 ) {
				Journal_Printf( stream, "(from BC)", currEqNum );
			}
			else {
				Journal_Printf( stream, "(eq num %d)", currEqNum );
			}	
			Journal_Printf( stream, "\n", currEqNum );
		}
	}
}

/* TODO: can't assume all swarms have particles of type integrationPoint anymore.
   should check that the given swarm does have I.P for the rest of these functions.*/
double FeVariable_IntegrateElement_AxisIndependent( 
		void* feVariable, void* _swarm, 
		Element_DomainIndex dElement_I, Dimension_Index dim, 
		Axis axis0, Axis axis1, Axis axis2 ) 
{
	FeVariable*          self               = (FeVariable*)         feVariable;
	Swarm*               swarm              = (Swarm*)              _swarm;
	FiniteElement_Mesh*  feMesh             = self->feMesh;
	FiniteElement_Mesh*  geometryMesh       = self->geometryMesh;
	ElementType*         elementType;
	ElementType*         geometryElementType;
	Cell_LocalIndex      cell_I;
	Particle_InCellIndex cParticle_I;
	Particle_InCellIndex cellParticleCount;
	IntegrationPoint*    particle;
	double               detJac;
	double               integral;
	double               value;
	
	/* Initialise Summation of Integral */
	integral = 0.0;

	/* Use feVariable's mesh as geometry mesh if one isn't passed in */
	elementType = FeMesh_ElementTypeAt( feMesh, dElement_I );
	geometryElementType = FeMesh_ElementTypeAt( geometryMesh, dElement_I );

	/* Determine number of particles in element */
	cell_I = CellLayout_MapElementIdToCellId( swarm->cellLayout, dElement_I );
	cellParticleCount = swarm->cellParticleCountTbl[ cell_I ];

	/* Loop over all particles in element */
	for( cParticle_I = 0 ; cParticle_I < cellParticleCount ; cParticle_I++ ) {
		/* Get Pointer to particle */
		particle = (IntegrationPoint*) Swarm_ParticleInCellAt( swarm, cell_I, cParticle_I );

		/* Interpolate Value of Field at Particle */
		FeVariable_InterpolateWithinElement( feVariable, dElement_I, particle->xi, &value );

		/* Calculate Determinant of Jacobian */
		detJac = ElementType_JacobianDeterminant_AxisIndependent( 
				geometryElementType, geometryMesh, dElement_I, particle->xi, dim, axis0, axis1, axis2 );

		/* Sum Integral */
		integral += detJac * particle->weight * value;
	}
	
	return integral;
}

double FeVariable_Integrate( void* feVariable, void* _swarm ) {
	FeVariable*          self               = (FeVariable*)         feVariable;
	Swarm*               swarm              = (Swarm*)              _swarm;
	FiniteElement_Mesh*  feMesh             = self->feMesh;
	Element_LocalIndex   lElement_I;
	Element_LocalIndex   elementLocalCount  = feMesh->elementLocalCount;
	double               integral, integralGlobal;
	
	/* Initialise Summation of Integral */
	integral = 0.0;

	/* Loop over all local elements */
	for ( lElement_I = 0 ; lElement_I < elementLocalCount ; lElement_I++ ) 
		integral += FeVariable_IntegrateElement( self, swarm, lElement_I );
	
	/* Gather and sum integrals from other processors */
	MPI_Allreduce( &integral, &integralGlobal, 1, MPI_DOUBLE, MPI_SUM, self->communicator );

	return integralGlobal;
}

double FeVariable_AverageTopLayer( void* feVariable, void* swarm, Axis layerAxis ) {
	FeVariable*                self               = (FeVariable*)         feVariable;
	FiniteElement_Mesh*        feMesh             = self->feMesh;
	ElementLayout*             elementLayout      = feMesh->layout->elementLayout;
	IJKTopology*               elementTopology    = (IJKTopology*) elementLayout->topology;
	Index                      layerIndex         = elementTopology->size[ layerAxis ] - 1;

	return FeVariable_AverageLayer( self, swarm, layerAxis, layerIndex );
}

double FeVariable_AverageBottomLayer( void* feVariable, void* swarm, Axis layerAxis ) {
	FeVariable*                self               = (FeVariable*)         feVariable;

	return FeVariable_AverageLayer( self, swarm, layerAxis, 0 );
}

double FeVariable_AverageLayer( void* feVariable, void* swarm, Axis layerAxis, Index layerIndex ) {
	FeVariable*                self               = (FeVariable*)         feVariable;
	Axis                       aAxis              = ( layerAxis == I_AXIS ? J_AXIS : I_AXIS );
	Axis                       bAxis              = ( layerAxis == K_AXIS ? J_AXIS : K_AXIS );
	Dimension_Index            dim                = self->dim;
	FiniteElement_Mesh*        feMesh             = self->feMesh;
	BlockGeometry*             geometry           = (BlockGeometry*)feMesh->layout->elementLayout->geometry;
	double                     integral;
	double                     layerThickness     = 0.0;
	
	integral = FeVariable_IntegrateLayer( self, swarm, layerAxis, layerIndex );

	/* Calculate Layer Thickness */
	layerThickness = ((ParallelPipedHexaEL*) feMesh->layout->elementLayout)->elementLengthEachDim[ layerAxis ];

	integral /= layerThickness * ( geometry->max[ aAxis ] - geometry->min[ aAxis ] );
	if ( dim == 3 )
		integral /= geometry->max[ bAxis ] - geometry->min[ bAxis ];

	return integral;
}


double FeVariable_IntegrateLayer_AxisIndependent( 
		void* feVariable, void* _swarm,
		Axis layerAxis, Index layerIndex, Dimension_Index dim, 
		Axis axis0, Axis axis1, Axis axis2 ) 
{ 
	FeVariable*                self               = (FeVariable*)         feVariable;
	Swarm*                     swarm              = (Swarm*)              _swarm;
	FiniteElement_Mesh*        feMesh             = self->feMesh;
	ElementLayout*             elementLayout      = feMesh->layout->elementLayout;
	IJKTopology*               elementTopology    = (IJKTopology*) elementLayout->topology;
	Element_LocalIndex         lElement_I;
	Element_GlobalIndex        gElement_I;
	IJK                        elementIJK;
	Element_LocalIndex         elementLocalCount  = feMesh->elementLocalCount;
	double                     elementIntegral;
	double                     integral;
	double                     integralGlobal;
	Bool                       elementG2LBuiltTemporarily = False;

	Journal_DPrintf( self->debug, "In %s() for FeVariable \"%s\":\n", __func__, self->name );

	/* Modification: build the temporary global tables for speed here. Important for parallel runs. */
	if ( (self->feMesh->buildTemporaryGlobalTables == True) && (self->feMesh->elementG2L == 0) ) {
		self->feMesh->elementG2L = MeshDecomp_BuildElementGlobalToLocalMap( self->feMesh->layout->decomp );
		elementG2LBuiltTemporarily = True;
	}

	/* Initialise Sumation of Integral */
	integral = 0.0;

	Stream_Indent( self->debug );
	for ( gElement_I = 0 ; gElement_I < feMesh->elementGlobalCount ; gElement_I++ ) {
		IJK_1DTo3D( elementTopology, gElement_I, elementIJK );

		/* Check if element is in layer plane */
		if ( elementIJK[ layerAxis ] != layerIndex )
			continue;

		/* Check if element is local */
		lElement_I = Mesh_ElementMapGlobalToLocal( self->feMesh, gElement_I );

		if ( lElement_I >= elementLocalCount ) 
			continue;

		elementIntegral = FeVariable_IntegrateElement_AxisIndependent( self, swarm, lElement_I, dim, axis0, axis1, axis2 );
		Journal_DPrintfL( self->debug, 2, "Integral of element %d was %f\n", lElement_I, elementIntegral );
		integral += elementIntegral;
	}
	Stream_UnIndent( self->debug );

	if ( elementG2LBuiltTemporarily ) {
		Memory_Free( self->feMesh->elementG2L );
		self->feMesh->elementG2L = NULL;
	}


	/* Gather and sum integrals from other processors */
	MPI_Allreduce( &integral, &integralGlobal, 1, MPI_DOUBLE, MPI_SUM, self->communicator );

	Journal_DPrintf( self->debug, "Calculated global integral of layer %d in Axis %d was %f\n", layerIndex, layerAxis, integralGlobal );
	return integralGlobal;
}


double FeVariable_AveragePlane( void* feVariable, Axis planeAxis, double planeHeight ) {
	FeVariable*                self               = (FeVariable*)         feVariable;
	double                     integral;
	Axis                       aAxis              = ( planeAxis == I_AXIS ? J_AXIS : I_AXIS );
	Axis                       bAxis              = ( planeAxis == K_AXIS ? J_AXIS : K_AXIS );
	Dimension_Index            dim                = self->dim;
	BlockGeometry*             geometry           = (BlockGeometry*) self->feMesh->layout->elementLayout->geometry;
	
	integral = FeVariable_IntegratePlane( self, planeAxis, planeHeight );

	integral /= geometry->max[ aAxis ] - geometry->min[ aAxis ];
	if ( dim == 3 )
		integral /= geometry->max[ bAxis ] - geometry->min[ bAxis ];

	return integral;
}

double FeVariable_IntegratePlane( void* feVariable, Axis planeAxis, double planeHeight ) {
	FeVariable*                self               = (FeVariable*)         feVariable;
	FiniteElement_Mesh*        feMesh             = self->feMesh;
	MeshDecomp*                meshDecomp         = feMesh->layout->decomp;
	ElementLayout*             elementLayout      = feMesh->layout->elementLayout;
	IJKTopology*               elementTopology    = (IJKTopology*) elementLayout->topology;
	IJK                        planeIJK;
	Element_LocalIndex         lElement_I;
	Element_GlobalIndex        gElement_I;
	Element_LocalIndex         elementLocalCount  = feMesh->elementLocalCount;
	Axis                       aAxis              = ( planeAxis == I_AXIS ? J_AXIS : I_AXIS );
	Axis                       bAxis              = ( planeAxis == K_AXIS ? J_AXIS : K_AXIS );
	double                     integral;
	/* Swarm Stuff */
	Swarm*                     swarm;
	Bool                       dimExists[]        = { False, False, False };
	ExtensionManager_Register* extensionMgr_Register;
	SingleCellLayout*          singleCellLayout;
	GaussParticleLayout*       gaussParticleLayout;
	Particle_Index             lParticle_I;
	IntegrationPoint*          particle;
	/* Plane location stuff */
	double                     storedXi_J_AXIS;
	Coord                      planeCoord;
	Element_NodeIndex          elementNodeCount;
	Coord**                    globalNodeCoordPtrs;
	ElementType*               elementType;
	double                     planeXi           = -1;
	double                     planeXiGlobal;
	Index                      planeLayer        = 0;
	Index                      planeLayerGlobal;
	Particle_InCellIndex       particlesPerDim[] = {2,2,2};

	/* Find Elements which plane cuts through */
	memcpy( planeCoord, Mesh_CoordAt( feMesh, 0 ), sizeof( Coord ) );
	planeCoord[ planeAxis ] = planeHeight;
	lElement_I = elementLayout->elementWithPoint( elementLayout, meshDecomp, planeCoord, feMesh, 
						      EXCLUSIVE_UPPER_BOUNDARY, 0, NULL );

	if ( lElement_I < elementLocalCount ) {
		Coord planeXiCoord;
		gElement_I = _MeshDecomp_Element_LocalToGlobal1D( meshDecomp, lElement_I );
		IJK_1DTo3D( elementTopology, gElement_I, planeIJK );
		planeLayer = planeIJK[ planeAxis ];
		
		/* Find Local Coordinate of plane */
		elementNodeCount = feMesh->elementNodeCountTbl[lElement_I];
		globalNodeCoordPtrs = Memory_Alloc_Array( Coord*, elementNodeCount, "globalNodeCoordPtrs" );
		Mesh_GetNodeCoordPtrsOfElement( feMesh, lElement_I, globalNodeCoordPtrs );

		elementType = FeMesh_ElementTypeAt( feMesh, lElement_I );
		ElementType_ConvertGlobalCoordToElLocal( elementType, elementLayout,
			(const Coord**) globalNodeCoordPtrs, planeCoord, planeXiCoord );
		planeXi = planeXiCoord[ planeAxis ];
		Memory_Free( globalNodeCoordPtrs );
	}
	
	/* Should be broadcast */
	MPI_Allreduce( &planeXi,    &planeXiGlobal, 1, MPI_DOUBLE, MPI_MAX, self->communicator );
	MPI_Allreduce( &planeLayer, &planeLayerGlobal, 1, MPI_UNSIGNED, MPI_MAX, self->communicator );

	/* Create Swarm in plane */
	extensionMgr_Register = ExtensionManager_Register_New();
	dimExists[ aAxis ] = True;
	if (self->dim == 3)
		dimExists[ bAxis ] = True;
	
	singleCellLayout = SingleCellLayout_New( "cellLayout", dimExists, NULL, NULL );
	particlesPerDim[ planeAxis ] = 1;
	gaussParticleLayout = GaussParticleLayout_New( "particleLayout", self->dim - 1, particlesPerDim );
	swarm = Swarm_New( 
			"gaussSwarm",
			singleCellLayout, 
			gaussParticleLayout,
			self->dim,
			sizeof(IntegrationPoint), 
			extensionMgr_Register, 
			NULL,
			self->communicator );
	Build( swarm, NULL, False );

	/* Change Positions of the particles */
	Initialise( swarm, NULL, False );
	for ( lParticle_I = 0 ; lParticle_I < swarm->particleLocalCount ; lParticle_I++ ) {
		particle = (IntegrationPoint*) Swarm_ParticleAt( swarm, lParticle_I );

		storedXi_J_AXIS = particle->xi[ J_AXIS ];
		particle->xi[ aAxis ]     = particle->xi[ I_AXIS ];
		particle->xi[ bAxis ]     = storedXi_J_AXIS;
		particle->xi[ planeAxis ] = planeXiGlobal;
	}
	
	integral = FeVariable_IntegrateLayer_AxisIndependent( self, swarm, planeAxis, planeLayerGlobal, 
			self->dim - 1, aAxis, bAxis, planeAxis );

	/* Delete */
	Stg_Class_Delete( swarm );
	Stg_Class_Delete( gaussParticleLayout );
	Stg_Class_Delete( singleCellLayout );
	Stg_Class_Delete( extensionMgr_Register );
	
	return integral;
}


void FeVariable_ImportExportInfo_Delete( void* ptr ) {
	/* Nothing to do - the ObjectAdaptor will take care of deleting the actual struct itself */
}

void FeVariable_ImportExportInfo_Print ( void* ptr, struct Stream* stream ) {
	FeVariable_ImportExportInfo*      importExportInfo = (FeVariable_ImportExportInfo*)ptr;

	Journal_Printf( stream, "readNodalValuesFromFile (func ptr): %p", importExportInfo->readNodalValuesFromFile );
	Journal_Printf( stream, "saveNodalValuesToFile (func ptr): %p", importExportInfo->saveNodalValuesToFile );
}

void* FeVariable_ImportExportInfo_Copy( 
		void*					ptr, 
		void*					dest,
		Bool					deep,
		Name					nameExt, 
		struct PtrMap*				ptrMap ) {
	FeVariable_ImportExportInfo*	 feVariable_ImportExportInfo = NULL;		
	/* TODO */
	assert( 0 );
	return feVariable_ImportExportInfo;		
}		


void FeVariable_SaveToFile( void* feVariable, const char* prefixStr, unsigned int timeStep ) {
	FeVariable*                       self = (FeVariable*)feVariable;
	FeVariable_ImportExportInfo*      importExportInfo;
	/* put in an extra loop here that first checks if the feVariable should be
	saved to file. */
	if (self->isCheckpointedAndReloaded) {
	
		importExportInfo = Stg_ObjectList_Get( FeVariable_FileFormatImportExportList, 
								self->exportFormatType );
		assert( importExportInfo );
		importExportInfo->saveNodalValuesToFile( feVariable, prefixStr, timeStep );
	}
		
}


void FeVariable_SaveNodalValuesToFile_StgFEM_Native( void* feVariable, const char* prefixStr, unsigned int timeStep ) {
	FeVariable*        self = (FeVariable*)feVariable;
	char*              filename;
	Node_LocalIndex    lNode_I = 0;
	Node_GlobalIndex   gNode_I = 0;
	double*            coord;
	Dof_Index          dof_I;
	Dof_Index          dofAtEachNodeCount;
	FILE*              outputFile;
	double             variableValues[MAX_FIELD_COMPONENTS];	
	const int          FINISHED_WRITING_TAG = 100;
	int                confirmation = 0;
	int                myRank = self->feMesh->layout->decomp->rank; 
	MPI_Status         status;
	MPI_Comm           comm = self->feMesh->layout->decomp->communicator;
	
	/*                                                prefix             self->name       . 00000 .  dat \0 */
	filename = Memory_Alloc_Array_Unnamed( char, strlen(prefixStr) + strlen(self->name) + 1 + 5 + 1 + 3 + 1 );
	sprintf( filename, "%s%s.%.5u.dat", prefixStr, self->name, timeStep );

	/* wait for go-ahead from process ranked lower than me, to avoid competition writing to file */
	if ( myRank != 0 ) {
		MPI_Recv( &confirmation, 1, MPI_INT, myRank - 1, FINISHED_WRITING_TAG, comm, &status );
	}	

	if ( myRank == 0 ) {
		outputFile = fopen( filename, "w" );
	}
	else {
		outputFile = fopen( filename, "a" );
	}	

	/* Note: assumes same number of dofs at each node */
	dofAtEachNodeCount = self->fieldComponentCount;

	for ( lNode_I = 0; lNode_I < self->feMesh->nodeLocalCount; lNode_I++ ) {
		gNode_I = self->feMesh->nodeL2G[lNode_I];
		fprintf( outputFile, "%u ", gNode_I );
		coord = self->feMesh->nodeCoord[lNode_I];
		fprintf( outputFile, "%.15g %.15g %.15g ", coord[0], coord[1], coord[2] );
		FeVariable_GetValueAtNode( self, lNode_I, variableValues );
		for ( dof_I = 0; dof_I < dofAtEachNodeCount; dof_I++ ) {
			fprintf( outputFile, "%.15g ", variableValues[dof_I] );
		}	
		fprintf( outputFile, "\n" );
	}
	
	fclose( outputFile );
	
	/* send go-ahead from process ranked lower than me, to avoid competition writing to file */
	if ( myRank != (self->feMesh->layout->decomp->nproc-1) ) {
		MPI_Ssend( &confirmation, 1, MPI_INT, myRank + 1, FINISHED_WRITING_TAG, comm );
	}	
	
	Memory_Free( filename );
}


void FeVariable_ReadFromFile( void* feVariable, const char* prefixStr, unsigned int timeStep ) {
	FeVariable*        self = (FeVariable*)feVariable;
	FeVariable_ImportExportInfo*     importExportInfo;
	/* Note: if we wish to support changing grid resolution half way through, then we'll need to make this function smarter,
	and probably save some mesh info as part of the CP so we can work out how to subsample (probably can then just create
	a temporary feVariable using the old size mesh, then keep calling InterpolateValueAt on it whereever our nodes are
	located, though this may not work easily in parallel ),. - PatrickSunter 9 Jun 2006
	*/

	/* Put an extra check in here to check whether the feVariable should be read from a file */
	if (self->isCheckpointedAndReloaded) {
		importExportInfo = Stg_ObjectList_Get( FeVariable_FileFormatImportExportList, 
							self->importFormatType );
		assert( importExportInfo );
		importExportInfo->readNodalValuesFromFile( feVariable, prefixStr, timeStep );
	}
}

#define MAX_LINE_LENGTH_DEFINE 1024

void FeVariable_ReadNodalValuesFromFile_StgFEM_Native( void* feVariable, const char* prefixStr, unsigned int timeStep ) {
	FeVariable*        self = (FeVariable*)feVariable;
	char*              filename;
	Node_LocalIndex    lNode_I = 0;
	Node_GlobalIndex   gNode_I = 0;
	Dof_Index          dof_I;
	Dof_Index          dofAtEachNodeCount;
	FILE*              inputFile;
	double             variableVal;
	char               lineString[MAX_LINE_LENGTH_DEFINE];
	const unsigned int MAX_LINE_LENGTH = MAX_LINE_LENGTH_DEFINE;
	Processor_Index    proc_I=0;
	Dimension_Index    dim_I=0;
	BlockGeometry*     geometry = (BlockGeometry*)self->feMesh->layout->elementLayout->geometry;
	Coord              localGeometryMin;
	Coord              localGeometryMax;
	Bool               nodeG2LBuiltTemporarily = False;
	
	/* Necessary for now since we need to update the geometry min and max - see comment below */
	geometry = Stg_CheckType( geometry, BlockGeometry );

	/*                                                prefix            self->name        . 00000 .  dat \0 */
	filename = Memory_Alloc_Array_Unnamed( char, strlen(prefixStr) + strlen(self->name) + 1 + 5 + 1 + 3 + 1 );
	sprintf( filename, "%s%s.%.5u.dat", prefixStr, self->name, timeStep );

	/* TODO May need/want to change to MPI file stuff */
	
	/* This loop used to stop 2 processors trying to open the file at the same time, which
	  * seems to cause problems */
	for ( proc_I = 0; proc_I < self->feMesh->layout->decomp->nproc; proc_I++ ) {
		MPI_Barrier( self->feMesh->layout->decomp->communicator );

		if ( proc_I == self->feMesh->layout->decomp->rank ) {	
			inputFile = fopen( filename, "r" );
			/* Do the following since in parallel on some systems, the file
			 * doesn't get re-opened at the start automatically. */
			rewind( inputFile );

			if ( False == inputFile ) {
				Stream*    errorStr = Journal_Register( Error_Type, self->type );
				Journal_Printf( errorStr, "Error- in %s(), for feVariable \"%s\": Couldn't find checkpoint file with "
					"prefix \"%s\", timestep %d - thus full filename \"%s\" - aborting.\n", __func__, self->name,
					prefixStr, timeStep, filename );
				exit(EXIT_FAILURE);	
			}
	
		}
	}

	dofAtEachNodeCount = self->fieldComponentCount;

	/* Need to re-set the geometry here, in case we're loading from a checkpoint that had compression/squashing BCs,
		and hence ended up with a smaller mesh than the original */
	for ( dim_I = 0; dim_I < 3; dim_I++ ) {
		localGeometryMin[dim_I] = HUGE_VAL;
		localGeometryMax[dim_I] = -HUGE_VAL;
	}

	/* Modification: build the temporary global tables for speed here. Important for parallel runs. */
	if ( (self->feMesh->buildTemporaryGlobalTables == True) && (self->feMesh->nodeG2L == 0) ) {
		self->feMesh->nodeG2L = MeshDecomp_BuildNodeGlobalToLocalMap( self->feMesh->layout->decomp );
		nodeG2LBuiltTemporarily = True;
	}

	while ( !feof(inputFile) ) {
		fscanf( inputFile, "%u ", &gNode_I );
		lNode_I = Mesh_NodeMapGlobalToLocal( self->feMesh, gNode_I );
		if ( lNode_I != Mesh_Node_Invalid( self->feMesh ) ) {
			/* Note: until we have proper mesh geometry, topology etc checkpointing, we re-load the 
			node co-ords from the feVariable file - and also update the geometry */
			fscanf( inputFile, "%lg %lg %lg ", &self->feMesh->nodeCoord[lNode_I][0], &self->feMesh->nodeCoord[lNode_I][1],
				&self->feMesh->nodeCoord[lNode_I][2] );

			for ( dim_I = 0; dim_I < 3; dim_I++ ) {
				if ( self->feMesh->nodeCoord[lNode_I][dim_I] < localGeometryMin[dim_I] ) {
					localGeometryMin[dim_I] = self->feMesh->nodeCoord[lNode_I][dim_I];
				}
				else if ( self->feMesh->nodeCoord[lNode_I][dim_I] > localGeometryMax[dim_I] ) {
					localGeometryMax[dim_I] = self->feMesh->nodeCoord[lNode_I][dim_I];
				}
			}
			
			for ( dof_I = 0; dof_I < dofAtEachNodeCount; dof_I++ ) {
				fscanf( inputFile, "%lg ", &variableVal );
				DofLayout_SetValueDouble( self->dofLayout, lNode_I, dof_I, variableVal );
			}
		}
		else {
			fgets( lineString, MAX_LINE_LENGTH, inputFile );
		}
	}
	fclose( inputFile );

	if ( nodeG2LBuiltTemporarily ) {
		Memory_Free( self->feMesh->nodeG2L );
		self->feMesh->nodeG2L = NULL;
	}

	/* Note: this is a bit of a hack - we really should be loading the mesh from a separate checkpoint file
	anyway - but for now, we'll just check that its a CornerNL type eg velocity, since BodyNL node points
	are at the centroids, and would hence stuff the Geometry object */
	if ( self->feMesh->layout->nodeLayout->type == CornerNL_Type ) {
		/* TODO: separate into re-usable function */
		/* Since we could be loading in parallel, need to find global min and max of geometry */
		for ( dim_I = 0; dim_I < 3; dim_I++ ) {
			MPI_Allreduce( localGeometryMin, geometry->min, 3, MPI_DOUBLE, MPI_MIN, 
				self->feMesh->layout->decomp->communicator );
			MPI_Allreduce( localGeometryMax, geometry->max, 3, MPI_DOUBLE, MPI_MAX, 
				self->feMesh->layout->decomp->communicator );
		}
	}

	if ( Stg_Class_IsInstance( self->feMesh->layout->elementLayout, ParallelPipedHexaEL_Type ) ) {
		ParallelPipedHexaEL_UpdateGeometryPartitionInfo( self->feMesh->layout->elementLayout, self->feMesh->layout->decomp );
	}
}
