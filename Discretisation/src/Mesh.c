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
** $Id: Mesh.c 656 2006-10-18 06:45:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include "units.h"
#include "types.h"
#include "shortcuts.h"
#include "Mesh.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "Element.h"
#include "ElementType.h"
#include "ConstantElementType.h"
#include "BilinearElementType.h"
#include "TrilinearElementType.h"
#include "LinearTriangleElementType.h"
#include "ElementType_Register.h"

const Type FiniteElement_Mesh_Type = "FiniteElement_Mesh";

void* FiniteElement_Mesh_DefaultNew( Name name ) {
	return _FiniteElement_Mesh_New(
		sizeof(FiniteElement_Mesh),
		FiniteElement_Mesh_Type,
		_FiniteElement_Mesh_Delete, 
		_FiniteElement_Mesh_Print,
		_FiniteElement_Mesh_Copy,
		FiniteElement_Mesh_DefaultNew,
		_FiniteElement_Mesh_Construct,
		_FiniteElement_Mesh_Build,
		_FiniteElement_Mesh_Initialise,
		_FiniteElement_Mesh_Execute, 
		_FiniteElement_Mesh_Destroy,
		name,
		False,
		_Mesh_Node_IsLocal1D,
		_Mesh_Node_IsShadow1D,
		_Mesh_Element_IsLocal1D,
		_Mesh_Element_IsShadow1D,
		NULL, 
		0,
		0, 
		NULL,
		NULL,
		NULL);
}

FiniteElement_Mesh* FiniteElement_Mesh_New(
		Name					name,
		void*					layout,
		SizeT					nodeSize,
		SizeT					elementSize,
		void*					extensionMgr_Register,
		void*					elementType_Register,
		Dictionary*				dictionary )
{
	return _FiniteElement_Mesh_New(
		sizeof(FiniteElement_Mesh),
		FiniteElement_Mesh_Type,
		_FiniteElement_Mesh_Delete, 
		_FiniteElement_Mesh_Print,
		_FiniteElement_Mesh_Copy,
		FiniteElement_Mesh_DefaultNew,
		_FiniteElement_Mesh_Construct,
		_FiniteElement_Mesh_Build,
		_FiniteElement_Mesh_Initialise,
		_FiniteElement_Mesh_Execute, 
		_FiniteElement_Mesh_Destroy, 
		name,
		True,
		_Mesh_Node_IsLocal1D,
		_Mesh_Node_IsShadow1D,
		_Mesh_Element_IsLocal1D,
		_Mesh_Element_IsShadow1D,
		layout, 
		nodeSize,
		elementSize, 
		extensionMgr_Register,
		elementType_Register,
		dictionary );
}

void FiniteElement_Mesh_Init(
		FiniteElement_Mesh*			self,
		Name					name,
		void*					layout,
		SizeT					nodeSize,
		SizeT					elementSize,
		void*					extensionMgr_Register,
		void*					elementType_Register,
		Dictionary*				dictionary )
{
	/* General info */
	self->type = FiniteElement_Mesh_Type;
	self->_sizeOfSelf = sizeof(FiniteElement_Mesh);
	self->_deleteSelf = False;
	
	/* Virtual info */
	self->_delete = _FiniteElement_Mesh_Delete;
	self->_print = _FiniteElement_Mesh_Print;
	self->_copy = _FiniteElement_Mesh_Copy;
	self->_defaultConstructor = FiniteElement_Mesh_DefaultNew;
	self->_construct = _FiniteElement_Mesh_Construct;
	self->_build = _FiniteElement_Mesh_Build;
	self->_initialise = _FiniteElement_Mesh_Initialise;
	self->_execute = _FiniteElement_Mesh_Execute;
	self->_destroy = _FiniteElement_Mesh_Destroy;
	self->nodeIsLocal = _Mesh_Node_IsLocal1D;
	self->nodeIsShadow = _Mesh_Node_IsShadow1D;
	self->elementIsLocal = _Mesh_Element_IsLocal1D;
	self->elementIsShadow = _Mesh_Element_IsShadow1D;

	_Stg_Class_Init( (Stg_Class*)self );
	_Stg_Object_Init( (Stg_Object*)self, name, NON_GLOBAL );
	_Stg_Component_Init( (Stg_Component*)self );
	_Mesh_Init( (Mesh*)self, layout, nodeSize, elementSize, extensionMgr_Register );
	
	/* FiniteElement_Mesh info */
	_FiniteElement_Mesh_Init( self, elementType_Register );
}


FiniteElement_Mesh* _FiniteElement_Mesh_New( 
		SizeT					_sizeOfSelf,
		Type					type,
		Stg_Class_DeleteFunction*			_delete,
		Stg_Class_PrintFunction*			_print,
		Stg_Class_CopyFunction*			_copy, 
		Stg_Component_DefaultConstructorFunction*	_defaultConstructor,
		Stg_Component_ConstructFunction*	_construct,
		Stg_Component_BuildFunction*		_build,
		Stg_Component_InitialiseFunction*		_initialise,
		Stg_Component_ExecuteFunction*		_execute,
		Stg_Component_DestroyFunction*		_destroy,
		Name					name,
		Bool					initFlag,
		Mesh_Node_IsLocalFunction*		_nodeIsLocal,
		Mesh_Node_IsShadowFunction*		_nodeIsShadow,
		Mesh_Element_IsLocalFunction*		_elementIsLocal,
		Mesh_Element_IsShadowFunction*		_elementIsShadow,
		void*					layout,
		SizeT					nodeSize,
		SizeT					elementSize, 
		void*					extensionMgr_Register,
		void*					elementType_Register,
		Dictionary*				dictionary )
{
	FiniteElement_Mesh*			self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(FiniteElement_Mesh) );
	/* Construct using parent */
	self = (FiniteElement_Mesh*)_Mesh_New( _sizeOfSelf, type, _delete, _print, _copy, _defaultConstructor, _construct,
			_build, _initialise, _execute, _destroy, name, initFlag, _nodeIsLocal, _nodeIsShadow, _elementIsLocal, _elementIsShadow,
		layout, nodeSize, elementSize, extensionMgr_Register, dictionary );

	/* General info */
	
	/* Virtual functions */
	
	/* FiniteElement_Mesh info */
	if( initFlag ){
		_FiniteElement_Mesh_Init( self, elementType_Register );
	}
	
	return self;
}

void _FiniteElement_Mesh_Init(
		FiniteElement_Mesh*			self,
		void*					elementType_Register )
{
	/* General and Virtual info should already be set */
	
	/* FiniteElement_Mesh info */
	self->isConstructed = True;
	self->elementType_Register = (ElementType_Register*)elementType_Register;

	self->debug = Stream_RegisterChild( StgFEM_Discretisation_Debug, self->type );

	/* Check the mesh is appropriately configured for F.E. problems */
	Mesh_ActivateElementNodeTbl( self );
	Mesh_ActivateNodeElementTbl( self );
	Mesh_ActivateNodeLocalToGlobalMap( self );
	if ( self->layout->decomp->allowPartitionOnElement == True ) {
		Stream* feMeshError = Journal_Register( ErrorStream_Type, self->type );
		if ( self->isBuilt ) {
			Journal_Printf( feMeshError, "Error: Mesh already built using illegal partition on element." );
			assert( 0 );
		}
		else {
			Journal_Printf( feMeshError, "Warning: setting Mesh to only partition on node." );
			self->layout->decomp->allowPartitionOnElement = False;
			self->layout->decomp->allowPartitionOnNode = True;
		}
	}
}

void _FiniteElement_Mesh_Delete( void* feMesh ) {
	FiniteElement_Mesh* self = (FiniteElement_Mesh*)feMesh;
	
	Journal_DPrintf( self->debug, "In %s - for %s\n", __func__, self->name );
	/* Stg_Class_Delete parent*/
	_Mesh_Delete( self );
	
	/* mesh, elementType_Register are purposely not deleted */
}

void _FiniteElement_Mesh_Print( void* feMesh, Stream* stream ) {
	FiniteElement_Mesh* self = (FiniteElement_Mesh*)feMesh;
	
	/* General info */
	Journal_Printf( stream, "FiniteElement_Mesh (ptr): %p\n", self );
	
	/* Print parent */
	_Mesh_Print( self, stream );
	
	/* Virtual info */
	
	/* FiniteElement_Mesh info */
	Journal_Printf( stream, "Available element types:\n" );
	Print( self->elementType_Register, stream );
}


void* _FiniteElement_Mesh_Copy( void* feMesh, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	FiniteElement_Mesh*	self = (FiniteElement_Mesh*)feMesh;
	FiniteElement_Mesh*	newFeMesh;
	PtrMap*			map = ptrMap;
	Bool			ownMap = False;
	
	if( !map ) {
		map = PtrMap_New( 10 );
		ownMap = True;
	}
	
	newFeMesh = (FiniteElement_Mesh*)_Mesh_Copy( self, dest, deep, nameExt, map );
	
	newFeMesh->elementType_Register = self->elementType_Register;
	
	newFeMesh->debug = self->debug;
	
	if( ownMap ) {
		Stg_Class_Delete( map );
	}
	
	return (void*)newFeMesh;
}

void _FiniteElement_Mesh_Construct( void* feMesh, Stg_ComponentFactory* cf, void* data ) 
{
	FiniteElement_Mesh *self = (FiniteElement_Mesh*)feMesh;
	void *elementType_Register = NULL;

	_Mesh_Construct( self, cf, data );
	
	elementType_Register = Stg_ObjectList_Get( cf->registerRegister, "ElementType_Register" );
	assert( elementType_Register );
	_FiniteElement_Mesh_Init( (FiniteElement_Mesh*)self, elementType_Register );
}

void _FiniteElement_Mesh_Build( void* feMesh, void* data ) {
	FiniteElement_Mesh* self = (FiniteElement_Mesh*)feMesh;
	
	Journal_DPrintf( self->debug, "In %s - for %s\n", __func__, self->name );
	_Mesh_Build( self, data );

}


void _FiniteElement_Mesh_Initialise( void* feMesh, void* data ) {
	FiniteElement_Mesh* self = (FiniteElement_Mesh*)feMesh;
	
	Journal_DPrintf( self->debug, "In %s - for %s\n", __func__, self->name );
	_Mesh_Initialise( self, data );
	Stream_IndentBranch( StgFEM_Debug );

	/* TODO: generalise into an entry point so it can easily be over-ridden */
	_FiniteElement_Mesh_SetupElementTypes_Default( self );
	
	#if DEBUG
	if ( Stream_IsPrintableLevel( self->debug, 3 ) ) {
		_Mesh_PrintCoords( (Mesh*) self, self->debug );
	}
	#endif
	Stream_UnIndentBranch( StgFEM_Debug );
}


void _FiniteElement_Mesh_SetupElementTypes_Default( FiniteElement_Mesh* self ) {
	ElementLayout*     elementLayout         = self->layout->elementLayout;
	Type               elementLayoutType     = elementLayout->type;
	Type               nodeLayoutType        = self->layout->nodeLayout->type;
	ElementType_Index  elementTypeIndexToUse = (ElementType_Index) -1;
	Type               elementTypeToUse      = NULL;
	Element_LocalIndex element_I             = 0;
	Stream*            error                 = Journal_Register( ErrorStream_Type, self->type );
	
	Journal_DPrintf( self->debug, "In %s - for %s\n", __func__, self->name );
	Stream_IndentBranch( StgFEM_Debug );

	/* TODO: this will need to be generalised when we allow meshes that are composed of different elements */

	if ( nodeLayoutType == BodyNL_Type ) {
		/* A Body (1 node at centre of each element) NodeLayout implies Constant shape func FE element types */
		elementTypeIndexToUse = ElementType_Register_GetIndex( self->elementType_Register, ConstantElementType_Type );
		Journal_DPrintf( self->debug, "Detected %s node layout -> setting all my element types to %s.\n",
			BodyNL_Type, ConstantElementType_Type );
	}
	else {
		if ( strstr( elementLayoutType, IrregEL_Type ) ) {
			/* Irregular topology: only have linear triangular so far */
			if ( nodeLayoutType == CornerNL_Type ) {
				elementTypeToUse = LinearTriangleElementType_Type;
				elementTypeIndexToUse = ElementType_Register_GetIndex( self->elementType_Register, LinearTriangleElementType_Type );
			}	
		}		
		else if ( strstr( elementLayoutType, HexaEL_Type ) ) {
			/* Hexahedral 3D elements: */
			if ( ((HexaEL*)elementLayout)->dim == 2 ) {
				if ( nodeLayoutType == CornerNL_Type ) {
					elementTypeToUse = BilinearElementType_Type;
					elementTypeIndexToUse = ElementType_Register_GetIndex( self->elementType_Register, BilinearElementType_Type );
				}
			}
			else {
				elementTypeToUse = TrilinearElementType_Type;
				if ( nodeLayoutType == CornerNL_Type ) {
					elementTypeIndexToUse = ElementType_Register_GetIndex( self->elementType_Register, TrilinearElementType_Type );
				}
			}
		}
		else {
			Journal_Firewall( 0, error, "Error: Don't know what FE Element type to set for element "
				"layout/node layout combo %s/%s. Please override this Hook with code to set the FE "
				"element types for your non-standard mesh.\n", elementLayoutType, nodeLayoutType );
		}

		Journal_DPrintf( self->debug, "Detected %s el. layout, %s node layout -> setting all element types to %s "
			"(index %u).\n", elementLayoutType, nodeLayoutType, elementTypeToUse, elementTypeIndexToUse );
	}		

	Journal_Firewall( (elementTypeIndexToUse != (unsigned int)-1), error,
		"Stuffup: should have worked out which element type to use by now.\n" );

	for( element_I = 0; element_I < self->elementDomainCount; element_I++ ) {
		FeMesh_ElementAt( self, element_I )->elementType_I = elementTypeIndexToUse;
		/* TODO: stick to default FE cell index till we sort that design issue out... */
		FeMesh_ElementAt( self, element_I )->cell_I = 0;
	}
	Stream_UnIndentBranch( StgFEM_Debug );

}


void _FiniteElement_Mesh_Execute( void* feMesh, void* data ) {
	/* Do nothing */
}

void _FiniteElement_Mesh_Destroy( void* feMesh, void* data ) {
	/* Do nothing */
}

double* _FiniteElement_Mesh_CoordAt( void* feMesh, Node_DomainIndex node_dI ) {
	FiniteElement_Mesh* self = (FiniteElement_Mesh*)feMesh;
	
	return FiniteElement_Mesh_CoordAt( self, node_dI );
}

Node* _FiniteElement_Mesh_NodeAt( void* feMesh, Node_DomainIndex node_dI ) {
	FiniteElement_Mesh* self = (FiniteElement_Mesh*)feMesh;
	
	return FiniteElement_Mesh_NodeAt( self, node_dI );
}

FiniteElement_Element* _FiniteElement_Mesh_ElementAt( void* feMesh, Element_DomainIndex element_I ) {
	FiniteElement_Mesh* self = (FiniteElement_Mesh*)feMesh;
	
	return FiniteElement_Mesh_ElementAt( self, element_I );
}

ElementType* _FiniteElement_Mesh_ElementTypeAt( void* feMesh, Element_DomainIndex element_I ) {
	FiniteElement_Mesh* self = (FiniteElement_Mesh*)feMesh;
	
	return FiniteElement_Mesh_ElementTypeAt( self, element_I );
}

Variable* FiniteElement_Mesh_RegisterNodeCoordsAsVariables( void* feMesh, void* _variable_Register, Variable** variableList ) {
	FiniteElement_Mesh* self              = (FiniteElement_Mesh*)feMesh;
	Variable_Register*  variable_Register = (Variable_Register*) _variable_Register;
	Variable*           variable;
	Name                variableName;
	Name                variableNameX;
	Name                variableNameY;
	Name                variableNameZ;
	
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
		3, 
		&self->nodeDomainCount, 
		(void**)&self->nodeCoord, 
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

Variable* FiniteElement_Mesh_RegisterLocalNodeCoordsAsVariables( void* feMesh, void* _variable_Register, Variable** variableList ) {
	FiniteElement_Mesh* self              = (FiniteElement_Mesh*)feMesh;
	Variable_Register*  variable_Register = (Variable_Register*) _variable_Register;
	Variable*           variable;
	Name                variableName;
	Name                variableNameX;
	Name                variableNameY;
	Name                variableNameZ;
	
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
		3, 
		&self->nodeLocalCount, 
		(void**)&self->nodeCoord, 
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

void FiniteElement_Mesh_CalcGlobalCoordFromLocalCoord( void* feMesh, Dimension_Index dim, Element_LocalIndex lElement_I, Coord xi, Coord globalCoord ) {
	FiniteElement_Mesh*        self                = (FiniteElement_Mesh*) feMesh;
	ElementType*               elementType         = FeMesh_ElementTypeAt( self, lElement_I );
	Node_Index                 elementNodeCount    = elementType->nodeCount;
	Node_Index                 elLocalNode_I;
	Node_Index                 lNode_I;
	double*                    evaluatedShapeFuncs;
	double                     shapeFunc;
	double*                    nodeCoord;
	
	/* Evaluate Shape Functions at this point */
	evaluatedShapeFuncs = Memory_Alloc_Array( double, elementNodeCount, "evaluatedShapeFuncs" );
	ElementType_EvaluateShapeFunctionsAt( elementType, xi, evaluatedShapeFuncs );
	
	/* See FEM/BEM notes pg. 9
	x(Xi) = \Sigma phi_n(xi) * x_n  */
	memset( globalCoord, 0, sizeof(Coord) );
	for ( elLocalNode_I = 0; elLocalNode_I < elementNodeCount; elLocalNode_I++ ) {
		shapeFunc = evaluatedShapeFuncs[elLocalNode_I];
		lNode_I   = Mesh_Element_Node_I( self, lElement_I, elLocalNode_I );
		nodeCoord = Mesh_CoordAt( self, lNode_I );

		globalCoord[ 0 ] += shapeFunc * nodeCoord[ 0 ];
		globalCoord[ 1 ] += shapeFunc * nodeCoord[ 1 ];
		if ( dim == 3 )
			globalCoord[ 2 ] += shapeFunc * nodeCoord[ 2 ];
	}	

	Memory_Free( evaluatedShapeFuncs );
}


void FiniteElement_Mesh_CalcLocalCoordFromGlobalCoord(
		void* feMesh,
		Element_DomainIndex elementCoordIn,
		const Coord globalCoord,
		Coord elLocalCoord )
{
	FiniteElement_Mesh*     self                = (FiniteElement_Mesh*) feMesh;
	Node_LocalIndex		currElementNodeCount=0;
	Coord**			globalNodeCoordPtrs=NULL;
	ElementType*		elementType = NULL;

	// TODO: would be good to do an optional cautious check that the provided global 
	// coord is actually in the provided element 

	/* convert global coordinate to local co-ordinates of element the coord is in */
	currElementNodeCount = self->elementNodeCountTbl[elementCoordIn];
	globalNodeCoordPtrs = Memory_Alloc_Array( Coord*, currElementNodeCount, "globalNodeCoordPtrs" );
	Mesh_GetNodeCoordPtrsOfElement( self, elementCoordIn, globalNodeCoordPtrs );

	elementType = FeMesh_ElementTypeAt( self, elementCoordIn );
	ElementType_ConvertGlobalCoordToElLocal( elementType, self->layout->elementLayout,
		(const Coord**) globalNodeCoordPtrs, globalCoord, elLocalCoord );

	Memory_Free( globalNodeCoordPtrs );
}
