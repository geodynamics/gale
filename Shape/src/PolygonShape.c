/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org) ) 
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
** $Id: PolygonShape.c 4056 2007-03-29 04:55:51Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/Geometry/Geometry.h>

#include "types.h"
#include "ShapeClass.h"
#include "PolygonShape.h"

#include <assert.h>
#include <string.h>
#include <math.h>

/* Textual name of this class */
const Type PolygonShape_Type = "PolygonShape";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

PolygonShape* PolygonShape_New(
		Name                                  name,
		Dimension_Index                       dim,
		XYZ                                   centre, 
		double                                alpha,
		double                                beta,
		double                                gamma,
		Coord_List                            vertexList,
		Index                                 vertexCount,
		XYZ                                   start,
		XYZ                                   end,
	   Axis                                  perpendicularAxis	)
{ 
	PolygonShape* self = (PolygonShape*)_PolygonShape_DefaultNew( name );

   _Stg_Shape_Init( self, dim, centre, False, alpha, beta, gamma);
   _PolygonShape_Init( self, vertexList, vertexCount, start, end, perpendicularAxis );
	return self;
}

PolygonShape* _PolygonShape_New(  POLYGONSHAPE_DEFARGS  )
{
	PolygonShape* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(PolygonShape) );
	self = (PolygonShape*)_Stg_Shape_New(  STG_SHAPE_PASSARGS  );
	
	/* General info */
	return self;
}

void _PolygonShape_Init( void* polygon, Coord_List vertexList, Index vertexCount, XYZ start, XYZ end, Axis perpendicularAxis ) {
	PolygonShape* self = (PolygonShape*)polygon;
	
	self->vertexList = Memory_Alloc_Array( Coord, vertexCount, "vertexList" );
	memcpy( self->vertexList , vertexList, sizeof(Coord) * vertexCount );
	self->vertexCount = vertexCount;
	memcpy( self->start , start, sizeof(XYZ) );
	memcpy( self->end , end, sizeof(XYZ) );
	self->perpendicularAxis = perpendicularAxis;
}
	
/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _PolygonShape_Delete( void* polygon ) {
	PolygonShape*       self = (PolygonShape*)polygon;
	
	/* Delete parent */
	_Stg_Shape_Delete( self );
}


void _PolygonShape_Print( void* polygon, Stream* stream ) {
	PolygonShape* self = (PolygonShape*)polygon;
	
	/* Print parent */
	_Stg_Shape_Print( self, stream );
}



void* _PolygonShape_Copy( void* polygon, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	PolygonShape*	self = (PolygonShape*)polygon;
	PolygonShape*	newPolygonShape;
	
	newPolygonShape = (PolygonShape*)_Stg_Shape_Copy( self, dest, deep, nameExt, ptrMap );

	newPolygonShape->vertexList = Memory_Alloc_Array( Coord, self->vertexCount, "vertexList" );
	memcpy( newPolygonShape->vertexList , self->vertexList, sizeof(Coord) * self->vertexCount );

	newPolygonShape->vertexList  = self->vertexList;
	newPolygonShape->vertexCount = self->vertexCount;
	memcpy( newPolygonShape->start, self->start, sizeof(XYZ) );
	memcpy( newPolygonShape->end, self->end, sizeof(XYZ) );
	
	return (void*)newPolygonShape;
}

void* _PolygonShape_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                  _sizeOfSelf = sizeof(PolygonShape);
	Type                                                          type = PolygonShape_Type;
	Stg_Class_DeleteFunction*                                  _delete = _PolygonShape_Delete;
	Stg_Class_PrintFunction*                                    _print = _PolygonShape_Print;
	Stg_Class_CopyFunction*                                      _copy = _PolygonShape_Copy;
	Stg_Component_DefaultConstructorFunction*      _defaultConstructor = _PolygonShape_DefaultNew;
	Stg_Component_ConstructFunction*                        _construct = _PolygonShape_AssignFromXML;
	Stg_Component_BuildFunction*                                _build = _PolygonShape_Build;
	Stg_Component_InitialiseFunction*                      _initialise = _PolygonShape_Initialise;
	Stg_Component_ExecuteFunction*                            _execute = _PolygonShape_Execute;
	Stg_Component_DestroyFunction*                            _destroy = _PolygonShape_Destroy;
	Stg_Shape_IsCoordInsideFunction*                    _isCoordInside = _PolygonShape_IsCoordInside;
	Stg_Shape_CalculateVolumeFunction*                _calculateVolume = _PolygonShape_CalculateVolume;
	Stg_Shape_DistanceFromCenterAxisFunction*  _distanceFromCenterAxis = _PolygonShape_DistanceFromCenterAxis;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _PolygonShape_New(  POLYGONSHAPE_PASSARGS  );
}


void _PolygonShape_AssignFromXML( void* polygon, Stg_ComponentFactory* cf, void* data ) {
	PolygonShape*           self       = (PolygonShape*)polygon;
	Index                   vertexCount;
	Index                   vertex_I;
	Coord_List              vertexList;
	XYZ                     start;
	XYZ                     end;
	Axis                    perpendicularAxis;
	char*                   perpendicularAxisName;
	double*                 coord;
	Dictionary_Entry_Value* optionSet;
	Dictionary_Entry_Value* optionsList;
	Dictionary*             dictionary  = Dictionary_GetDictionary( cf->componentDict, self->name );
	Stream*                 stream      = cf->infoStream;
	Stream*                 errorStream = Journal_Register( Error_Type, (Name)self->type  );
	
	_Stg_Shape_AssignFromXML( self, cf, data );

	start[I_AXIS] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"startX", 0.0  );
	end[I_AXIS]   = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"endX", 0.0  );
	start[J_AXIS] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"startY", 0.0  );
	end[J_AXIS]   = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"endY", 0.0  );
	start[K_AXIS] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"startZ", 0.0  );
	end[K_AXIS]   = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"endZ", 0.0  );

	perpendicularAxisName = Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"perpendicularAxis", "z" );
	switch ( perpendicularAxisName[0]  ) {
		case 'x': case 'X': case 'i': case 'I': case '0':
			perpendicularAxis = I_AXIS; break;
		case 'y': case 'Y': case 'j': case 'J': case '1':
			perpendicularAxis = J_AXIS; break;
		case 'z': case 'Z': case 'k': case 'K': case '2':
			perpendicularAxis = K_AXIS; break;
		default:
			Journal_Firewall( False, Journal_Register( Error_Type, (Name)self->type  ),
					"Cannot understand perpendicularAxis '%s'\n", perpendicularAxisName );
	}
	if( self->dim == 3 && ( start[perpendicularAxis] == 0 && end[perpendicularAxis] == 0 ) ) {
		Journal_Firewall( False, Journal_Register( Error_Type, (Name)self->type  ),
		"Problem with %s.\n"
		"You've set the perpendicular axis to be %s, but you've not given the polygon any depth in that axis\n",
	        self->name, perpendicularAxisName );
	}	

	optionsList = Dictionary_Get( dictionary, (Dictionary_Entry_Key)"vertices" );
	
	vertexCount = Dictionary_Entry_Value_GetCount(optionsList );
	Journal_Firewall( vertexCount >= 3, errorStream, 
			"Too few vertices given in trying to build shape '%s' named '%s'.\n"
			"A polygon needs at least three vertices.\n",
			self->type, self->name );

	/* Allocate space */
	vertexList = Memory_Alloc_Array( Coord , vertexCount, "Vertex Array" );
	memset( vertexList, 0, vertexCount * sizeof(Coord) );

	Stream_Indent( stream );
	for ( vertex_I = 0 ; vertex_I < vertexCount ; vertex_I++) { 
		optionSet = Dictionary_Entry_Value_GetElement(optionsList, vertex_I );
		coord = vertexList[vertex_I];
		/* Read Vertex */
		if( perpendicularAxis != I_AXIS )
			coord[ I_AXIS ] = Dictionary_Entry_Value_AsDouble( Dictionary_Entry_Value_GetMember( optionSet, (Dictionary_Entry_Key)"x"));
		if( perpendicularAxis != J_AXIS  )
			coord[ J_AXIS ] = Dictionary_Entry_Value_AsDouble( Dictionary_Entry_Value_GetMember( optionSet, (Dictionary_Entry_Key)"y"));
		if( perpendicularAxis != K_AXIS  )
			coord[ K_AXIS ] = Dictionary_Entry_Value_AsDouble( Dictionary_Entry_Value_GetMember( optionSet, (Dictionary_Entry_Key)"z") );

		/* Print Position */
		Journal_PrintfL( stream, 2, "(%0.3g, %0.3g, %0.3g)\n", coord[I_AXIS], coord[J_AXIS], coord[K_AXIS] );
	}
	Stream_UnIndent( stream );

	_PolygonShape_Init( self, vertexList, vertexCount, start, end, perpendicularAxis );

	Memory_Free( vertexList );
}

void _PolygonShape_Build( void* polygon, void* data ) {
	PolygonShape*	self = (PolygonShape*)polygon;

	_Stg_Shape_Build( self, data );
}
void _PolygonShape_Initialise( void* polygon, void* data ) {
	PolygonShape*	self = (PolygonShape*)polygon;
	
	_Stg_Shape_Initialise( self, data );
}
void _PolygonShape_Execute( void* polygon, void* data ) {
	PolygonShape*	self = (PolygonShape*)polygon;
	
	_Stg_Shape_Execute( self, data );
}
void _PolygonShape_Destroy( void* polygon, void* data ) {
	PolygonShape*	self = (PolygonShape*)polygon;

	Coord_List     vertexList = self->vertexList;
	Memory_Free( vertexList );

	_Stg_Shape_Destroy( self, data );
}

/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/
	
/*---------------------------------------------------------------------------------------------------------------------
** Private Member functions
*/
/* Algorithm describe Paul Bourke's page http://astronomy.swin.edu.au/~pbourke/geometry/insidepoly/ (see solution 2)
 *
 * Algorithm works by summing the angles between the test coordinate and each pair of vertices that make up an edge 
 * in the polygon. An inside point will give an angle of 2pi and and outside point will give an angle of 0 */

Bool _PolygonShape_IsCoordInside( void* polygon, Coord coord ) {
	PolygonShape*        self                = (PolygonShape*) polygon;
	Index           vertexCount         = self->vertexCount;
	Coord_List      vertexList          = self->vertexList;
	Axis            perpendicularAxis   = self->perpendicularAxis;
	XYZ             vectorToStartVertex = { 0.0, 0.0, 0.0 };
	XYZ             vectorToEndVertex   = { 0.0, 0.0, 0.0 };
	XYZ             crossproduct        = { 0.0, 0.0, 0.0 };
	double          currAngle;
	double          totalAngle          = 0.0;
	Index           vertex_I;
	double*         startVertex;
	double*         endVertex;
       Coord           newCoord;

	/* Transform coordinate into canonical reference frame */
	Stg_Shape_TransformCoord( self, coord, newCoord );

	/* Check to make sure that the coordinate is within startZ and endZ in 3D */
	if ( self->dim == 3 && ( newCoord[ perpendicularAxis ] < self->start[perpendicularAxis] || newCoord[ perpendicularAxis ] > self->end[perpendicularAxis] ))
		return False;	

	for ( vertex_I = 0 ; vertex_I < vertexCount ; vertex_I++ ) {
		/* Get vertices of current edge */
		startVertex = vertexList[ vertex_I ];
		endVertex   = vertexList[ (vertex_I + 1) % vertexCount ];

		/* Work out vectors */
		StGermain_VectorSubtraction( vectorToStartVertex, newCoord, startVertex, 3 );
		StGermain_VectorSubtraction( vectorToEndVertex,   newCoord, endVertex,  3 );

		vectorToStartVertex[ perpendicularAxis ] = 0;
		vectorToEndVertex[ perpendicularAxis ] = 0;

		/* Work out angle - just by doing dot product - will always be positive */
		currAngle = StGermain_AngleBetweenVectors( vectorToStartVertex, vectorToEndVertex, 3 );

		/* Work out 'sign' of angle but working out cross product */
		StGermain_VectorCrossProduct( crossproduct, vectorToEndVertex, vectorToStartVertex );

		if ( crossproduct[ perpendicularAxis ] > 0.0 )
			totalAngle += currAngle;
		else
			totalAngle -= currAngle;
	}


	/* work out whether the coord is within the polygon */
	if ( fabs( totalAngle ) < M_PI )
		return False;
	else  
		return True;
}


double _PolygonShape_CalculateVolume( void* polygon ) {
	assert( 0 );
	return 0.0;
}

void _PolygonShape_DistanceFromCenterAxis( void* shape, Coord coord, double* disVec ){
	Stg_Shape* self = (Stg_Shape*)shape;
	Journal_Firewall( False, Journal_Register( Error_Type, (Name)self->type  ),
	"Error in function %s: This functions hasn't been implemented.", 
	"Please inform underworld-dev@vpac.org you've received this error.\n", __func__ );
}


