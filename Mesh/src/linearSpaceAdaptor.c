/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
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
** $Id: LinearSpaceAdaptor.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>
#include <StGermain/StGermain.h>

#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>

#include "types.h"
#include "shortcuts.h"
#include "Grid.h"
#include "Decomp.h"
#include "Sync.h"
#include "MeshTopology.h"
#include "IGraph.h"
#include "MeshClass.h"
#include "MeshGenerator.h"
#include "MeshAdaptor.h"
#include "linearSpaceAdaptor.h"


typedef double (LinearSpaceAdaptor_DeformFunc)( LinearSpaceAdaptor* self, Mesh* mesh, 
					    unsigned* globalSize, unsigned vertex, unsigned* vertexInds );


/* Textual name of this class */
const Type LinearSpaceAdaptor_Type = "LinearSpaceAdaptor";



/* My private functions */
void   LinearSpaceAdaptor_FillTable( segment* table, unsigned size  );
double LinearSpaceAdaptor_MapPoint( segment* table, unsigned size, double x );


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

LinearSpaceAdaptor* LinearSpaceAdaptor_New( Name name, AbstractContext* context ) {
   LinearSpaceAdaptor* self = _LinearSpaceAdaptor_New( sizeof(LinearSpaceAdaptor), 
				    LinearSpaceAdaptor_Type, 
				    _LinearSpaceAdaptor_Delete, 
				    _LinearSpaceAdaptor_Print, 
				    NULL, 
				    (void* (*)(Name))_LinearSpaceAdaptor_New, 
				    _LinearSpaceAdaptor_AssignFromXML, 
				    _LinearSpaceAdaptor_Build, 
				    _LinearSpaceAdaptor_Initialise, 
				    _LinearSpaceAdaptor_Execute, 
				    _LinearSpaceAdaptor_Destroy, 
				    name, 
				    NON_GLOBAL, 
				    _MeshGenerator_SetDimSize, 
				    LinearSpaceAdaptor_Generate );

   _MeshGenerator_Init( (MeshGenerator*)self, context );
   _MeshAdaptor_Init( (MeshAdaptor*)self );
	_LinearSpaceAdaptor_Init( self );

   return self;
}

LinearSpaceAdaptor* _LinearSpaceAdaptor_New( COMPRESSIONADAPTOR_DEFARGS ) {
	LinearSpaceAdaptor* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(LinearSpaceAdaptor) );
	self = (LinearSpaceAdaptor*)_MeshAdaptor_New( MESHADAPTOR_PASSARGS );

	/* Virtual info */
	return self;
}

void _LinearSpaceAdaptor_Init( LinearSpaceAdaptor* self ) {}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _LinearSpaceAdaptor_Delete( void* adaptor ) {
	LinearSpaceAdaptor*	self = (LinearSpaceAdaptor*)adaptor;

	/* Delete the parent. */
	_MeshAdaptor_Delete( self );
}

void _LinearSpaceAdaptor_Print( void* adaptor, Stream* stream ) {
	LinearSpaceAdaptor*	self = (LinearSpaceAdaptor*)adaptor;
	
	/* Set the Journal for printing informations */
	Stream* adaptorStream;
	adaptorStream = Journal_Register( InfoStream_Type, "LinearSpaceAdaptorStream" );

	/* Print parent */
	Journal_Printf( stream, "LinearSpaceAdaptor (ptr): (%p)\n", self );
	_MeshAdaptor_Print( self, stream );
}

void _LinearSpaceAdaptor_AssignFromXML( void* adaptor, Stg_ComponentFactory* cf, void* data ) {
	LinearSpaceAdaptor*	self = (LinearSpaceAdaptor*)adaptor;
	Dictionary_Entry_Value* optionsList;
	Dictionary_Entry_Value* optionSet;
	Index                   segmentCount;
	Index                   segment_I;
	segment*                seg;
	Dictionary*             dictionary  = Dictionary_GetDictionary( cf->componentDict, self->name );
	AbstractContext*        context;

	assert( self );
	assert( cf );

	context = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Context", AbstractContext, True, data );
	if( !context )
		context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data );
	self->loadFromCheckPoint = context->loadFromCheckPoint;

	if( self->loadFromCheckPoint )
	  return;

	/* Call parent construct. */
	_MeshAdaptor_AssignFromXML( self, cf, data );

	self->minX = Dictionary_Entry_Value_AsDouble( Dictionary_Get( cf->rootDict, "minX" ) );
	self->maxX = Dictionary_Entry_Value_AsDouble( Dictionary_Get( cf->rootDict, "maxX" ) );
	self->minY = Dictionary_Entry_Value_AsDouble( Dictionary_Get( cf->rootDict, "minY" ) );
	self->maxY = Dictionary_Entry_Value_AsDouble( Dictionary_Get( cf->rootDict, "maxY" ) );
	self->minZ = Dictionary_Entry_Value_AsDouble( Dictionary_Get( cf->rootDict, "minZ" ) );
	self->maxZ = Dictionary_Entry_Value_AsDouble( Dictionary_Get( cf->rootDict, "maxZ" ) );

	/* Read maping functions - X axis*/
        optionsList = Dictionary_Get( dictionary, "mappingFunctionX" );

	if( optionsList ) {
	  segmentCount = Dictionary_Entry_Value_GetCount(optionsList);
	  self->nSegmentsx = segmentCount;

	  self->tablex = Memory_Alloc_Array( segment , segmentCount, "mapping table x" );
	  memset( self->tablex, 0, segmentCount * sizeof(segment) );

	  for ( segment_I = 0 ; segment_I < segmentCount ; segment_I++) { 
	    optionSet = Dictionary_Entry_Value_GetElement(optionsList, segment_I );
	    seg = &(self->tablex[segment_I]);
	    seg->x = Dictionary_Entry_Value_AsDouble( Dictionary_Entry_Value_GetMember( optionSet, "point" ) );
	    seg->y = Dictionary_Entry_Value_AsDouble( Dictionary_Entry_Value_GetMember( optionSet, "mappedTo" ) );
	  }
	  LinearSpaceAdaptor_FillTable( self->tablex, segmentCount );
	} else {
	  self->nSegmentsx = 0;
	}


	/* Read maping functions - Y axis*/
        optionsList = Dictionary_Get( dictionary, "mappingFunctionY" );
	
	if( optionsList ) {
	  segmentCount = Dictionary_Entry_Value_GetCount(optionsList);
	  self->nSegmentsy = segmentCount;

	  self->tabley = Memory_Alloc_Array( segment , segmentCount, "mapping table y" );
	  memset( self->tabley, 0, segmentCount * sizeof(segment) );

	  for ( segment_I = 0; segment_I < segmentCount; segment_I++) { 
	    optionSet = Dictionary_Entry_Value_GetElement(optionsList, segment_I );
	    seg = &(self->tabley[segment_I]);
	    seg->x = Dictionary_Entry_Value_AsDouble( Dictionary_Entry_Value_GetMember( optionSet, "point" ) );
	    seg->y = Dictionary_Entry_Value_AsDouble( Dictionary_Entry_Value_GetMember( optionSet, "mappedTo" ) );
	  }
	  LinearSpaceAdaptor_FillTable( self->tabley, segmentCount );
	} else {
 	  self->nSegmentsy = 0;
	}


	/* Read maping functions - Z axis*/
        optionsList = Dictionary_Get( dictionary, "mappingFunctionZ" );
	
	if( optionsList ) {
	  segmentCount = Dictionary_Entry_Value_GetCount(optionsList);
	  self->nSegmentsz = segmentCount;

	  self->tablez = Memory_Alloc_Array( segment , segmentCount, "mapping table x" );
	  memset( self->tablez, 0, segmentCount * sizeof(segment) );

	  for ( segment_I = 0 ; segment_I < segmentCount ; segment_I++) { 
	    optionSet = Dictionary_Entry_Value_GetElement(optionsList, segment_I );
	    seg = &(self->tablez[segment_I]);
	    seg->x = Dictionary_Entry_Value_AsDouble( Dictionary_Entry_Value_GetMember( optionSet, "point" ) );
	    seg->y = Dictionary_Entry_Value_AsDouble( Dictionary_Entry_Value_GetMember( optionSet, "mappedTo" ) );
	  }
	  LinearSpaceAdaptor_FillTable( self->tablez, segmentCount );
	} else {
 	  self->nSegmentsz = 0;
	}

	_LinearSpaceAdaptor_Init( self );
}

void _LinearSpaceAdaptor_Build( void* adaptor, void* data ) {
	_MeshAdaptor_Build( adaptor, data );
}

void _LinearSpaceAdaptor_Initialise( void* adaptor, void* data ) {
	_MeshAdaptor_Initialise( adaptor, data );
}

void _LinearSpaceAdaptor_Execute( void* adaptor, void* data ) {
}

void _LinearSpaceAdaptor_Destroy( void* adaptor, void* data ) {
   LinearSpaceAdaptor*		self = (LinearSpaceAdaptor*)adaptor;

   Memory_Free( self->tablex );
   Memory_Free( self->tabley );
   Memory_Free( self->tablez );

   _MeshAdaptor_Destroy( self, data );
}

/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void LinearSpaceAdaptor_Generate( void* adaptor, void* _mesh, void* data ) {

  LinearSpaceAdaptor*		self = (LinearSpaceAdaptor*)adaptor;
  Mesh*				mesh = (Mesh*)_mesh;
  const Sync*			sync;
  double			x;
  Index   			n_i;

  /* Build base mesh, which is assumed to be cartesian. */
  MeshGenerator_Generate( self->generator, mesh, data );

  if( self->loadFromCheckPoint )
    return;

  /* Loop over domain nodes. */
  sync = IGraph_GetDomain( mesh->topo, MT_VERTEX );
  
  if( self->nSegmentsx > 0 ) {
    for( n_i = 0; n_i < Sync_GetNumDomains( sync ); n_i++ ) {
      /* get the x coord */
      x = mesh->verts[n_i][I_AXIS];      

      if( x != self->minX && x != self->maxX ) {	      
	/* normalize the x coord*/
	x = (x - self->minX) / (self->maxX - self->minX);
	/* map the normalized point */
	x = LinearSpaceAdaptor_MapPoint( self->tablex, self->nSegmentsx, x  );
	/* denormalize mapped point */
	x = (self->maxX - self->minX)*x + self->minX;
	/* move the node */
	mesh->verts[n_i][I_AXIS] = x;
      }
    }
  }

  if( self->nSegmentsy > 0 ) {
    for( n_i = 0; n_i < Sync_GetNumDomains( sync ); n_i++ ) {
      /* get the y coord */
      x = mesh->verts[n_i][J_AXIS];      

      if( x != self->minY && x != self->maxY ) {
	/* normalize the y coord*/
	x = (x - self->minY) / (self->maxY - self->minY);
	/* map the normalized point */
	x = LinearSpaceAdaptor_MapPoint( self->tabley, self->nSegmentsy, x  );
	/* denormalize mapped point */
	x = (self->maxY - self->minY)*x + self->minY;
	/* move the node */
	mesh->verts[n_i][J_AXIS] = x;
      }
    }
  }  

  if( self->nSegmentsz > 0 ) {
    for( n_i = 0; n_i < Sync_GetNumDomains( sync ); n_i++ ) {
      /* get the x coord */
      x = mesh->verts[n_i][K_AXIS];

      if( x != self->minZ && x != self->maxZ ) {	      
	/* normalize the z coord*/
	x = (x - self->minZ) / (self->maxZ - self->minZ);
	/* map the normalized point */
	x = LinearSpaceAdaptor_MapPoint( self->tablez, self->nSegmentsz, x  );
	/* denormalize mapped point */
	x = (self->maxZ - self->minZ)*x + self->minZ;
	/* move the node */
	mesh->verts[n_i][K_AXIS] = x;
      }
      /*	    
	FILE* fh;
	fh = fopen( "coords.txt", "a+");
	if( fh != NULL )
	fprintf( fh, "%g\t%g\n", mesh->verts[n_i][I_AXIS], mesh->verts[n_i][J_AXIS] );		
	fclose( fh );*/
    }
  }
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void LinearSpaceAdaptor_FillTable( segment* table, unsigned size  ) {

  Index i;

  table[0].p  = table[0].y / table[0].x;
  table[0].y0 = 0; 
  for( i = 1; i < size; i++ ){
    table[i].p  = (table[i].y - table[i-1].y) / (table[i].x - table[i-1].x);
    table[i].y0 = table[i-1].y;
  }
}

double LinearSpaceAdaptor_MapPoint( segment* table, unsigned size, double x  ) {

  Index i;

  if( x < table[0].x )
    return x * table[0].p + table[0].y0;

  for( i = 1; i < size; i++ ){
    if( x < table[i].x )
      return (x - table[i-1].x) * table[i].p + table[i].y0;
  }
  return 0;
}

