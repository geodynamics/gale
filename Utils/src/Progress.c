/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd,
** 110 Victoria Street, Melbourne, 3053, Australia.
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
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/Base/Foundation/Foundation.h>
#include <StGermain/Base/IO/IO.h>
#include "types.h"
#include "Progress.h"


void Progress_PrintStatus( void* self );
Bool Progress_CalcStatus( Progress* self );


const Type Progress_Type = "Progress";


Progress* Progress_New() {
	/* Variables set in this function */
	SizeT                      _sizeOfSelf = sizeof(Progress);
	Type                              type = Progress_Type;
	Stg_Class_DeleteFunction*      _delete = _Progress_Delete;
	Stg_Class_PrintFunction*        _print = _Progress_Print;
	Stg_Class_CopyFunction*          _copy = NULL;

   return _Progress_New(  PROGRESS_PASSARGS  );
}


Progress* _Progress_New(  PROGRESS_DEFARGS  )
{
   Progress* self;

   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof(Progress) );
   self = (Progress*)_Stg_Class_New(  STG_CLASS_PASSARGS  );
   _Progress_Init( self );

   return self;
}


void _Progress_Init( void* _self ) {
   Progress* self = (Progress*)_self;

   MPI_Comm_rank( MPI_COMM_WORLD, &self->rank );
   self->title = NULL;
   self->printTitle = True;
   self->preStr = NULL;
   self->width = 40;
   self->start = 0;
   self->end = 0;
   self->pos = 0;
   self->perc = 0;
   self->nBars = 0;
   self->nSpaces = 0;
}


void _Progress_Delete( void* _self ) {
   Progress* self = (Progress*)_self;

   if( self->preStr )
      MemFree( self->preStr );
   if( self->title )
      MemFree( self->title );
   _Stg_Class_Delete( self );
}


void _Progress_Print( void* _self, struct Stream* stream) {
}


void Progress_SetStream( void* _self, Stream* strm ) {
   Progress* self = (Progress*)_self;

   self->strm = strm;
}


void Progress_SetTitle( void* _self, Name str ) {
   Progress* self = (Progress*)_self;

   if( self->title )
      MemFree( self->title );
   self->title = StG_Strdup( str );
}


void Progress_SetPrefix( void* _self, Name str ) {
   Progress* self = (Progress*)_self;

   if( self->preStr )
      MemFree( self->preStr );
   self->preStr = StG_Strdup( str );
}


void Progress_SetRange( void* _self, int start, int end ) {
   Progress* self = (Progress*)_self;

   assert( start <= end );
   self->start = start;
   self->end = end;
   Progress_Restart( self );
}


void Progress_Restart( void* _self ) {
   Progress* self = (Progress*)_self;

   self->printTitle = True;
   self->pos = self->start;
   self->perc = 0;
   self->nBars = 0;
   self->nSpaces = 0;
   Progress_CalcStatus( self );
}


void Progress_Update( void* _self ) {
   Progress* self = (Progress*)_self;

   if( self->rank != 0 || !self->strm )
      return;

   if( self->printTitle && self->title ) {
      Journal_Printf( self->strm, "%s\n", self->title );
      self->printTitle = False;
   }

   Progress_PrintStatus( self );
}


void Progress_Increment( void* _self ) {
   Progress* self = (Progress*)_self;

   self->pos++;
   if( Progress_CalcStatus( self ) )
      Progress_Update( self );
}


void Progress_PrintStatus( void* _self ) {
   Progress* self = (Progress*)_self;
   int ii;

   if( self->rank != 0 || !self->strm )
      return;

   assert( self->width >= 8 );
   if( self->preStr )
      Journal_Printf( self->strm, "%s", self->preStr );
   Journal_Printf( self->strm, "[%3d%%|", self->perc );
   for( ii = 0; ii < self->nBars; ii++ )
      Journal_Printf( self->strm, "=" );
   for( ii = 0; ii < self->nSpaces; ii++ )
      Journal_Printf( self->strm, " " );
   Journal_Printf( self->strm, "]" );
   if( self->perc == 100 )
      Journal_Printf( self->strm, "\n" );
   else
      Journal_Printf( self->strm, "\r" );
}

Bool Progress_CalcStatus( Progress* self ) {
   int oldPerc, oldnBars;
   float frac;

   oldPerc = self->perc;
   oldnBars = self->nBars;

   if( self->start != self->end )
      frac = (float)(self->pos - self->start) / (float)(self->end - self->start);
   else
      frac = 0.0;
   self->perc = (int)(frac * 100);
   self->nBars = (int)(frac * (float)(self->width - 7));
   self->nSpaces = self->width - 7 - self->nBars;

   return (self->perc != oldPerc || self->nBars != oldnBars) ? True : False;
}


