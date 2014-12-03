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
/** \file
 ** <b>Role:</b>
 **	Abstract class faciliting how class inheritance is done.
 **
 ** <b>Assumptions:</b>
 **	None
 **
 ** <b>Comments:</b>
 **	None
 **
 ** $Id: ISet.h 3904 2006-12-14 00:52:06Z LukeHodkinson $
 **
 **/

#ifndef __StGermain_Base_Foundation_TestBegin_h__
#define __StGermain_Base_Foundation_TestBegin_h__

#include <mpi.h>
#include "debug.h"

#define TestBegin( name ) 				\
   Bool test##name( TestSuite* suite ) {		\
      Bool _passed = True;

#define TestEnd						\
   return _passed;					\
   }

#define TestBarrier					\
   do {							\
      Bool _result;					\
      MPI_Allreduce( &_passed, &_result, 1, MPI_INT, 	\
		     MPI_MIN, MPI_COMM_WORLD );		\
      _passed = _result;				\
      if( !_passed ) goto done;				\
   } while( 0 )

#define TestFail					\
   _passed = False;					\
   TestBarrier

#ifndef NDEBUG

#define TestTrue( expr )				\
   do {							\
      assert_jmpEnabled = True;				\
      if( !setjmp( assert_env ) ) {			\
	 if( !(expr) ) { TestFail; }			\
	 assert_jmpEnabled = False;			\
	 TestBarrier;					\
      }							\
      else { TestFail; }				\
   } while( 0 )

#define TestAssert( stmnt )				\
   do {							\
      assert_jmpEnabled = True;				\
      if( !setjmp( assert_env ) ) {			\
	 stmnt;						\
	 TestFail;					\
      }							\
      else {						\
	 assert_jmpEnabled = False;			\
	 TestBarrier;					\
      }							\
   } while( 0 )

#define TestNoAssert( stmnt )				\
   do {							\
      assert_jmpEnabled = True;				\
      if( !setjmp( assert_env )) {			\
	 stmnt;						\
	 assert_jmpEnabled = False;			\
	 TestBarrier;					\
      }							\
      else { TestFail; }				\
   } while( 0 )

#else

#define TestTrue( expr )				\
   do {							\
      if( !(expr) ) { TestFail; }			\
      TestBarrier;					\
   } while( 0 )

#define TestAssert( stmnt ) stmnt
#define TestNoAssert( stmnt ) stmnt

#endif

#endif /* __StGermain_Base_Foundation_TestBegin_h__ */
