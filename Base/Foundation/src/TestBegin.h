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
*/
/** \file
**  Role:
**
** Assumptions:
**
** Invariants:
**
** Comments:
**
** $Id: TestBegin.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StGermain_Base_Foundation_TestBegin_h__
#define __StGermain_Base_Foundation_TestBegin_h__

#include <mpi.h>

#define TestBegin( name ) Bool test##name { Bool ut_passed = False;
#define TestEnd done:; return passed; }
#define TestFail passed = False; goto done

#define TestTrue( expr )			\
   assert_jmpEnabled = True;			\
   if( !setjmp( assert_env ) ) {		\
      if( !(expr) ) { TestFail; }		\
      assert_jmpEnabled = False;		\
   }						\
   else { TestFail; }

#define TestAssert( stmnt )			\
   assert_jmpEnabled = True;			\
   if( !setjmp( assert_env ) ) {		\
      stmnt;					\
      TestFail;					\
   }						\
   assert_jmpEnabled = False;

#define TestNoAssert( stmnt )			\
   assert_jmpEnabled = True;			\
   if( !setjmp( assert_env )) {			\
      stmnt;					\
      assert_jmpEnabled = False;		\
   }						\
   else { TestFail; }

#endif __StGermain_Base_Foundation_TestBegin_h__
