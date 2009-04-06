/*
**  Copyright 2008 Luke Hodkinson
**
**  This file is part of pcu.
**
**  pcu is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  Foobar is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef pcu_runner_h
#define pcu_runner_h

void pcu_runner_init( int argc, char* argv[] );
void pcu_runner_finalise();
void pcu_runner_run( pcu_listener_t* lsnr );
void _pcu_runner_addSuite( const char* name,
			   void (initfunc)( pcu_suite_t* ) );

#define pcu_runner_addSuite( suite, initfunc )  \
   _pcu_runner_addSuite( #suite, initfunc )

#endif
