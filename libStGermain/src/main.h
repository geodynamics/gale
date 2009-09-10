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
**	Non as yet.
**
** Comments:
**	None as yet.
**
** $Id: Init.h 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StGermain_main_h__
#define __StGermain_main_h__

	/** The StGermain main - the context life cycle */
	void stgMain( Dictionary* dictionary, MPI_Comm CommWorld );

	/** The StGermain main initialisation */
	Stg_ComponentFactory* stgMainInit( Dictionary* dictionary, MPI_Comm communicator, void* _context );

	/** Initialise the context, from a particular XML file. This saves the user manipulating
	  * an IO_Handler and dictionary to get the data into the context. Useful for test code. */
	Stg_ComponentFactory* stgMainInitFromXML( char* xmlInputFilename, MPI_Comm communicator, void* _context );

	/** The StGermain main loop */
	void stgMainLoop( Stg_ComponentFactory* cf );

	/** The StGermain main destruction */
	void stgMainDestroy( Stg_ComponentFactory* cf );

	/** Add a toolbox to the "import" list in the dictionary */
	void stgImportToolbox( Dictionary* dictionary, char* toolboxName );
	
#endif /* __StGermain_main_h__ */
