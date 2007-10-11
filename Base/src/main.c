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
**	Kent Humphries, Software Engineer, VPAC. (kenth@vpac.org)
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
** $Id: main.c 532 2006-04-04 00:21:59Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include "StGermain.h"

#include "main.h"

#include <stdio.h>


void stgMainLoop( Dictionary* dictionary, MPI_Comm CommWorld ) {
	AbstractContext*		context = NULL;
	
	/* Construction phase -----------------------------------------------------------------------------------------------*/
	context = _AbstractContext_New( 
			sizeof(AbstractContext),
	       	        AbstractContext_Type,
	                _AbstractContext_Delete,
	                _AbstractContext_Print,
	                NULL,
	                NULL,
	                _AbstractContext_Construct,
	                _AbstractContext_Build,
	                _AbstractContext_Initialise,
	                _AbstractContext_Execute,
	                _AbstractContext_Destroy,
	                "context",
	                True,
	                NULL,
	                0,
	                0,
	                CommWorld,
	                dictionary );

	/* Construction phase -----------------------------------------------------------------------------------------------*/
	Stg_Component_Construct( context, 0 /* dummy */, &context, True );
	
	/* Building phase ---------------------------------------------------------------------------------------------------*/
	Stg_Component_Build( context, 0 /* dummy */, False );
	
	/* Initialisaton phase ----------------------------------------------------------------------------------------------*/
	Stg_Component_Initialise( context, 0 /* dummy */, False );
	
	/* Run (Solve) phase ------------------------------------------------------------------------------------------------*/
	AbstractContext_Dump( context );
	Stg_Component_Execute( context, 0 /* dummy */, False );

	/* Destruct phase ---------------------------------------------------------------------------------------------------*/
	Stg_Component_Destroy( context, 0 /* dummy */, False );
	Stg_Class_Delete( context );
}

void stgImportToolbox( Dictionary* dictionary, char* toolboxName ) {
	Dictionary_Entry_Value* dev = Dictionary_Get( dictionary, "import" );
	
	if( !dev ) {
		dev = Dictionary_Entry_Value_NewList();
		Dictionary_Add( dictionary, "import", dev );
	}
	
	Dictionary_Entry_Value_AddElement( dev, Dictionary_Entry_Value_FromString( toolboxName ) );
}
