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
** $Id: Init.c 4153 2007-07-26 02:25:22Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdarg.h>

#include "Base/Foundation/Foundation.h"
#include "Base/IO/IO.h"
#include "Base/Container/Container.h"

#include "types.h"
#include "shortcuts.h"
#include "Init.h"
#include "Stg_Component.h"
#include "Stg_ComponentRegister.h"
#include "Stg_ComponentFactory.h"
#include "HierarchyTable.h"
#include "CallGraph.h"

#include <stdio.h>

Bool BaseAutomation_Init( int* argc, char** argv[] ) 
{
	Journal_Printf( Journal_Register( DebugStream_Type, "Context" ), "In: %s\n", __func__ ); /* DO NOT CHANGE OR REMOVE */

	/** Initializing the Component Register singleton */
	stgComponentRegister = Stg_ComponentRegister_New( );
	
	/** Initializing the Hierarchy Table singleton */
	stgHierarchyTable = HierarchyTable_New();

	/** Initializing the Call Graph singleton */
	stgCallGraph = Stg_CallGraph_New();

	/** Initializing the ComponentRegister singleton */
	
	/** Register Parents for All Classes */
	RegisterParent( Stg_ComponentFactory_Type,           Stg_Class_Type );
	RegisterParent( Stg_ComponentRegister_Type,          Stg_Class_Type );
	RegisterParent( Stg_Component_Type,                  Stg_Object_Type );
	RegisterParent( HierarchyTable_Type,             HashTable_Type );
	RegisterParent( Stg_CallGraph_Type,              Stg_Class_Type );

	return True;
}


