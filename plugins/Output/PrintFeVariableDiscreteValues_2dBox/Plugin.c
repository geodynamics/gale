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
** $Id: Plugin.c 964 2007-10-11 08:03:06Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include "Plugin.h"

#include <stdlib.h>
#include <string.h>

const Type StgFEM_PrintFeVariableDiscreteValues_2dBox_Type = "StgFEM_PrintFeVariableDiscreteValues_2dBox";

const char* PRINT_FE_VARIABLE_DISCRETE_VALUES_2D_BOX_TAG = "PrintFeVariableDiscreteValues_2dBox";

void _StgFEM_PrintFeVariableDiscreteValues_2dBox_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
	FiniteElementContext* context;

	context = (FiniteElementContext*)Stg_ComponentFactory_ConstructByName( cf, (Name)"context", FiniteElementContext, True, data  ); 
	/* Add extensions to nodes, elements and the context */

	/* Add extensions to functionality (entry points) */ 
	ContextEP_Append( context, AbstractContext_EP_Dump, PrintFeVariableDiscreteValues_2dBox );
}

void* _StgFEM_PrintFeVariableDiscreteValues_2dBox_DefaultNew( Name name ) {
	return Codelet_New(
			StgFEM_PrintFeVariableDiscreteValues_2dBox_Type,
			_StgFEM_PrintFeVariableDiscreteValues_2dBox_DefaultNew,
			_StgFEM_PrintFeVariableDiscreteValues_2dBox_AssignFromXML,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
}

Index StgFEM_PrintFeVariableDiscreteValues_2dBox_Register( PluginsManager* pluginsManager ) {
	Journal_DPrintf( StgFEM_Debug, "In: %s( void* )\n", __func__ );

	return PluginsManager_Submit( pluginsManager, StgFEM_PrintFeVariableDiscreteValues_2dBox_Type, (Name)"0", _StgFEM_PrintFeVariableDiscreteValues_2dBox_DefaultNew );
}


void PrintFeVariableDiscreteValues_2dBox( void* _context ) {
	FiniteElementContext*			context = (FiniteElementContext* )_context;
	FeVariable*				currFeVar;
	Stream*					stream;
	Name					currFeVarName;
	Dictionary_Entry_Value*			feVarList=NULL;
	Dictionary_Entry_Value*			currFvParam=NULL;
	Index					feVar_I=0;
	Index					numFeVarsToPrint=0;
	Stream*					warningStr = Journal_Register( Error_Type, (Name)CURR_MODULE_NAME  );
	
	stream = Journal_Register( Info_Type, (Name)CURR_MODULE_NAME  );

	feVarList = Dictionary_Get( context->dictionary, (Dictionary_Entry_Key)(char*)PRINT_FE_VARIABLE_DISCRETE_VALUES_2D_BOX_TAG );
	if ( NULL == feVarList  ) {
		Journal_Printf( warningStr, "Warning - in %s: Plugin \"%s\" loaded, but no \"%s\" tag found "
			"in dictionary. Not printing any FE vars.\n", __func__, CURR_MODULE_NAME, PRINT_FE_VARIABLE_DISCRETE_VALUES_2D_BOX_TAG );
		return;
	}

	numFeVarsToPrint = Dictionary_Entry_Value_GetCount( feVarList );
	if ( 0 == numFeVarsToPrint ) {
		Journal_Printf( warningStr, "Warning - in %s: Plugin \"%s\" loaded, \"%s\" list found, "
			"but list has 0 entries.\n",
			__func__, CURR_MODULE_NAME, PRINT_FE_VARIABLE_DISCRETE_VALUES_2D_BOX_TAG );
		return;
	}


	for ( feVar_I=0; feVar_I < numFeVarsToPrint; feVar_I++ ) {
		currFvParam = Dictionary_Entry_Value_GetElement( feVarList, feVar_I );
		currFeVarName = Dictionary_Entry_Value_AsString( currFvParam );
		currFeVar = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, currFeVarName );
		if ( NULL == currFeVar ) {
			Journal_Printf( warningStr, "Warning - in %s: You requested printing the values of feVariable "
				"\"%s\", but it doesn't exist in the context's field variable register. Skipping.\n",
				__func__, currFeVarName );
			Journal_Printf( warningStr, "(Field Vars currently registered are: " );	
			Stg_ObjectList_PrintAllEntryNames( context->fieldVariable_Register->objects, warningStr );
			Journal_Printf( warningStr, ")\n" );	
			continue;
		}
		Journal_Printf( stream, "%s Values (at end of timestep %d):\n", currFeVarName, context->timeStep );
		FeVariable_PrintLocalDiscreteValues_2dBox( currFeVar, stream );
		Journal_Printf( stream, "\n" );
	}
}


