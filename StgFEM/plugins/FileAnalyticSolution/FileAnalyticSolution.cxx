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
** $Id: FileAnalyticSolution.c 989 2007-12-18 13:57:56Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

/* This is taken from Mirko Velic's Analytic Stokes Flow solution */

const Type FileAnalyticSolution_Type = "FileAnalyticSolution";

typedef struct { 
	__AnalyticSolution 
	char* referencePath;
} FileAnalyticSolution;

void FileAnalyticSolution_DummyFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* value ) {
	abort();
}



void _FileAnalyticSolution_AssignFromXML( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	FileAnalyticSolution *   self = (FileAnalyticSolution*)analyticSolution;
	Dictionary_Entry_Value*  varList;
	Index                    var_I;
	char*                    varName;
	FeVariable*              feVarToTest;

	_AnalyticSolution_AssignFromXML( self, cf, data );
	varList = Dictionary_Get( cf->rootDict, (Dictionary_Entry_Key)self->name  );
	Journal_Firewall( varList != NULL, Journal_Register( Error_Type, (Name)self->type  ), 
		"Error- in %s(): Can't find list in XML '%s'\n", __func__, self->name );

	for ( var_I = 0; var_I < Dictionary_Entry_Value_GetCount( varList ); ++var_I ) {
		varName = Dictionary_Entry_Value_AsString( Dictionary_Entry_Value_GetElement( varList, var_I ) );
		feVarToTest = Stg_ComponentFactory_ConstructByName( cf, (Name)varName, FeVariable, True, data  );
	
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, feVarToTest, FileAnalyticSolution_DummyFunction );
	}


	/* Get reference path */
	self->referencePath = Dictionary_GetString_WithDefault( cf->rootDict, "referencePath", "./expected/" );
}

void _FileAnalyticSolution_Initialise( void* analyticSolution, void* data ) {
	FileAnalyticSolution *self = (FileAnalyticSolution*)analyticSolution;
	FeVariable*           analyticFeVariable;
	Index                 analyticFeVariableCount = Stg_ObjectList_Count( self->analyticFeVariableList );
	Index                 analyticFeVariable_I;


	assert( analyticFeVariableCount == Stg_ObjectList_Count( self->analyticFeVariableFuncList ) );

	/* Assign values to all analytic fields */
	for ( analyticFeVariable_I = 0 ; analyticFeVariable_I < analyticFeVariableCount ; analyticFeVariable_I++ ) {
		Stg_Component_Initialise( Stg_ObjectList_At( self->feVariableList, analyticFeVariable_I ), data, False ) ;
		Stg_Component_Initialise( Stg_ObjectList_At( self->analyticFeVariableList, analyticFeVariable_I ), data, False ) ;
		Stg_Component_Initialise( Stg_ObjectList_At( self->errorMagnitudeFieldList, analyticFeVariable_I ), data, False ) ;
		Stg_Component_Initialise( Stg_ObjectList_At( self->relativeErrorMagnitudeFieldList, analyticFeVariable_I ), data, False ) ;

		/* Initialise values from file */
		analyticFeVariable = (FeVariable*) Stg_ObjectList_At( self->analyticFeVariableList, analyticFeVariable_I );
		FeVariable_ReadFromFile( analyticFeVariable, self->referencePath, 0 );
	}
}

void* _FileAnalyticSolution_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(FileAnalyticSolution);
	Type                                                      type = FileAnalyticSolution_Type;
	Stg_Class_DeleteFunction*                              _delete = _AnalyticSolution_Delete;
	Stg_Class_PrintFunction*                                _print = _AnalyticSolution_Print;
	Stg_Class_CopyFunction*                                  _copy = _AnalyticSolution_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _FileAnalyticSolution_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _FileAnalyticSolution_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _AnalyticSolution_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _FileAnalyticSolution_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _AnalyticSolution_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _AnalyticSolution_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _AnalyticSolution_New(  ANALYTICSOLUTION_PASSARGS  );
}

/* This function is automatically run by StGermain when this plugin is loaded. The name must be "<plugin-name>_Register". */
Index StgFEM_FileAnalyticSolution_Register( PluginsManager* pluginsManager ) {
	/* A plugin is only properly registered once it returns the handle provided when submitting a codelet to StGermain. */
	return PluginsManager_Submit( pluginsManager, FileAnalyticSolution_Type, (Name)"0", _FileAnalyticSolution_DefaultNew  );
}


