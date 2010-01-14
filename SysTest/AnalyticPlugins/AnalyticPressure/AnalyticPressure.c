/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2005, Monash Cluster Computing 
** All rights reserved.
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
** 		* Redistributions of source code must retain the above copyright notice, 
** 			this list of conditions and the following disclaimer.
** 		* Redistributions in binary form must reproduce the above copyright 
**			notice, this list of conditions and the following disclaimer in the 
**			documentation and/or other materials provided with the distribution.
** 		* Neither the name of the Monash University nor the names of its contributors 
**			may be used to endorse or promote products derived from this software 
**			without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
** THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
** PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS 
** BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
** CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
** SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
** HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
** LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
** OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
**
** Contact:
*%		Louis Moresi - Louis.Moresi@sci.monash.edu.au
*%
** Contributors:
*+		Julian Giordani
*+		Mirko Velic
*+		Robert Turnbull
*+		Vincent Lemiale
*+		Louis Moresi
*+		David May
*+		David Stegman
*+		Patrick Sunter
** $Id: solA.c 636 2006-09-08 03:07:06Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "AnalyticPressure.h"

const Type AnalyticPressure_Type = "PICellerator_AnalyticPressure";

void _PICellerator_AnalyticPressure_PressureFunction( void* _self, double* coord, double* pressure ) {
  AnalyticPressure*  self    = (AnalyticPressure*)_self;
  double             density = self->density;
  double             gravity = self->gravity;
  double             y;

  /* Find coordinate of node */
  y = -coord[ J_AXIS ] + (self->minY+self->maxY)/2.0;

  *pressure = y * density * gravity;
}

void _AnalyticPressure_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data ) {
  AnalyticPressure* 		self = (AnalyticPressure*) _self;

  /* Construct Parent */
  _FieldTest_AssignFromXML( self, cf, data );

    /* Create Analytic Fields */
  self->density  = Stg_ComponentFactory_GetDouble( cf, "layer", (Dictionary_Entry_Key)"density", 0.0  );
  self->gravity  = Stg_ComponentFactory_GetRootDictDouble( cf, (Dictionary_Entry_Key)"gravity", 0  );
  self->maxY     = Stg_ComponentFactory_GetRootDictDouble( cf, (Dictionary_Entry_Key)"maxY", 0  );
  self->minY     = Stg_ComponentFactory_GetRootDictDouble( cf, (Dictionary_Entry_Key)"minY", 0  );
}

void _AnalyticPressure_Build( void* _self, void* data ) {
  AnalyticPressure* 		self = (AnalyticPressure*) _self;

  _FieldTest_Build( self, data );

  /* here we assign the memory and the func ptr for analytic sols */
  self->_analyticSolutionList = Memory_Alloc_Array_Unnamed( FieldTest_AnalyticSolutionFunc*, 1 );
  /* this order MUST be consistent with the xml file definition */
  self->_analyticSolutionList[0] = _PICellerator_AnalyticPressure_PressureFunction;
}

void* _AnalyticPressure_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(AnalyticPressure);
	Type                                                      type = AnalyticPressure_Type;
	Stg_Class_DeleteFunction*                              _delete = _FieldTest_Delete;
	Stg_Class_PrintFunction*                                _print = _FieldTest_Print;
	Stg_Class_CopyFunction*                                  _copy = _FieldTest_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _AnalyticPressure_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _AnalyticPressure_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _AnalyticPressure_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _FieldTest_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _FieldTest_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _FieldTest_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

  return _FieldTest_New(  FIELDTEST_PASSARGS  );
}

Index PICellerator_AnalyticPressure_Register( PluginsManager* pluginsManager ) {
  return PluginsManager_Submit( pluginsManager, AnalyticPressure_Type, (Name)"0", _AnalyticPressure_DefaultNew  );
}


