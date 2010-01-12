/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**   Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**   Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**   Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**   Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**   Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**   Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
** Role:
**   Tests the ElementTypeRegisterSuite
**
** $Id: testElementTypeRegister.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/Discretisation/Discretisation.h>
#include "ElementTypeRegisterSuite.h"

typedef struct {
	ElementType_Register*	etReg;
} ElementTypeRegisterSuiteData;

void ElementTypeRegisterSuite_Setup( ElementTypeRegisterSuiteData* data ) {
	data->etReg = ElementType_Register_New( "elementType_Register" );

	Journal_Enable_AllTypedStream( False );

	_ElementType_Register_Init( data->etReg );

	ElementType_Register_Add( data->etReg, (ElementType*)ConstantElementType_New( "constant" ) );
	ElementType_Register_Add( data->etReg, (ElementType*)BilinearElementType_New( "bilinear" ) );
	ElementType_Register_Add( data->etReg, (ElementType*)TrilinearElementType_New( "trilinear" ) );
	ElementType_Register_Add( data->etReg, (ElementType*)Biquadratic_New( "biquadratic" ) );
	ElementType_Register_Add( data->etReg, (ElementType*)Triquadratic_New( "triquadratic" ) );
}

void ElementTypeRegisterSuite_Teardown( ElementTypeRegisterSuiteData* data ) {
	Stg_Class_Delete( data->etReg );
}

void ElementTypeRegisterSuite_Test( ElementTypeRegisterSuiteData* data ) {
	/* Variables set in this function */
	SizeT                                                                            _sizeOfSelf = sizeof(ConstantElementType);
	Type                                                                                    type = "TestElementType_0";
	Stg_Class_DeleteFunction*                                                            _delete = _ConstantElementType_Delete;
	Stg_Class_PrintFunction*                                                              _print = _ConstantElementType_Print;
	Stg_Class_CopyFunction*                                                                _copy = NULL;
	Stg_Component_DefaultConstructorFunction*                                _defaultConstructor = ConstantElementType_DefaultNew;
	Stg_Component_ConstructFunction*                                                  _construct = _ConstantElementType_AssignFromXML;
	Stg_Component_BuildFunction*                                                          _build = _ConstantElementType_Build;
	Stg_Component_InitialiseFunction*                                                _initialise = _ConstantElementType_Initialise;
	Stg_Component_ExecuteFunction*                                                      _execute = _ConstantElementType_Execute;
	Stg_Component_DestroyFunction*                                                      _destroy = _ConstantElementType_Destroy;
	Name                                                                                    name = "TestElementType_0_Name";
	AllocationType                                                            nameAllocationType = NON_GLOBAL;
	ElementType_EvaluateShapeFunctionsAtFunction*                      _evaluateShapeFunctionsAt = _ConstantElementType_SF_allNodes;
	ElementType_EvaluateShapeFunctionLocalDerivsAtFunction*  _evaluateShapeFunctionLocalDerivsAt = _ConstantElementType_SF_allLocalDerivs_allNodes;
	ElementType_ConvertGlobalCoordToElLocalFunction*                _convertGlobalCoordToElLocal = _ConstantElementType_ConvertGlobalCoordToElLocal;
	ElementType_JacobianDeterminantSurfaceFunction*                  _jacobianDeterminantSurface = _ElementType_JacobianDeterminantSurface;
	ElementType_SurfaceNormalFunction*                                            _surfaceNormal = _ElementType_SurfaceNormal;

	ElementType* elType;
	//unsigned	numTypes	= data->etReg->count;
	unsigned	newIndex;
	unsigned	testIndex;
	
	/* manually create extra types to test the list re-sizing */
	newIndex = ElementType_Register_Add( data->etReg, _ElementType_New(  ELEMENTTYPE_PASSARGS  ) );
	pcu_check_true( newIndex == data->etReg->count - 1 );
	
	type = "TestElementType_1";
	name = "TestElementType_1_Name";

	newIndex = ElementType_Register_Add( data->etReg, _ElementType_New(  ELEMENTTYPE_PASSARGS  ) );
	pcu_check_true( newIndex == data->etReg->count - 1 );

	type = "TestElementType_2";
	name = "TestElementType_2_Name";

	newIndex = ElementType_Register_Add( data->etReg, _ElementType_New(  ELEMENTTYPE_PASSARGS  ) );
	pcu_check_true( newIndex == data->etReg->count - 1 );

	testIndex = ElementType_Register_GetIndex( data->etReg, ConstantElementType_Type );
	elType    = ElementType_Register_At( data->etReg, testIndex );
	pcu_check_true( !strcmp( elType->type, ConstantElementType_Type ) );

	testIndex = ElementType_Register_GetIndex( data->etReg, BilinearElementType_Type );
	elType    = ElementType_Register_At( data->etReg, testIndex );
	pcu_check_true( !strcmp( elType->type, BilinearElementType_Type ) );

	testIndex = ElementType_Register_GetIndex( data->etReg, TrilinearElementType_Type );
	elType    = ElementType_Register_At( data->etReg, testIndex );
	pcu_check_true( !strcmp( elType->type, TrilinearElementType_Type ) );

	testIndex = ElementType_Register_GetIndex( data->etReg, Biquadratic_Type );
	elType    = ElementType_Register_At( data->etReg, testIndex );
	pcu_check_true( !strcmp( elType->type, Biquadratic_Type ) );

	testIndex = ElementType_Register_GetIndex( data->etReg, Triquadratic_Type );
	elType    = ElementType_Register_At( data->etReg, testIndex );
	pcu_check_true( !strcmp( elType->type, Triquadratic_Type ) );
}

void ElementTypeRegisterSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, ElementTypeRegisterSuiteData );
   pcu_suite_setFixtures( suite, ElementTypeRegisterSuite_Setup, ElementTypeRegisterSuite_Teardown );
   pcu_suite_addTest( suite, ElementTypeRegisterSuite_Test );
}


