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
** $Id: Init.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdarg.h>
#include "Base/Foundation/Foundation.h"
#include "Base/IO/IO.h"
#include "Base/Container/Container.h"
#include "Base/Automation/Automation.h"
#include "Base/Extensibility/Extensibility.h"

#include "types.h"
#include "Variable.h"
#include "Variable_Register.h"
#include "VariableCondition.h"
#include "VariableCondition_Register.h"
#include "ConditionFunction_Register.h"
#include "ConditionFunction.h"
#include "CompositeVC.h"
#include "DynamicVC.h"
#include "VariableAllVC.h"
#include "SetVC.h"
#include "VariableDumpStream.h"
#include "AbstractContext.h"
#include "ContextEntryPoint.h"
#include "Init.h"

#include <stdio.h>

Bool BaseContext_Init( int* argc, char** argv[] ) {
	Stream* typedStream;
	
	Journal_Printf( Journal_Register( DebugStream_Type, (Name)"Context"  ), "In: %s\n", __func__ ); /* DO NOT CHANGE OR REMOVE */

	condFunc_Register = ConditionFunction_Register_New();
	variableCondition_Register = VariableCondition_Register_New();
	VariableCondition_Register_Add( variableCondition_Register, SetVC_Type, SetVC_Factory );
	VariableCondition_Register_Add( variableCondition_Register, CompositeVC_Type, CompositeVC_Factory );
	VariableCondition_Register_Add( variableCondition_Register, VariableAllVC_Type, VariableAllVC_Factory );

	typedStream = VariableDumpStream_New( VariableDumpStream_Type );
	Stream_Enable( typedStream, False );
	Stream_SetLevel( typedStream, 1 );
	Stream_SetFile( typedStream, stJournal->stdOut );
	
	Journal_RegisterTypedStream( typedStream );
	
	/** Adding default constructors of various components to the Stg_ComponentRegister */
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), Variable_Type, (Name)"0", (Stg_Component_DefaultConstructorFunction*)_Variable_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister( ), CompositeVC_Type, "0", (Stg_Component_DefaultConstructorFunction*)_CompositeVC_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), SetVC_Type, (Name)"0", (Stg_Component_DefaultConstructorFunction*)_SetVC_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister( ), VariableAllVC_Type, "0", (Stg_Component_DefaultConstructorFunction*)_VariableAllVC_DefaultNew );
	Stg_ComponentRegister_Add( Stg_ComponentRegister_Get_ComponentRegister(), DynamicVC_Type, (Name)"0", (Stg_Component_DefaultConstructorFunction*)_DynamicVC_DefaultNew  );

	/** Register Parents for All Classes */
	RegisterParent( Variable_Type, Stg_Component_Type );
	RegisterParent( VariableCondition_Register_Type, Stg_Class_Type );
	RegisterParent( VariableDumpStream_Type, CStream_Type );
	RegisterParent( Variable_Register_Type, Stg_Class_Type );
	RegisterParent( VariableCondition_Type, Stg_Component_Type );
	RegisterParent( ConditionFunction_Type, Stg_Class_Type );
	RegisterParent( ConditionFunction_Register_Type, Stg_Class_Type );
	RegisterParent( CompositeVC_Type, VariableCondition_Type );
	RegisterParent( DynamicVC_Type, VariableCondition_Type );
	RegisterParent( VariableAllVC_Type, VariableCondition_Type );
	RegisterParent( AbstractContext_Type, Stg_Component_Type );
	RegisterParent( ContextEntryPoint_Type, EntryPoint_Type );
	
	return True;
}


