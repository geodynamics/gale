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
**	Context types.
**
** Assumptions:
**	None as yet.
**
** Comments:
**	None as yet.
**
** $Id: types.h 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StGermain_Base_Context_types_h__
#define __StGermain_Base_Context_types_h__
	
	/* types/classes */
	typedef struct Codelet				Codelet;
	typedef struct _SetVC				SetVC;
	typedef struct _CompositeVC			CompositeVC;
	typedef struct DynamicVC			DynamicVC;
	typedef struct _VariableAllVC_Entry		VariableAllVC_Entry;
	typedef struct _VariableAllVC			VariableAllVC;
	typedef struct _Variable			Variable;
	typedef struct _Variable_Register		Variable_Register;
	typedef struct VariableDumpStream		VariableDumpStream;
	typedef struct _VariableCondition		VariableCondition;
	typedef struct _VariableCondition_Register	VariableCondition_Register;
	typedef struct _ConditionFunction		ConditionFunction;
	typedef struct _ConditionFunction_Register	ConditionFunction_Register;
	typedef Index					VariableAllVC_Entry_Index;
	typedef Index					ConditionFunction_Index;

	/* Variable_Register types */
	typedef Index					Variable_Set_Index;
	typedef Index					Variable_Index;
	typedef Index					Dof_Index;
	
	/* VariableCondition_Register types */
	typedef struct _VariableCondition_Register_Entry VariableCondition_Register_Entry;
	
	/* VariableCondition types */
	typedef enum
	{
		VC_ValueType_Double = 1,
		VC_ValueType_Int,
		VC_ValueType_Short,
		VC_ValueType_Char,
		VC_ValueType_Ptr,
		VC_ValueType_DoubleArray,
		VC_ValueType_CFIndex,
		VC_ValueType_Equation
	} VariableCondition_ValueType;
	
	typedef Index					VariableCondition_Index;
	typedef struct _VariableCondition_Value		VariableCondition_Value;
	typedef struct _VariableCondition_Tuple		VariableCondition_Tuple;
	typedef Index					VariableCondition_ValueIndex;
	typedef Index					VariableCondition_VariableIndex;

	typedef struct _SetVC_Entry			SetVC_Entry;
	typedef Index					SetVC_Entry_Index;
	
	/* CompositeVC types */
	typedef Index					CompositeVC_ItemIndex;

	/* Context types/classes */
	typedef struct AbstractContext		AbstractContext;
	typedef struct ContextEntryPoint	ContextEntryPoint;
	
	typedef AbstractContext			Context;
	
	/* AbstractContext types */
	typedef struct Context_CallInfo	Context_CallInfo;

	typedef Stg_ObjectList			Pointer_Register;
	typedef Stg_ObjectList			Register_Register;
	
#endif /* __StGermain_Base_Context_types_h__ */
