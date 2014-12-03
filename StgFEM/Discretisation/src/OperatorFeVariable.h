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
*/
/** \file
**  Role:
**	The OperatorFeVariable class is used for performing some kind of useful operation on 1 or more FeVariables, to produce another FeVariable. E.G. Taking a gradient, additions, subtraction, etc.
**
** Assumptions:
**
** Comments:
**
**	$Id: OperatorFeVariable.h 843 2007-05-21 22:07:31Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_Discretisation_OperatorFeVariable_h__
#define __StgFEM_Discretisation_OperatorFeVariable_h__

	#define MAX_FIELD_COMPONENTS 9

	/** Textual name of this class */
	extern const Type OperatorFeVariable_Type;

	/** Class contents */
	#define __OperatorFeVariable \
		/* Parent info */ \
		__FeVariable \
		\
		/* Other info */ \
		Name				operatorName; \
		Operator*		_operator; \
		Index				feVariableCount; \
		FeVariable**	feVariableList; \
		Bool				useGradient;

	struct OperatorFeVariable { __OperatorFeVariable };	



	/** Shortcut constructors */
	OperatorFeVariable* OperatorFeVariable_NewUnary(
		Name				name,
		DomainContext*	context,
		void*				_feVariable,
		Name				operatorName );

	OperatorFeVariable* OperatorFeVariable_NewBinary(
		Name				name,
		DomainContext*	context,
		void*				_feVariable1,
		void*				_feVariable2,
		Name				operatorName );
	
	OperatorFeVariable* OperatorFeVariable_NewUnary_OwnOperator( 
		Name				name,
		DomainContext*	context,
		void*				_feVariable,
		Operator*		ownOperator );	
	
	/* Public Constructor */
	void* _OperatorFeVariable_DefaultNew( Name name );

	OperatorFeVariable* OperatorFeVariable_New( 
		Name														name,
		DomainContext*											context,
		void*														feMesh,
      void*														geometryMesh,
      DofLayout*												dofLayout,                                                                                 
      void*														bcs,
      void*														ics,
      void*														linkedDofInfo,
      void*														templateFeVariable,    
		FeVariable_InterpolateWithinElementFunction*	interpolateWithinElement,
		FeVariable_GetValueAtNodeFunction*				getValueAtNode,
		Name														operatorName,
		Operator*												ownOperator,
		Index														feVariableCount,
		FeVariable**											feVariableList,
		Dimension_Index										dim,
		Bool														isCheckpointedAndReloaded,
		MPI_Comm													communicator,
		FieldVariable_Register*								fieldVariable_Register );

	/** Private Constructor */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define OPERATORFEVARIABLE_DEFARGS \
                FEVARIABLE_DEFARGS

	#define OPERATORFEVARIABLE_PASSARGS \
                FEVARIABLE_PASSARGS

	OperatorFeVariable* _OperatorFeVariable_New(  OPERATORFEVARIABLE_DEFARGS  );

	void _OperatorFeVariable_Init( void* feVariable, Name operatorName, Index feVariableCount, FeVariable** feVariableList, Operator* ownOperator );

	/* 'Stg_Class' Virtual Implementations */
	void _OperatorFeVariable_Delete( void* variable );

	void _OperatorFeVariable_Print( void* _feVariable, Stream* stream );

	void* _OperatorFeVariable_Copy( const void* feVariable, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	/* 'Stg_Component' Virtual Implementations */
	void _OperatorFeVariable_AssignFromXML( void* feVariable, Stg_ComponentFactory* cf, void* data );

	void _OperatorFeVariable_Build( void* feVariable, void* data );

	void _OperatorFeVariable_Execute( void* variable, void* data );

	void _OperatorFeVariable_Destroy( void* variable, void* data );

	void _OperatorFeVariable_Initialise( void* variable, void* data );

	/* 'FeVariable Virtual Implementations */
	void _OperatorFeVariable_InterpolateWithinElement( void* feVariable, Element_DomainIndex dElement_I, Coord coord, double* value );

	void _OperatorFeVariable_GetValueAtNode( void* feVariable, Node_DomainIndex dNode_I, double* value);

	InterpolationResult _OperatorFeVariable_InterpolateValueAt( void* variable, Coord globalCoord, double* value );

	void _OperatorFeVariable_SyncShadowValues( void* feVariable );

	/** Private Functions */
	Bool _OperatorFeVariable_CheckIfValidToInterpolateInShadowSpace( OperatorFeVariable* self );

	void _OperatorFeVariable_SetFunctions( void* feVariable );

	void OperatorFeVariable_UnaryInterpolationFunc( void* feVariable, Element_DomainIndex dElement_I, Coord coord, double* value );

	void OperatorFeVariable_BinaryInterpolationFunc( void* feVariable, Element_DomainIndex dElement_I, Coord coord, double* value );

	void OperatorFeVariable_GradientInterpolationFunc( void* feVariable, Element_DomainIndex dElement_I, Coord coord, double* value );

	void OperatorFeVariable_UnaryValueAtNodeFunc( void* feVariable, Node_DomainIndex dNode_I, double* value );

	void OperatorFeVariable_BinaryValueAtNodeFunc( void* feVariable, Node_DomainIndex dNode_I, double* value );

	void OperatorFeVariable_GradientValueAtNodeFunc( void* feVariable, Node_DomainIndex dNode_I, double* value );

#endif /* __StgFEM_Discretisation_OperatorFeVariable_h__ */

