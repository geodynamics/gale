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
**	Encapsulates the submission information to StGermain and running of "operator" functions (including mapping of operator 
**		parameters to the order in the data array to expect them).
**
** Assumptions:
**
** Comments:
**
** $Id: $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_Assembly_OperatorFunction_h__
#define __StgFEM_Assembly_OperatorFunction_h__

	typedef union {
		void*				_ptr;
		double				_double;
		int				_int;
		unsigned			_unsigned;
	} OperatorFunction_Datum;

	typedef void (OperatorFunction_ApplyFunc) ( 
		OperatorFunction*		self, /* so the integration swarm can be accessed */
		StiffnessMatrix*		sm, 
		Element_LocalIndex		lElement_I, 
		double**			elsm, 
		OperatorFunction_Datum*		data );
	
	typedef void (OperatorFunction_ApplyRHSFunc) ( 
		OperatorFunction*		self, /* so the integration swarm can be accessed */
		StiffnessMatrix*		sm, 
		Element_LocalIndex		lElement_I, 
		double*				elvec, 
		OperatorFunction_Datum*		data );
	

	extern const Type OperatorFunction_Type;
	
	#define __OperatorFunction \
		/* General info */ \
		__Stg_Component \
		\
		/* Virtual info */ \
		\
		/* As a OperatorFunction... */ \
		/*OperatorFunction_ApplyFunc*	apply;*/ \
		/*Name				name; */\
		Index				_dataCount; \
		Name*				_dataName; \
		Bool*				_dataIsComponent; \
		Type*				_dataType; \
		Bool*				_dataIsRequired; \
		/*void**			_values;*/   /* Storage space for param values */\
		/* ...dave - 22.01.09 */ \
		OperatorFunction_ApplyFunc*	applyMatrix; \
		OperatorFunction_ApplyFunc*	applyNeumann; \
		OperatorFunction_ApplyRHSFunc*	applyRHS; \
		Swarm*				integrationSwarm; \
		Swarm*				borderSwarm; \
		Index				rowPos; \
		Index				colPos; \
		StiffnessMatrix*		sm; \
		OperatorFunction_Datum*		_values;   /* Storage space for param values */\
		
	struct _OperatorFunction { __OperatorFunction };
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructor
	*/
	
	OperatorFunction*	OperatorFunction_New( 
					OperatorFunction_ApplyFunc*	applyMatrix, 
					OperatorFunction_ApplyFunc*	applyNeumann, 
					OperatorFunction_ApplyRHSFunc*	applyRHS, 
					Name				name,
					Index				dataCount,
					Name*				dataName,
					Bool*				dataIsComponent,
					Type*				dataType,
					Bool*				dataIsRequired,
					void**				values ); /* Assumes ownership */
	
	void			OperatorFunction_Init( 
					OperatorFunction*		self, 
					OperatorFunction_ApplyFunc*	applyMatrix, 
					OperatorFunction_ApplyFunc*	applyNeumann, 
					OperatorFunction_ApplyRHSFunc*	applyRHS, 
					/*Name				name,*/
					Index				dataCount,
					Name*				dataName,
					Bool*				dataIsComponent,
					Type*				dataType,
					Bool*				dataIsRequired,
					void**				values ); /* Assumes ownership */
	
	/*
	OperatorFunction*	_OperatorFunction_New( 
					SizeT				_sizeOfSelf, 
					Type				type,
					Stg_Class_DeleteFunction*	_delete,
					Stg_Class_PrintFunction*	_print, 
					Stg_Class_CopyFunction*		_copy, 
					OperatorFunction_ApplyFunc*	apply, 
					Name				name,
					Index				dataCount,
					Name*				dataName,
					Bool*				dataIsComponent,
					Type*				dataType,
					Bool*				dataIsRequired,
					void**				values );*/ /* Assumes ownership */
	OperatorFunction*	_OperatorFunction_New( 
					SizeT						_sizeOfSelf, 
					Type						type,
					Stg_Class_DeleteFunction*			_delete,
					Stg_Class_PrintFunction*			_print, 
					Stg_Class_CopyFunction*				_copy, 
					Stg_Component_DefaultConstructorFunction*	_defaultConstructor,
					Stg_Component_ConstructFunction*		_construct,
					Stg_Component_BuildFunction*			_build,
					Stg_Component_InitialiseFunction*		_initialise,
					Stg_Component_ExecuteFunction*			_execute,
					Stg_Component_DestroyFunction*			_destroy,
					OperatorFunction_ApplyFunc*			applyMatrix, 
					OperatorFunction_ApplyFunc*			applyNeumann, 
					OperatorFunction_ApplyRHSFunc*			applyRHS, 
					Index						dataCount,
					Name*						dataName,
					Bool*						dataIsComponent,
					Type*						dataType,
					Bool*						dataIsRequired,
					void**						values, /* Assumes ownership */
					Name						name );
	
	void			_OperatorFunction_Init(
					void*				operatorFunction, 
					OperatorFunction_ApplyFunc*	applyMatrix, 
					OperatorFunction_ApplyFunc*	applyNeumann, 
					OperatorFunction_ApplyRHSFunc*	applyRHS, 
					/*Name				name,*/
					Index				dataCount,
					Name*				dataName,
					Bool*				dataIsComponent,
					Type*				dataType,
					Bool*				dataIsRequired,
					void**				values ); /* Assumes ownership */
	
	/** The purpose of this constructor is to create a new OperatorFunction object based on information provided in by a 
	     StGermain meta file description for the operator. The provided name is the look up for meta details in the meta 
	     run-time database */
	OperatorFunction*	OperatorFunction_NewFromMeta( 
					OperatorFunction_ApplyFunc*	applyMatrix,
				        OperatorFunction_ApplyFunc*	applyNeumann,
					OperatorFunction_ApplyRHSFunc*	applyRHS,	
					Name				name );
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** General virtual functions
	*/

	void	_OperatorFunction_Delete( void* operatorFunction );
	
	void	_OperatorFunction_Print( void* operatorFunction, Stream* stream );
		
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/
	void*	_OperatorFunction_DefaultNew( Name name );
	void	_OperatorFunction_Construct( void* operatorFunction, Stg_ComponentFactory* cf, void* data );	
	void	_OperatorFunction_Build( void* operatorFunction, void* data );
	void	_OperatorFunction_Initialise( void* operatorFunction, void* data );
	void	_OperatorFunction_Execute( void* operatorFunction, void* data );
	void	_OperatorFunction_Destroy( void* operatorFunction, void* data );
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Build functions
	*/
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** As a Operatorfunction...
	*/


	/* Apply/run the operator. As macro for speed - anticipated use is in assembly inner loops */
	/*#define OperatorFunction_Apply( self, sm, lElement_I, elsm, param ) \
		(self)->apply( (sm), (lElement_I), (elsm), (param) )*/
	#define OperatorFunction_ApplyMatrix( self, sm, lElement_I, elsm, param ) \
		(self)->applyMatrix( (self), (sm), (lElement_I), (elsm), (param) )
	#define OperatorFunction_ApplyNeumann( self, sm, lElement_I, elsm, param ) \
		(self)->applyNeumann( (self), (sm), (lElement_I), (elsm), (param) )
	#define OperatorFunction_ApplyRHS( self, sm, lElement_I, elsm, param ) \
		(self)->applyRHS( (self), (sm), (lElement_I), (elsm), (param) )

	/** \fn OperatorFunction_Apply( void* operatorFunction, StiffnessMatrix* sm, Element_LocalIndex lElement_I, double** elsm, void* param )
	    \brief Apply/run the operator.
	*/
	
	/** The purpose of this function is to create the "data list" to pass to the operator function, based on the operator's 
	     definition in the input file/model. The data list items are raw values for use by the operator function, listed in the
	     order defined on creating this OperatorFunction (i.e. in calling OperatorFunction_New). Where this object is created
	     from a meta file's definition, the order is the order in the meta file, parameter then uses/dependancies. */
	OperatorFunction_Datum* OperatorFunction_CreateDataFromDictionary( 
		void*			operatorFunction, 
		Stg_ComponentFactory*	cf, 
		Dictionary*		dictionary );

#endif /* __StgFEM_Assembly_OperatorFunction_h__ */
