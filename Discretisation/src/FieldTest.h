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
**
** Assumptions:
**
** Comments:
**
** $Id $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_Discretisation_FieldTest_h__
#define __StgFEM_Discretisation_FieldTest_h__
	
	typedef void (FieldTest_AnalyticSolutionFunc) (void* fieldTest, double* coord, double* value );

	/** Textual name of this class */
	extern const Type FieldTest_Type;
	
	/** FieldTest information */
	#define __FieldTest \
		/* Parent info */ \
		__Stg_Component \
		/* Virtual info */ \
		/* FieldTest info */ \
		FeVariable**				numericFieldList;	\
		FeVariable**				referenceFieldList;	\
		FeVariable**				errorFieldList;		\
		OperatorFeVariable**			referenceMagFieldList;	\
		OperatorFeVariable**			errorMagFieldList;	\
		unsigned				fieldCount;		\
		Swarm*					numericSwarm;		\
		Name*					swarmVarNameList;	\
		unsigned				swarmVarCount;		\
		FeMesh*					constantMesh;		\
		FeMesh*					elementMesh;		\
		Bool					normalise;		\
		double**				gAnalyticSq;		\
		double**				gErrorSq;		\
		double**				gError;			\
		double**				gErrorNorm;		\
		Name*					referenceSolnFileList;	\
		Name					referenceSolnPath;	\
		Swarm*                 			integrationSwarm;       \
		DomainContext*				context;		\
		LiveComponentRegister*			LCRegister;		\
		Bool					referenceSolnFromFile;	\
		/* must fill this array in the analytic solution plugin so */	\
		/* that the correct analytic solution func is applied for the*/ \
		/* correct feVariable */			 		\
		Index*					analyticSolnForFeVarKey;\
		FieldTest_AnalyticSolutionFunc**	_analyticSolutionList;	\
		/* store any constants required by the analytic solution */     \
		/* plugin in this array (to be initialised in the plugin) */    \
		void**					data;			\
		

	/** Brings together and manages the life-cycle of a a mesh and all the 
	info relevant to it for the Finite Element Method. See Mesh.h for more. */
	struct FieldTest { __FieldTest };
	
	/* --- Constructors / Destructor --- */
	
	/** Create a new FieldTest and initialise */

	void* _FieldTest_DefaultNew( Name name );
	
	/* Creation implementation / Virtual constructor */
	FieldTest* _FieldTest_New(
		SizeT                                       _sizeOfSelf,
		Type                                        type,
		Stg_Class_DeleteFunction*                   _delete,
		Stg_Class_PrintFunction*                    _print,
		Stg_Class_CopyFunction*                     _copy, 
		Stg_Component_DefaultConstructorFunction*   _defaultConstructor,
		Stg_Component_ConstructFunction*            _construct,
		Stg_Component_BuildFunction*                _build,
		Stg_Component_InitialiseFunction*           _initialise,
		Stg_Component_ExecuteFunction*              _execute,
		Stg_Component_DestroyFunction*              _destroy,
		Name                                        name );			

	/* Stg_Class_Delete a FieldTest construst */
	void _FieldTest_Delete( void* fieldTest );

	/* --- Virtual Function implementations --- */
	
	/* Print the contents of an FieldTest construct */
	void _FieldTest_Print( void* fieldTest, Stream* stream );
	
	/* Copy */
	#define FieldTest_Copy( self ) \
		(FieldTest*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define FieldTest_DeepCopy( self ) \
		(FieldTest*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
	
	void* _FieldTest_Copy( void* fieldTest, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	/* Build implementation */
	void _FieldTest_Build( void* fieldTest, void* data );
	
	/* Construct implementation */
	void _FieldTest_Construct( void* fieldTest, Stg_ComponentFactory* cf, void* data );
	
	/* Initialisation implementation */
	void _FieldTest_Initialise( void* fieldTest, void* data );
	
	/* This function is called when the 'Update' phase happens */
	void FieldTest_Update( void* fieldTest ) ;

	/* Execution implementation */
	void _FieldTest_Execute( void* fieldTest, void* data );
	
	/* Destruct Implementation */
	void _FieldTest_Destroy( void* fieldTest, void* data );

	/* --- Public Functions --- */
	void FieldTest_BuildReferenceField( void* fieldTest, Index field_I );
	void FieldTest_BuildErrField( void* fieldTest, Index field_I );

	void FieldTest_LoadReferenceSolutionFromFile( FeVariable* referenceField, Name referenceSolnName, Name referenceSolnPath, DomainContext* context );
	void FieldTest_CalculateAnalyticSolutionForField( void* fieldTest, Index field_I );

	void FieldTest_GenerateErrFields( void* fieldTest, void* data );

	void FieldTest_ElementErrAnalyticFromField( void* fieldTest, Index field_I, Index lElement_I, double* elErrorSq, double* elNormSq );
	void FieldTest_ElementErrReferenceFromField( void* fieldTest, Index field_I, Index lElement_I, double* elErrorSq, double* elNormSq );
	void FieldTest_ElementErrAnalyticFromSwarm( void* fieldTest, Index field_I, Index lElement_I, double* elErrorSq, double* elNormSq );
	void FieldTest_ElementErrReferenceFromSwarm( void* fieldTest, Index field_I, Index lElement_I, double* elErrorSq, double* elNormSq );

	void FieldTest_AddAnalyticSolutionFuncToListAtIndex( void* fieldTest, Index func_I, FieldTest_AnalyticSolutionFunc* func, Index field_I );
#endif /* __StgFEM_Discretisation_FieldTest_h__ */
