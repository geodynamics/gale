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

#ifndef __StgFEM_Discretisation_AnalyticSolution_h__
#define __StgFEM_Discretisation_AnalyticSolution_h__
	
	typedef void (AnalyticSolution_SolutionFunction) (void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* value );

	/** Textual name of this class */
	extern const Type AnalyticSolution_Type;
	
	/** AnalyticSolution information */
	#define __AnalyticSolution \
		/* Parent info */ \
		__Stg_Component \
		/* Virtual info */ \
		/* AnalyticSolution info */ \
		Stg_ObjectList*        feVariableList;             \
		Stg_ObjectList*        analyticFeVariableList;     \
		Stg_ObjectList*        analyticFeVariableFuncList; \
		Stg_ObjectList*        errorMagnitudeFieldList;    \
		Stg_ObjectList*        relativeErrorMagnitudeFieldList;    \
		Stg_ObjectList*        streamList;                 \
		double*                toleranceList;              \
		Swarm*                 integrationSwarm;           \
		LiveComponentRegister* LC_Register;                \
		DomainContext*         context;                    \
		AnalyticSolution_SolutionFunction* _getAnalyticVelocity; \
		AnalyticSolution_SolutionFunction* _getAnalyticPressure; \
		AnalyticSolution_SolutionFunction* _getAnalyticTotalStress; \
		AnalyticSolution_SolutionFunction* _getAnalyticStrainRate; \
		
	/** Brings together and manages the life-cycle of a a mesh and all the 
	info relevant to it for the Finite Element Method. See Mesh.h for more. */
	struct AnalyticSolution { __AnalyticSolution };
	
	/* --- Constructors / Destructor --- */
	
	/** Create a new AnalyticSolution and initialise */

	void* _AnalyticSolution_DefaultNew( Name name );
	
	/* Creation implementation / Virtual constructor */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define ANALYTICSOLUTION_DEFARGS \
                STG_COMPONENT_DEFARGS

	#define ANALYTICSOLUTION_PASSARGS \
                STG_COMPONENT_PASSARGS

	AnalyticSolution* _AnalyticSolution_New(  ANALYTICSOLUTION_DEFARGS  );			

	/* Stg_Class_Delete a AnalyticSolution construst */
	void _AnalyticSolution_Delete( void* analyticSolution );

	/* --- Virtual Function implementations --- */
	
	/* Print the contents of an AnalyticSolution construct */
	void _AnalyticSolution_Print( void* analyticSolution, Stream* stream );
	
	/* Copy */
	#define AnalyticSolution_Copy( self ) \
		(AnalyticSolution*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define AnalyticSolution_DeepCopy( self ) \
		(AnalyticSolution*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
	
	void* _AnalyticSolution_Copy( void* analyticSolution, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	/* Build implementation */
	void _AnalyticSolution_Build( void* analyticSolution, void* data );
	
	/* Construct implementation */
	void _AnalyticSolution_AssignFromXML( void* analyticSolution, Stg_ComponentFactory* cf, void* data );
	
	/* Initialisation implementation */
	void _AnalyticSolution_Initialise( void* analyticSolution, void* data );
	
	/* This function is called when the 'Update' phase happens */
	void AnalyticSolution_Update( void* analyticSolution ) ;

	/* Execution implementation */
	void _AnalyticSolution_Execute( void* analyticSolution, void* data );
	
	/* Destruct Implementation */
	void _AnalyticSolution_Destroy( void* analyticSolution, void* data );

	/* --- Public Functions --- */
	void AnalyticSolution_Test( void* analyticSolution, Index analyticFeVariable_I ) ;

	void AnalyticSolution_TestAll( void* analyticSolution, void* data ) ;
	
	void AnalyticSolution_PutAnalyticSolutionOntoNodes( void* analyticSolution, Index analyticFeVariable_I ) ;

	void AnalyticSolution_RegisterFeVariableWithAnalyticFunction(
		void*											analyticSolution,
		FeVariable*									feVariable,
		AnalyticSolution_SolutionFunction*	solutionFunction );

	FeVariable* AnalyticSolution_RegisterFeVariableFromCF(
		void*											analyticSolution,
		char*											fieldName,
		AnalyticSolution_SolutionFunction*	solutionFunction,
		Stg_ComponentFactory*					cf,
		Bool											isEssential,
		void*											data );

	void AnalyticSolution_BuildAllAnalyticFields( void* analyticSolution, void* data );

	FeVariable* AnalyticSolution_CreateAnalyticField( void* analyticSolution, FeVariable* feVariable ) ;

	FeVariable* AnalyticSolution_CreateAnalyticVectorField(
		void*											analyticSolution,
		FeVariable*									vectorField,
		AnalyticSolution_SolutionFunction*	solutionFunction ) ;

	FeVariable* AnalyticSolution_CreateAnalyticSymmetricTensorField( void* analyticSolution, FeVariable* vectorField ) ;

	FeVariable* AnalyticSolution_GetFeVariableFromAnalyticFeVariable( void* analyticSolution, FeVariable* analyticFeVariable ) ;

	InterpolationResult AnalyticSolution_InterpolateValueFromNormalFeVariable(
		void*			analyticSolution,
		FeVariable*	analyticFeVariable,
		double*		coord,
		double*		value ) ;

	AnalyticSolution* AnalyticSolution_GetAnalyticSolution();
	
#endif /* __StgFEM_Discretisation_AnalyticSolution_h__ */

