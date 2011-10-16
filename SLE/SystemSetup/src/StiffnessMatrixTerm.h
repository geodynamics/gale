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
** $Id: StiffnessMatrixTerm.h 733 2007-02-07 00:55:26Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_SLE_SystemSetup_StiffnessMatrixTerm_h__
#define __StgFEM_SLE_SystemSetup_StiffnessMatrixTerm_h__

	typedef void (StiffnessMatrixTerm_AssembleElementFunction)	(
			void*                             stiffnessMatrixTerm, 
			StiffnessMatrix*                  stiffnessMatrix, 
			Element_LocalIndex                lElement_I,
			SystemLinearEquations*            sle,
			FiniteElementContext*             context,
			double**                          elStiffMatToAdd );
	
	
	/* Textual name of this class */
	extern const Type StiffnessMatrixTerm_Type;
	
	/* StiffnessMatrix information */
	#define __StiffnessMatrixTerm  \
		/* General info */ \
		__Stg_Component \
		\
		FiniteElementContext*				     context;                  \
		/* Virtual info */ \
		StiffnessMatrixTerm_AssembleElementFunction*         _assembleElement;         \
		\
		/* General info */ \
		Stream*                                              debug;                    \
		Swarm*                                               integrationSwarm;         \
		Stg_Component*                                       extraInfo;                \
		StiffnessMatrix*                                     stiffnessMatrix;          \
		/* Data for GNx storage */ \
	        double                   **GNx; /* store globalDerivative ptr here */ \
	        double                   *N; /* store array for shape functions here */ \
		uint                      max_nElNodes;  /* holds the maxNumNodes per element */ 
	
	struct StiffnessMatrixTerm { __StiffnessMatrixTerm };
	
	/* Creation implementation / Virtual constructor */
	StiffnessMatrixTerm* StiffnessMatrixTerm_New(
		Name                                                 name,
		FiniteElementContext*				                    context,
		StiffnessMatrix*                                     stiffnessMatrix,
		Swarm*                                               integrationSwarm,
		Stg_Component*                                       extraInfo );

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define STIFFNESSMATRIXTERM_DEFARGS \
                STG_COMPONENT_DEFARGS, \
                StiffnessMatrixTerm_AssembleElementFunction*  _assembleElement

	#define STIFFNESSMATRIXTERM_PASSARGS \
                STG_COMPONENT_PASSARGS, \
	        _assembleElement

	StiffnessMatrixTerm* _StiffnessMatrixTerm_New(  STIFFNESSMATRIXTERM_DEFARGS  );
	
	/* 'Stg_Class' Virtual Functions */
	void _StiffnessMatrixTerm_Delete( void* stiffnessMatrixTerm );
	void _StiffnessMatrixTerm_Print( void* stiffnessMatrixTerm, Stream* stream );
	#define StiffnessMatrixTerm_Copy( self ) \
		(StiffnessMatrixTerm*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define StiffnessMatrixTerm_DeepCopy( self ) \
		(StiffnessMatrixTerm*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _StiffnessMatrixTerm_Copy( const void* stiffnessMatrixTerm, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	/* 'Stg_Component' Virtual Functions */
	void* _StiffnessMatrixTerm_DefaultNew( Name name );
	void _StiffnessMatrixTerm_AssignFromXML( void* stiffnessMatrixTerm, Stg_ComponentFactory* cf, void* data );
	void _StiffnessMatrixTerm_Build( void* stiffnessMatrixTerm, void* data );
	void _StiffnessMatrixTerm_Initialise( void* stiffnessMatrixTerm, void* data );
	void _StiffnessMatrixTerm_Execute( void* stiffnessMatrixTerm, void* data );
	void _StiffnessMatrixTerm_Destroy( void* stiffnessMatrixTerm, void* data );

	void _StiffnessMatrixTerm_Init(
		void*                                                stiffnessMatrixTerm,
		FiniteElementContext*				                    context,
		StiffnessMatrix*                                     stiffnessMatrix,
		Swarm*                                               integrationSwarm,
		Stg_Component*                                       extraInfo );
	
	void StiffnessMatrixTerm_AssembleElement( 
			void*                             stiffnessMatrixTerm, 
			StiffnessMatrix*                  stiffnessMatrix, 
			Element_LocalIndex                lElement_I,
			SystemLinearEquations*            sle,
			FiniteElementContext*             context,
			double**                          elStiffMatToAdd );

	void _StiffnessMatrixTerm_AssembleElement( 
			void*                             stiffnessMatrixTerm, 
			StiffnessMatrix*                  stiffnessMatrix, 
			Element_LocalIndex                lElement_I,
			SystemLinearEquations*            sle,
			FiniteElementContext*             context,
			double**                          elStiffMatToAdd ) ;

	void StiffnessMatrixTerm_SetAssembleElementFunction( void* stiffnessMatrixTerm, StiffnessMatrixTerm_AssembleElementFunction* assembleElementFunction ) ;

#endif /* __StgFEM_SLE_SystemSetup_StiffnessMatrixTerm_h__ */

