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
** $Id: ForceTerm.h 656 2006-10-18 06:45:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_SLE_SystemSetup_ForceTerm_h__
#define __StgFEM_SLE_SystemSetup_ForceTerm_h__

	typedef void (ForceTerm_AssembleElementFunction)	(
		void*						forceTerm, 
		ForceVector*			forceVector, 
		Element_LocalIndex	lElement_I,
		double*					elForceVecToAdd );
	
	
	/* Textual name of this class */
	extern const Type ForceTerm_Type;
	
	/* StiffnessMatrix information */
	#define __ForceTerm  \
		/* General info */ \
		__Stg_Component \
		\
		FiniteElementContext*					context; \
		/* Virtual info */ \
		ForceTerm_AssembleElementFunction*	_assembleElement; \
		\
		/* General info */ \
		Stream*										debug; \
		Swarm*										integrationSwarm; \
		Stg_Component*								extraInfo;
	
	struct ForceTerm { __ForceTerm };


	
	/* Creation implementation / Virtual constructor */
	ForceTerm* ForceTerm_New(
		Name							name,
		FiniteElementContext*	context,
		ForceVector*				forceVector,
		Swarm*						integrationSwarm,
		Stg_Component*				extraInfo );

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define FORCETERM_DEFARGS \
                STG_COMPONENT_DEFARGS, \
                ForceTerm_AssembleElementFunction*  _assembleElement

	#define FORCETERM_PASSARGS \
                STG_COMPONENT_PASSARGS, \
	        _assembleElement

	ForceTerm* _ForceTerm_New(  FORCETERM_DEFARGS  );
	
	void _ForceTerm_Init(
		void*							forceTerm,
		FiniteElementContext*	context,
		ForceVector*				forceVector,
		Swarm*						integrationSwarm,
		Stg_Component*				extraInfo );

	/* 'Stg_Class' Virtual Functions */
	void _ForceTerm_Delete( void* forceTerm );

	void _ForceTerm_Print( void* forceTerm, Stream* stream );

	#define ForceTerm_Copy( self ) \
		(ForceTerm*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define ForceTerm_DeepCopy( self ) \
		(ForceTerm*)Stg_Class_Copy( self, NULL, True, NULL, NULL )

	void* _ForceTerm_Copy( void* forceTerm, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	/* 'Stg_Component' Virtual Functions */
	void* _ForceTerm_DefaultNew( Name name );

	void _ForceTerm_AssignFromXML( void* forceTerm, Stg_ComponentFactory* cf, void* data );

	void _ForceTerm_Build( void* forceTerm, void* data );

	void _ForceTerm_Initialise( void* forceTerm, void* data );

	void _ForceTerm_Execute( void* forceTerm, void* data );

	void _ForceTerm_Destroy( void* forceTerm, void* data );
	
	void ForceTerm_AssembleElement( 
		void*						forceTerm, 
		ForceVector*			forceVector, 
		Element_LocalIndex	lElement_I,
		double*					elForceVecToAdd );

	void _ForceTerm_AssembleElement( 
		void*						forceTerm, 
		ForceVector*			forceVector, 
		Element_LocalIndex	lElement_I,
		double*					elForceVecToAdd ) ;

	void ForceTerm_SetAssembleElementFunction( void* forceTerm, ForceTerm_AssembleElementFunction* assembleElementFunction ) ;

#endif /* __StgFEM_SLE_SystemSetup_ForceTerm_h__ */

