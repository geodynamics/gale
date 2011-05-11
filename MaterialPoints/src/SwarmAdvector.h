/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
** Copyright (c) 2005-2006, Monash Cluster Computing, Building 28, Monash University Clayton Campus,
**	Victoria, 3800, Australia
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**
** Contributors:
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Patrick D. Sunter, Software Engineer, VPAC. (patrick@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	David Stegman, Postdoctoral Fellow, Monash University. (david.stegman@sci.monash.edu.au)
**	Wendy Sharples, PhD Student, Monash University (wendy.sharples@sci.monash.edu.au)
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
** $Id: SwarmAdvector.h 374 2006-10-12 08:59:41Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __PICellerator_MaterialPoints_SwarmAdvector_h__
#define __PICellerator_MaterialPoints_SwarmAdvector_h__

	/* Textual name of this class */
	extern const Type SwarmAdvector_Type;

	/* SwarmAdvector information */
	#define __SwarmAdvector \
		/* General info */ \
		__TimeIntegrand \
		/* Virtual Info */\
		/* Other Info */\
		MaterialPointsSwarm*                  swarm;                \
		FeVariable*                           velocityField;        \
		PeriodicBoundariesManager*            periodicBCsManager;   \

	struct SwarmAdvector { __SwarmAdvector };
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/
	SwarmAdvector* SwarmAdvector_New(
		Name                                       name,
		DomainContext*                             context,
		TimeIntegrator*                            timeIntegrator,
		FeVariable*                                velocityField,
		Bool                                       allowFallbackToFirstOrder,
		MaterialPointsSwarm*                       swarm,
		PeriodicBoundariesManager*                 periodicBCsManager );

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define SWARMADVECTOR_DEFARGS \
                TIMEINTEGRAND_DEFARGS

	#define SWARMADVECTOR_PASSARGS \
                TIMEINTEGRAND_PASSARGS

	SwarmAdvector* _SwarmAdvector_New(  SWARMADVECTOR_DEFARGS  );

	void _SwarmAdvector_Init( 
		SwarmAdvector*                             self,
		FeVariable*                                velocityField,
		MaterialPointsSwarm*                       swarm,
		PeriodicBoundariesManager*                 periodicBCsManager );

	void _SwarmAdvector_Delete( void* materialSwarm );
	void _SwarmAdvector_Print( void* materialSwarm, Stream* stream );
	#define SwarmAdvector_Copy( self ) \
		(SwarmAdvector*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define SwarmAdvector_DeepCopy( self ) \
		(SwarmAdvector*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _SwarmAdvector_Copy( const void* materialSwarm, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	void* _SwarmAdvector_DefaultNew( Name name ) ;
void _SwarmAdvector_AssignFromXML( void* shape, Stg_ComponentFactory* cf, void* data ) ;
	void _SwarmAdvector_Build( void* materialSwarm, void* data ) ;
	void _SwarmAdvector_Initialise( void* materialSwarm, void* data ) ;
	void _SwarmAdvector_Execute( void* materialSwarm, void* data );
	void _SwarmAdvector_Destroy( void* materialSwarm, void* data ) ;
	Bool _SwarmAdvector_TimeDeriv( void* swarmAdvector, Index array_I, double* timeDeriv ) ;
	void _SwarmAdvector_Intermediate( void* swarmAdvector, Index array_I ) ;
	
		
	/*---------------------------------------------------------------------------------------------------------------------
	** Private functions
	*/
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Entry Point Hooks
	*/
	void SwarmAdvector_AdvectionSetup( TimeIntegrator* timeIntegrator, SwarmAdvector* self ) ;
	void SwarmAdvector_AdvectionFinish( TimeIntegrator* timeIntegrator, SwarmAdvector* self ) ;
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

#endif 

