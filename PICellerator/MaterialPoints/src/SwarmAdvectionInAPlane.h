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
**	Special swarm advector requested by Dave Stegman (hence AdvectorD) which doesn't advect in the K axis, even
**	for 3D problems.
**
** Assumptions:
**
** Comments:
**
** $Id: SwarmAdvectionInAPlane.h 189 2005-10-20 00:39:29Z RobertTurnbull $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __PICellerator_MaterialPoints_SwarmAdvectionInAPlane_h__
#define __PICellerator_MaterialPoints_SwarmAdvectionInAPlane_h__

	/* Textual name of this class */
	extern const Type SwarmAdvectionInAPlane_Type;

	/* SwarmAdvectionInAPlane information */
	#define __SwarmAdvectionInAPlane \
		/* General info */ \
		__SwarmAdvector 	/** Now inherits from SwarmAdvector class */ \
		int				whichaxis;								
		
	struct SwarmAdvectionInAPlane { __SwarmAdvectionInAPlane };
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/
	SwarmAdvectionInAPlane* SwarmAdvectionInAPlane_New(
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

	#define SWARMADVECTIONINAPLANE_DEFARGS \
                SWARMADVECTOR_DEFARGS

	#define SWARMADVECTIONINAPLANE_PASSARGS \
                SWARMADVECTOR_PASSARGS

	SwarmAdvectionInAPlane* _SwarmAdvectionInAPlane_New(  SWARMADVECTIONINAPLANE_DEFARGS  );

	void _SwarmAdvectionInAPlane_Init( 
		SwarmAdvectionInAPlane*                             self,
		int											whichaxis);

	#define SwarmAdvectionInAPlane_Copy( self ) \
		(SwarmAdvectionInAPlane*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define SwarmAdvectionInAPlane_DeepCopy( self ) \
		(SwarmAdvectionInAPlane*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _SwarmAdvectionInAPlane_Copy( const void* materialSwarm, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	void* _SwarmAdvectionInAPlane_DefaultNew( Name name ) ;
void _SwarmAdvectionInAPlane_AssignFromXML( void* shape, Stg_ComponentFactory* cf, void* data ) ;
	void _SwarmAdvectionInAPlane_Build( void* materialSwarm, void* data ) ;
	void _SwarmAdvectionInAPlane_Initialise( void* materialSwarm, void* data ) ;
	void _SwarmAdvectionInAPlane_Execute( void* materialSwarm, void* data );
	void _SwarmAdvectionInAPlane_Destroy( void* materialSwarm, void* data ) ;
	Bool _SwarmAdvectionInAPlane_TimeDeriv( void* swarmAdvector, Index array_I, double* timeDeriv ) ;
	void _SwarmAdvectionInAPlane_Intermediate( void* swarmAdvector, Index array_I ) ;
	
		
	/*---------------------------------------------------------------------------------------------------------------------
	** Private functions
	*/
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Entry Point Hooks
	*/
	void SwarmAdvectionInAPlane_AdvectionSetup( TimeIntegrator* timeIntegrator, SwarmAdvectionInAPlane* self ) ;
	void SwarmAdvectionInAPlane_AdvectionFinish( TimeIntegrator* timeIntegrator, SwarmAdvectionInAPlane* self ) ;
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

#endif 

