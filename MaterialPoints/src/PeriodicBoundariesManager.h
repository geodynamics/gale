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
**	Manages the interactions between particles and periodic boundary conditions.
**
** Assumptions:
**
** Comments:
**
** $Id: PeriodicBoundariesManager.h 456 2007-04-27 06:21:01Z LukeHodkinson $
*
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __PICellerator_PeriodicBoundaries_PeriodicBoundariesManager_h__
#define __PICellerator_PeriodicBoundaries_PeriodicBoundariesManager_h__
	
	/* Textual name of this class */
	extern const Type PeriodicBoundariesManager_Type;

	typedef struct PeriodicBoundary {
		Axis				axis; /* Which plane the BC is in */
		double			minWall;
		double			maxWall;
		unsigned int	particlesUpdatedMinEndCount;
		unsigned int	particlesUpdatedMaxEndCount;
	} PeriodicBoundary;

	#define __PeriodicBoundariesManager \
		__Stg_Component \
		PICelleratorContext*	context; \
		\
		Dictionary*				dictionary; \
		Mesh*						mesh; \
		Index						count; \
		Index						size; \
		Index						delta; \
		PeriodicBoundary*		boundaries; \
		Swarm*					swarm; \
		Stream*					debug; 

	struct PeriodicBoundariesManager { __PeriodicBoundariesManager };



	void* _PeriodicBoundariesManager_DefaultNew( Name name );

	PeriodicBoundariesManager* PeriodicBoundariesManager_New( 
		Name						name,
		PICelleratorContext*	context,
		Mesh*						mesh, 
		Swarm*					swarm,
		Dictionary*				dictionary );

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define PERIODICBOUNDARIESMANAGER_DEFARGS \
                STG_COMPONENT_DEFARGS

	#define PERIODICBOUNDARIESMANAGER_PASSARGS \
                STG_COMPONENT_PASSARGS

	PeriodicBoundariesManager* _PeriodicBoundariesManager_New(  PERIODICBOUNDARIESMANAGER_DEFARGS  );

	void _PeriodicBoundariesManager_Init(
		void*						periodicBCsManager,
		PICelleratorContext*	context,
		Mesh*						mesh, 
		Swarm*					swarm,
		Dictionary*				dictionary );
		
	void _PeriodicBoundariesManager_AssignFromXML( void* periodicBCsManager, Stg_ComponentFactory* cf, void* data );
	
	void _PeriodicBoundariesManager_Delete( void* context );

	void _PeriodicBoundariesManager_Print( void* context, Stream* stream );

	void* _PeriodicBoundariesManager_Copy( void* periodicBCsManager, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

	void _PeriodicBoundariesManager_Build( void* periodicBCsManager, void* data );

	void _PeriodicBoundariesManager_Initialise( void* periodicBCsManager, void* data );

	void _PeriodicBoundariesManager_Execute( void* periodicBCsManager, void* data );

	void _PeriodicBoundariesManager_Destroy( void* periodicBCsManager, void* data );

	void PeriodicBoundariesManager_AddPeriodicBoundary( void* periodicBCsManager, Axis axis );

	void PeriodicBoundariesManager_UpdateParticle( void* periodicBCsManager, Particle_Index lParticle_I );

#endif

