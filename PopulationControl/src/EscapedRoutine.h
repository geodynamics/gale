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
** $Id: EscapedRoutine.h 189 2005-10-20 00:39:29Z RobertTurnbull $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __PICellerator_PopulationControl_EscapedRoutine_h__
#define __PICellerator_PopulationControl_EscapedRoutine_h__

	typedef void (EscapedRoutine_SelectFunction)( void* escapedRoutine, void* _swarm );

	/* Textual name of this class */
	extern const Type EscapedRoutine_Type;

	/* EscapedRoutine information */
	#define __EscapedRoutine \
		/* General info */ \
		__Stg_Component \
		DomainContext*				context; \
		/* Virtual Info */\
		EscapedRoutine_SelectFunction*		_select; \
		/* Other Info */\
		Stream*                                    debug;                    \
		Dimension_Index                            dim;                      \
		/* Removal Info */  \
		Particle_Index                             particlesToRemoveCount;   \
		Particle_Index                             particlesToRemoveAlloced; \
		Particle_Index                             particlesToRemoveDelta;   \
		unsigned*				   particlesToRemoveList;    

	struct EscapedRoutine { __EscapedRoutine };

	/*---------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define ESCAPEDROUTINE_DEFARGS \
                STG_COMPONENT_DEFARGS, \
                EscapedRoutine_SelectFunction*  _select

	#define ESCAPEDROUTINE_PASSARGS \
                STG_COMPONENT_PASSARGS, \
	        _select

	EscapedRoutine* _EscapedRoutine_New(  ESCAPEDROUTINE_DEFARGS  );

	void* _EscapedRoutine_DefaultNew( Name name );
	
	/* Stg_Class_Delete EscapedRoutine implementation */
	void _EscapedRoutine_Delete( void* escapedRoutine );
	void _EscapedRoutine_Print( void* escapedRoutine, Stream* stream );
	#define EscapedRoutine_Copy( self ) \
		(EscapedRoutine*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define EscapedRoutine_DeepCopy( self ) \
		(EscapedRoutine*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _EscapedRoutine_Copy( void* escapedRoutine, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	void _EscapedRoutine_AssignFromXML( void* shape, Stg_ComponentFactory* cf, void* data ) ;
	void _EscapedRoutine_Build( void* escapedRoutine, void* data ) ;
	void _EscapedRoutine_Initialise( void* escapedRoutine, void* data ) ;
	void _EscapedRoutine_Execute( void* escapedRoutine, void* data );
	void _EscapedRoutine_Destroy( void* escapedRoutine, void* data ) ;
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Public member functions
	*/

	void EscapedRoutine_Select( void* escapedRoutine, void* _swarm ) ;
	void _EscapedRoutine_Select( void* escapedRoutine, void* _swarm );
	void EscapedRoutine_RemoveFromSwarm( void* escapedRoutine, void* _swarm ) ;

	void EscapedRoutine_InitialiseParticleList( void* escapedRoutine ) ;
	void EscapedRoutine_SetParticleToRemove( void* escapedRoutine, Swarm* swarm, Particle_Index lParticle_I ) ;

	void EscapedRoutine_SortParticleList( void* escapedRoutine ) ;
	void EscapedRoutine_RemoveParticles( void* escapedRoutine, Swarm* swarm ) ;

#endif 

