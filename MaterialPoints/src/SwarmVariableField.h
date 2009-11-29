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


#ifndef __PICellerator_Utils_SwarmVariableField_h__
#define __PICellerator_Utils_SwarmVariableField_h__

	/** Textual name of this class */
	extern const Type SwarmVariableField_Type;

	typedef void (SwarmVariableField_GetValueAtNodeFunction) ( void* swarmVariableField, Node_DomainIndex dNode_I, double* result );
	typedef void (SwarmVariableField_ValueAtParticleFunction) ( void* swarmVariableField, Particle_Index lParticle_I, double* result );

	/** SwarmVariableField class contents */
	#define __SwarmVariableField \
		/* General info */ \
		__ParticleFeVariable; \
		/* Virtual info */ \
		/* SwarmVariableField info */ \
		unsigned		dofCount;     		\
		MaterialPointsSwarm*	materialSwarm;  	\
		Name			swarmVarName;		\
		SwarmVariable*		swarmVar;     		\
		Variable_Register*	variable_Register;	\

	struct SwarmVariableField { __SwarmVariableField };
/*	
	SwarmVariableField* SwarmVariableField_New( 
		Name                                                name,
		MaterialPointsSwarm*                                swarm,
		Index                                               dofCount );
*/	
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define SWARMVARIABLEFIELD_DEFARGS \
                PARTICLEFEVARIABLE_DEFARGS

	#define SWARMVARIABLEFIELD_PASSARGS \
                PARTICLEFEVARIABLE_PASSARGS

	SwarmVariableField* _SwarmVariableField_New(  SWARMVARIABLEFIELD_DEFARGS  );
	
	void _SwarmVariableField_Delete( void* swarmVariable );
	void _SwarmVariableField_Print( void* swarmVariable, Stream* stream );

	void* _SwarmVariableField_DefaultNew( Name name ) ;
	void _SwarmVariableField_AssignFromXML( void* swarmVariable, Stg_ComponentFactory* cf, void* data ) ;
	void _SwarmVariableField_Build( void* swarmVariable, void* data ) ;
	void _SwarmVariableField_Initialise( void* swarmVariable, void* data ) ;
	void _SwarmVariableField_Execute( void* swarmVariable, void* data ) ;
	void _SwarmVariableField_Destroy( void* swarmVariable, void* data ) ;

	void _SwarmVariableField_ValueAt( void* swarmVariable, Particle_Index lParticle_I, double* value ) ;
	/*
	double _SwarmVariableField_GetMinGlobalMagnitude( void* swarmVariable ) ;
	double _SwarmVariableField_GetMaxGlobalMagnitude( void* swarmVariable ) ;
	*/
	void _SwarmVariableField_GetValueAtNode( void* swarmVariableField, Node_DomainIndex dNode_I, double* value );

	void _SwarmVariableField_ValueAtParticle( void* swarmVariableField, IntegrationPointsSwarm* swarm, Element_LocalIndex lElement_I, /*Particle_Index lParticle_I*/ IntegrationPoint* particle, double* value );
	
#endif

