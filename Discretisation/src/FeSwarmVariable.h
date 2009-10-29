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
**	$Id: FeSwarmVariable.h 654 2006-10-12 08:58:49Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Discretisation_FeSwarmVariable_h__
#define __Discretisation_FeSwarmVariable_h__

	/** Textual name of this class */
	extern const Type FeSwarmVariable_Type;

	/** FeSwarmVariable contents */
	#define __FeSwarmVariable \
		/* General info */ \
		__SwarmVariable \
		\
		/* Virtual info */ \
		FeVariable*                                                 feVariable;

	struct FeSwarmVariable { __FeSwarmVariable };	

	/* Public Constructor */

	/** Private Constructor */
	FeSwarmVariable* _FeSwarmVariable_New(
 		SizeT                                              _sizeOfSelf, 
		Type                                               type,
		Stg_Class_DeleteFunction*                          _delete,
		Stg_Class_PrintFunction*                           _print, 
		Stg_Class_CopyFunction*                            _copy, 
		Stg_Component_DefaultConstructorFunction*          _defaultConstructor,
		Stg_Component_ConstructFunction*                   _construct,
		Stg_Component_BuildFunction*                       _build,
		Stg_Component_InitialiseFunction*                  _initialise,
		Stg_Component_ExecuteFunction*                     _execute,
		Stg_Component_DestroyFunction*                     _destroy,
		SwarmVariable_ValueAtFunction*                     _valueAt,
		SwarmVariable_GetGlobalValueFunction*              _getMinGlobalMagnitude,
		SwarmVariable_GetGlobalValueFunction*              _getMaxGlobalMagnitude,
		Name                                               name );

	/* 'Stg_Class' Virtual Implementations */
	void _FeSwarmVariable_Delete( void* variable ) ;
	void _FeSwarmVariable_Print( void* _swarmVariable, Stream* stream ) ;
	void* _FeSwarmVariable_Copy( void* swarmVariable, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) ;
	
	/* 'Stg_Component' Virtual Implementations */
	void* _FeSwarmVariable_DefaultNew( Name name );
void _FeSwarmVariable_AssignFromXML( void* swarmVariable, Stg_ComponentFactory* cf, void* data ) ;
	void _FeSwarmVariable_Build( void* swarmVariable, void* data ) ;
	void _FeSwarmVariable_Execute( void* variable, void* data ) ;
	void _FeSwarmVariable_Destroy( void* variable, void* data ) ;
	void _FeSwarmVariable_Initialise( void* variable, void* data ) ;

	/* 'SwarmVariable Virtual Implementations */
	void _FeSwarmVariable_ValueAt( void* swarmVariable, Particle_Index lParticle_I, double* value ) ;
	double _FeSwarmVariable_GetMinGlobalMagnitude( void* swarmVariable ) ;
	double _FeSwarmVariable_GetMaxGlobalMagnitude( void* swarmVariable ) ;
	
	/** Private Functions */


#endif /* __Discretisation_FeSwarmVariable_h__ */
