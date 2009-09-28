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
**	Allows users to access Materialss based on their textual name,
**	or index.
**
** Assumptions:
**
** Comments:
**
**
** $Id: Materials_Register.h 280 2006-04-13 04:22:21Z AlanLo $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __PICellerator_MaterialPoints_Materials_Register_h__
#define __PICellerator_MaterialPoints_Materials_Register_h__
	
	
	extern const Type Materials_Register_Type;
	
	#define __Materials_Register \
		__NamedObject_Register

	struct Materials_Register { __Materials_Register };
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructor
	*/
	
	Materials_Register* Materials_Register_New( void );
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** General virtual functions
	*/
	void _Materials_Register_Delete( void* _materialsRegister ) ;
	void _Materials_Register_Print( void* _materialsRegister, Stream* stream ) ;
	void* _Materials_Register_Copy( void* _materialsRegister, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/
	#define Materials_Register_Add NamedObject_Register_Add
	#define Materials_Register_GetIndex NamedObject_Register_GetIndex 

	#define Materials_Register_GetByName( self, materialName ) \
		( (Material*) NamedObject_Register_GetByName( self, materialName ) ) 
	
	Material* Materials_Register_GetByIndex( Materials_Register* self, Index materialIndex );

	#define Materials_Register_GetCount( self ) \
		self->objects->count
	#define Materials_Register_PrintAllEntryNames NamedObject_Register_PrintAllEntryNames

	/** Setups material information for the given swarm */
	void Materials_Register_SetupSwarm( void* materialRegister, MaterialPointsSwarm* swarm );

	/** Adds an extension to each known material in the system.*/
	ExtensionInfo_Index Materials_Register_AddMaterialExtension( void* materialsRegister, Type type, SizeT extensionSize ) ;

	/** Assigns a material to each point of the swarm based on geometry setup of each material */	
	void _Materials_Register_LayoutGeometry( void* materialsRegister, void* _swarm );
	
	/** Assigns material property extension values to each material point based on the material they are given by using
	 *  the property variables.  */
	void Materials_Register_AssignParticleProperties( 
			void*                   materialRegister, 
			MaterialPointsSwarm*    swarm, 
			Variable_Register*      variableRegister );


	/** Assigns the value for the index a number taken from the dictionary with a key that has same name as the variable */
	void Variable_SetValueFromDictionary( void* _variable, Index index, Dictionary* dictionary );

	/** For the given index in each variable of the register, assign a numerical value from the dictionary with a key that has 
	 *  the same name as the variable. */
	void Variable_Register_SetAllVariablesFromDictionary( void* _variable_Register, Index index, Dictionary* dictionary );

	/** Private function useful for printing status updates during property assigning process. */
	void _Materials_Register_PrintParticleAssignUpdate(  
		void*                   materialRegister,
		MaterialPointsSwarm*    swarm,
		Particle_Index          lParticle_I,
		Stream*                 stream,
		Bool*                   firstStatusPrint );

#endif

