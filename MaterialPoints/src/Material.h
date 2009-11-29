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

#ifndef __PICellerator_MaterialPoints_Materials_h__
#define __PICellerator_MaterialPoints_Materials_h__
	
	/* Textual name of this class */
	extern const Type Material_Type;
	
   extern const Index UNDEFINED_MATERIAL;

	/* Material information */
	#define __Material \
		__Stg_Component \
		\
		PICelleratorContext*	context; \
		Dictionary*				dictionary; \
		Stg_Shape*				shape; \
		Material_Index			index; /**< The index inside the Materials_Register */ \
		ExtensionManager*		extensionMgr; 
		
	struct Material { __Material };

	/* Public Constructor */
	Material* Material_New(
		Name						name,
		PICelleratorContext*	context, 
		Stg_Shape*				shape,
		Dictionary*				materialDictionary,
		Materials_Register*	materialRegister );

	void* _Material_DefaultNew( Name name );

	/* Private Constructor */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define MATERIAL_DEFARGS \
                STG_COMPONENT_DEFARGS

	#define MATERIAL_PASSARGS \
                STG_COMPONENT_PASSARGS

	Material* _Material_New(  MATERIAL_DEFARGS  );

	void _Material_AssignFromXML( void* material, Stg_ComponentFactory* cf, void* data );

	void _Material_Init( 
		void*						material, 	
		PICelleratorContext*	context, 
		Stg_Shape*				shape,
		Dictionary*				materialDictionary,
		Materials_Register*	materialRegister );

	/** Virtual Functions */
	void _Material_Delete( void* material );

	void _Material_Print( void* material, Stream* stream );

	void* _Material_Copy( void* material, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

	void _Material_Build( void* material, void* data );

	void _Material_Initialise( void* material, void* data );

	void _Material_Execute( void* material, void* data );

	void _Material_Destroy( void* material, void* data );

	/** Performs a layout of this material onto the points of the swarm, assigning the material index if the global
	 *  coordinates fall within the shape */
	void Material_Layout( void* material, MaterialPointsSwarm* swarm ) ;

	/** Calculates the volume and centroid of the material by using the weights of an integration swarm.
	 *     @param  centroid Output argument where the calculated centroid is written
	 *     @result The volume of the material. */
	double Material_Volume( void* material, IntegrationPointsSwarm* swarm, Coord centroid ) ;
	
	/** Integrates an FeVariable over a particular material (not the whole domain).
	 *  The function will also return the volume of the material in case the purpose of the integral is to get an average
	 *     @param volumeGlobal Output argument where the global volume is written.
	 *     @param result       Output argument where integration result is written. It will integrate each component of the 
	 *                         field and store it in result. result must be an array of as many doubles as there are components
	 *                         in the field. */
	void Material_IntegrateField( 
		void*                   material, 
		IntegrationPointsSwarm* swarm, 
		FeVariable*             field, 
		double*                 volumeGlobal, 
		double*                 result );

#endif 

