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
 **     An IntegrationPointMapper which maps many MaterialPointsSwarms to one IntegrationPointsSwarm
 **
 ** Assumptions:
 **
 ** Comments:
 **
 ** $Id: ManyToOneMapper.h 189 2005-10-20 00:39:29Z RobertTurnbull $
 **
 **~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __PICellerator_MaterialPoints_ManyToOneMapper_h__
#define __PICellerator_MaterialPoints_ManyToOneMapper_h__

	/* Textual name of this class */
	extern const Type ManyToOneMapper_Type;

	/* ManyToOneMapper information */
	#define __ManyToOneMapper \
		__IntegrationPointMapper \
		\
		MaterialPointsSwarm**	materialSwarms; \
		Index							materialSwarmCount;

	struct ManyToOneMapper { __ManyToOneMapper };
	


	/*---------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define MANYTOONEMAPPER_DEFARGS \
                INTEGRATIONPOINTMAPPER_DEFARGS

	#define MANYTOONEMAPPER_PASSARGS \
                INTEGRATIONPOINTMAPPER_PASSARGS

	ManyToOneMapper* _ManyToOneMapper_New(  MANYTOONEMAPPER_DEFARGS  );

	void _ManyToOneMapper_Init( void* mapper, MaterialPointsSwarm** materialSwarms, Index materialSwarmCount );

	void _ManyToOneMapper_Delete( void* mapper );

	void _ManyToOneMapper_Print( void* mapper, Stream* stream );

	#define ManyToOneMapper_Copy( self ) \
		(ManyToOneMapper*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define ManyToOneMapper_DeepCopy( self ) \
		(ManyToOneMapper*) Stg_Class_Copy( self, NULL, True, NULL, NULL )

	void* _ManyToOneMapper_Copy( void* mapper, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	void* _ManyToOneMapper_DefaultNew( Name name );

	void _ManyToOneMapper_AssignFromXML( void* shape, Stg_ComponentFactory* cf, void* data );

	void _ManyToOneMapper_Build( void* mapper, void* data );

	void _ManyToOneMapper_Initialise( void* mapper, void* data );

	void _ManyToOneMapper_Execute( void* mapper, void* data );

	void _ManyToOneMapper_Destroy( void* mapper, void* data );

	MaterialPointsSwarm** ManyToOneMapper_GetMaterialPointsSwarms( void* mapper, Index* count );

#endif

