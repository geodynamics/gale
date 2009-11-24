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
 **	A OneToOne mapping between MaterialPointsSwarm and IntegrationPointsSwarm where the translation a direct
 **      global to local coordinates 
 **
 ** Assumptions:
 **
 ** Comments:
 **
 ** $Id: GaussCoincidentMapper.h 189 2005-10-20 00:39:29Z RobertTurnbull $
 **
 **~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __PICellerator_MaterialPoints_GaussCoincidentMapper_h__
#define __PICellerator_MaterialPoints_GaussCoincidentMapper_h__

	/* Textual name of this class */
	extern const Type GaussCoincidentMapper_Type;

	#define __GaussCoincidentMapper \
		__OneToOneMapper 

	struct GaussCoincidentMapper { __GaussCoincidentMapper };

	#define GAUSSCOINCIDENTMAPPER_DEFARGS \
		ONETOONEMAPPER_DEFARGS

	#define GAUSSCOINCIDENTMAPPER_PASSARGS \
		ONETOONEMAPPER_PASSARGS

	GaussCoincidentMapper* GaussCoincidentMapper_New(
		Name							name,
		PICelleratorContext*		context,
		IntegrationPointsSwarm*	integrationSwarm,
		MaterialPointsSwarm*		materialSwarm);	
	
	GaussCoincidentMapper* _GaussCoincidentMapper_New( GAUSSCOINCIDENTMAPPER_DEFARGS );

	void _GaussCoincidentMapper_Init( void* mapper );

	void _GaussCoincidentMapper_Delete( void* mapper );

	void _GaussCoincidentMapper_Print( void* mapper, Stream* stream );

	#define GaussCoincidentMapper_Copy( self ) \
		(GaussCoincidentMapper*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define GaussCoincidentMapper_DeepCopy( self ) \
		(GaussCoincidentMapper*) Stg_Class_Copy( self, NULL, True, NULL, NULL )

	void* _GaussCoincidentMapper_Copy( void* mapper, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	void* _GaussCoincidentMapper_DefaultNew( Name name );

	void _GaussCoincidentMapper_AssignFromXML( void* shape, Stg_ComponentFactory* cf, void* data );

	void _GaussCoincidentMapper_Build( void* mapper, void* data );

	void _GaussCoincidentMapper_Initialise( void* mapper, void* data );

	void _GaussCoincidentMapper_Execute( void* mapper, void* data );

	void _GaussCoincidentMapper_Destroy( void* mapper, void* data );
		
	void _GaussCoincidentMapper_Map( void* mapper );
	
#endif
