/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**
** This file may be distributed under the terms of the VPAC Public License
** as defined by VPAC of Australia and appearing in the file
** LICENSE.VPL included in the packaging of this file.
**
** This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
** WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
**
*/
/** \file
 **  Role:
 **	A OneToOne mapping between MaterialPointsSwarm and IntegrationPointsSwarm where the integration swarm uses a gauss
 **	layout so xi and weight calculation isn't needed. The materialPointsSwarm is assumed to use a BackgroundParticleLayout.
 **
 ** Assumptions:
 **
 ** Comments:
 **
 ** $Id: GaussMapper.h 189 2005-10-20 00:39:29Z RobertTurnbull $
 **
 **~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __PICellerator_MaterialPoints_GaussMapper_h__
#define __PICellerator_MaterialPoints_GaussMapper_h__

	/* Textual name of this class */
	extern const Type GaussMapper_Type;

	#define __GaussMapper \
		__OneToOneMapper 

	struct GaussMapper { __GaussMapper };

		

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define GAUSSMAPPER_DEFARGS \
                ONETOONEMAPPER_DEFARGS

	#define GAUSSMAPPER_PASSARGS \
                ONETOONEMAPPER_PASSARGS

	GaussMapper* _GaussMapper_New(  GAUSSMAPPER_DEFARGS  );

	void _GaussMapper_Init( void* mapper );

	void _GaussMapper_Delete( void* mapper );

	void _GaussMapper_Print( void* mapper, Stream* stream );

	#define GaussMapper_Copy( self ) \
		(GaussMapper*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define GaussMapper_DeepCopy( self ) \
		(GaussMapper*) Stg_Class_Copy( self, NULL, True, NULL, NULL )

	void* _GaussMapper_Copy( void* mapper, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	void* _GaussMapper_DefaultNew( Name name );

	void _GaussMapper_AssignFromXML( void* shape, Stg_ComponentFactory* cf, void* data );

	void _GaussMapper_Build( void* mapper, void* data );

	void _GaussMapper_Initialise( void* mapper, void* data );

	void _GaussMapper_Execute( void* mapper, void* data );

	void _GaussMapper_Destroy( void* mapper, void* data );
	
	void _GaussMapper_Map( void* mapper );
	
#endif

