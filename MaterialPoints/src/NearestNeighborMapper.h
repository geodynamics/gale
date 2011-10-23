/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
** Copyright (c) 2005-2006, Monash Cluster Computing, Building 28, Monash University Clayton Campus,
**	Victoria, 3800, Australia
** Copyright (c) 2011, California Institute of Technology
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
**      Walter Landry, CIG (walter@geodynamics.org)
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
**     An IntegrationPointMapper which maps the nearest neighbors of a
**     IntegrationPointsSwarm to an IntegrationPointsSwarm.
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __PICellerator_MaterialPoints_NearestNeighborMapper_h__
#define __PICellerator_MaterialPoints_NearestNeighborMapper_h__

	extern const Type NearestNeighborMapper_Type;

	/* NearestNeighborMapper information */
	#define __NearestNeighborMapper \
		__IntegrationPointMapper \
		\
		Stream*				errorStream; \
		IntegrationPointsSwarm*		swarm; \

	struct NearestNeighborMapper { __NearestNeighborMapper };

	#ifndef ZERO
	#define ZERO 0
	#endif

	#define NEARESTNEIGHBORMAPPER_DEFARGS \
                INTEGRATIONPOINTMAPPER_DEFARGS, \
                IntegrationPointsSwarm *_swarm

	#define NEARESTNEIGHBORMAPPER_PASSARGS \
                INTEGRATIONPOINTMAPPER_PASSARGS, \
                _swarm

NearestNeighborMapper* _NearestNeighborMapper_New( NEARESTNEIGHBORMAPPER_DEFARGS );

void _NearestNeighborMapper_Init( void* mapper, IntegrationPointsSwarm* swarm );

	void _NearestNeighborMapper_Delete( void* mapper );
	void _NearestNeighborMapper_Print( void* mapper, Stream* stream );
	#define NearestNeighborMapper_Copy( self ) \
		(NearestNeighborMapper*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define NearestNeighborMapper_DeepCopy( self ) \
		(NearestNeighborMapper*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _NearestNeighborMapper_Copy( const void* mapper, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	void* _NearestNeighborMapper_DefaultNew( Name name );
	void _NearestNeighborMapper_AssignFromXML( void* shape, Stg_ComponentFactory* cf, void* data );
	void _NearestNeighborMapper_Build( void* mapper, void* data ) ;
	void _NearestNeighborMapper_Initialise( void* mapper, void* data );
	void _NearestNeighborMapper_Execute( void* mapper, void* data );
	void _NearestNeighborMapper_Destroy( void* mapper, void* data );

        void _NearestNeighborMapper_Map(void* mapper);

        MaterialPointsSwarm** _NearestNeighborMapper_GetMaterialPointsSwarms
        ( void* mapper, Index* count );

        Material_Index _NearestNeighborMapper_GetMaterialIndexOn
        ( void* mapper, void* point );

        void* _NearestNeighborMapper_GetExtensionOn
        (void* mapper, void* point, ExtensionInfo_Index extHandle );

        double _NearestNeighborMapper_GetDoubleFromExtension
        (void* mapper,
         void* intPoint,
         ExtensionInfo_Index extHandle,
         int offs);

        double _NearestNeighborMapper_GetDoubleFromMaterial
        (void* mapper,
         void* intPoint,
         ExtensionInfo_Index extHandle,
         int offs);

        int NearestNeighbor_FindNeighbor(void* mapper,
                                         const Element_LocalIndex &lElement_I,
                                         const int &cell_I,
                                         double *xi, const int &dim);

        void NearestNeighbor_Replace(IntegrationPointsSwarm **swarm,
                                     IntegrationPoint **particle,
                                     int *particle_index,
                                     const Element_LocalIndex &lElement_I,
                                     const int &dim);

        inline void NearestNeighbor_Replace(IntegrationPointsSwarm **swarm,
                                            IntegrationPoint **particle,
                                            const Element_LocalIndex &lElement_I,
                                            const int &dim)
        {
          int temp(0);
          NearestNeighbor_Replace(swarm,particle,&temp,lElement_I,dim);
        }
#endif
