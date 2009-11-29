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
**	Instantiates the ParticleLayout abstract class to a manually distributed particle layout.
**
** Assumptions:
**	Cell is a right-angled cuboid.
**
** Comments:
**
** $Id: FileParticleLayout.h 3851 2006-10-12 08:57:22Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Domain_Swarm_FileParticleLayout_h__
#define __Domain_Swarm_FileParticleLayout_h__
	

	/* Textual name of this class */
	extern const Type FileParticleLayout_Type;
	
	/* FileParticleLayout information */
#ifdef READ_HDF5
	#define __FileParticleLayout \
		__GlobalParticleLayout \
		\
		Name                                             filename;    \
		Stream*                                          errorStream; \
		hid_t** fileData; \
		hid_t** fileSpace; \
		Index* lastParticleIndex; \
		hsize_t start[2]; \
		hsize_t count[2]; \
		/** number of files previous checkpoint stored across */ \
		Index                           checkpointfiles;		
#else
	#define __FileParticleLayout \
		__GlobalParticleLayout \
		\
		Name       filename;    \
		FILE*      file;        \
		Stream*    errorStream; \
		/** number of files previous checkpoint stored across */ \
		Index      checkpointfiles;
#endif

	struct FileParticleLayout { __FileParticleLayout };
	
	/* Create a new FileParticleLayout and initialise */
   FileParticleLayout* FileParticleLayout_New( Name name,
      AbstractContext* context, 
      CoordSystem      coordSystem,
      Bool             weightsInitialisedAtStartup,
      unsigned int     totalInitialParticles, 
      double           averageInitialParticlesPerCell, 
      Name             filename, 
      Index            checkpointfiles );

   /* Creation implementation / Virtual constructor */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define FILEPARTICLELAYOUT_DEFARGS \
                GLOBALPARTICLELAYOUT_DEFARGS, \
                Name          filename, \
                Index  checkpointfiles

	#define FILEPARTICLELAYOUT_PASSARGS \
                GLOBALPARTICLELAYOUT_PASSARGS, \
	        filename,        \
	        checkpointfiles

   FileParticleLayout* _FileParticleLayout_New(  FILEPARTICLELAYOUT_DEFARGS  );
	
	void _FileParticleLayout_Init( void* particleLayout, Name filename, Index checkpointfiles );
	
	/* 'Stg_Class' Stuff */
	void _FileParticleLayout_Delete( void* particleLayout );
	void _FileParticleLayout_Print( void* particleLayout, Stream* stream );
	#define FileParticleLayout_Copy( self ) \
		(FileParticleLayout*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define FileParticleLayout_DeepCopy( self ) \
		(FileParticleLayout*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _FileParticleLayout_Copy( void* particleLayout, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	/* 'Stg_Component' Stuff */
	void* _FileParticleLayout_DefaultNew( Name name ) ;
	void _FileParticleLayout_AssignFromXML( void* particleLayout, Stg_ComponentFactory *cf, void* data );
	void _FileParticleLayout_Build( void* particleLayout, void* data );
	void _FileParticleLayout_Initialise( void* particleLayout, void* data );
	void _FileParticleLayout_Execute( void* particleLayout, void* data );
	void _FileParticleLayout_Destroy( void* particleLayout, void* data );
	
	void _FileParticleLayout_SetInitialCounts( void* particleLayout, void* _swarm ) ;
	void _FileParticleLayout_InitialiseParticles( void* particleLayout, void* _swarm ) ;
	void _FileParticleLayout_InitialiseParticle( void* particleLayout, void* swarm, Particle_Index newParticle_I, void* particle);
	/* small routine to find out number of files fileParticleLayout is stored across, which maybe have been stored in the timeInfo checkpoint file */ 
	Index _FileParticleLayout_GetFileCountFromTimeInfoFile( void* context );

#endif /* __Domain_Swarm_FileParticleLayout_h__ */

