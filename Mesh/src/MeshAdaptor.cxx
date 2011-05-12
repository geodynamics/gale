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
** $Id: MeshAdaptor.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>
#include <StGermain/StGermain.h>

#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>

#include "types.h"
#include "shortcuts.h"
#include "Decomp.h"
#include "MeshTopology.h"
#include "MeshClass.h"
#include "MeshGenerator.h"
#include "MeshAdaptor.h"


/* Textual name of this class */
const Type MeshAdaptor_Type = "MeshAdaptor";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

MeshAdaptor* _MeshAdaptor_New(  MESHADAPTOR_DEFARGS  ) {
	MeshAdaptor* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(MeshAdaptor) );
	self = (MeshAdaptor*)_MeshGenerator_New(  MESHGENERATOR_PASSARGS  );

	/* Virtual info */

	return self;
}

void _MeshAdaptor_Init( MeshAdaptor* self ) {
	self->generator = NULL;
	self->srcMesh = NULL;
}

/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _MeshAdaptor_Delete( void* adaptor ) {
	MeshAdaptor*	self = (MeshAdaptor*)adaptor;

	/* Delete the parent. */
	_MeshGenerator_Delete( self );
}

void _MeshAdaptor_Print( void* adaptor, Stream* stream ) {
	MeshAdaptor*	self = (MeshAdaptor*)adaptor;
	
	/* Set the Journal for printing informations */
	Stream* adaptorStream;
	adaptorStream = Journal_Register( InfoStream_Type, (Name)"MeshAdaptorStream"  );

	/* Print parent */
	Journal_Printf( stream, "MeshAdaptor (ptr): (%p)\n", self );
	_Stg_Component_Print( self, stream );
}

void _MeshAdaptor_AssignFromXML( void* adaptor, Stg_ComponentFactory* cf, void* data ) {
	MeshAdaptor*	self = (MeshAdaptor*)adaptor;

	_MeshGenerator_AssignFromXML( self, cf, data );

	/* There could be either a generator or a mesh to use as a template.  Prefer the mesh. */
	self->srcMesh = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"sourceMesh", Mesh, False, data );
	if( !self->srcMesh  ) {
		/* Read the source generator. */
		self->generator = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"sourceGenerator", MeshGenerator, True, data  );
	}
}

void _MeshAdaptor_Build( void* _adaptor, void* data ) {
   MeshAdaptor*      self = (MeshAdaptor*)_adaptor;
   if(self->generator)
     Stg_Component_Build( self->generator, data, False );
   if(self->srcMesh)
     Stg_Component_Build( self->srcMesh, data, False );
   _MeshGenerator_Build( self, data );
}

void _MeshAdaptor_Initialise( void* _adaptor, void* data ) {
   MeshAdaptor*      self = (MeshAdaptor*)_adaptor;
   if(self->generator)
     Stg_Component_Initialise( self->generator, data, False );
   if(self->srcMesh)
     Stg_Component_Initialise( self->srcMesh, data, False );
   _MeshGenerator_Initialise( self, data );
   
}

void _MeshAdaptor_Execute( void* adaptor, void* data ) {
    _MeshGenerator_Execute( adaptor, data );
}

void _MeshAdaptor_Destroy( void* _adaptor, void* data ) {
   MeshAdaptor*      self = (MeshAdaptor*)_adaptor;
   if(self->generator)
     Stg_Component_Destroy( self->generator, data, False );
   if(self->srcMesh)
     Stg_Component_Destroy( self->srcMesh, data, False );
   _MeshGenerator_Destroy( self, data );
}

/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void MeshAdaptor_SetGenerator( void* adaptor, void* generator ) {
	MeshAdaptor*	self = (MeshAdaptor*)adaptor;

	self->generator = (MeshGenerator*)generator;
	if( self->generator )
		self->srcMesh = NULL;
}

void MeshAdaptor_SetSourceMesh( void* adaptor, void* mesh ) {
	MeshAdaptor*	self = (MeshAdaptor*)adaptor;

	self->srcMesh = (Mesh*)mesh;
	if( self->srcMesh )
		self->generator = NULL;
}

/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/


