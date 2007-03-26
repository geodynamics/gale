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
** $Id: testParticleCoords.c 376 2006-10-18 06:58:41Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/Voronoi/Voronoi.h>
#include <PICellerator/PopulationControl/PopulationControl.h>

#include <math.h>
#include <string.h>
#include <assert.h>

const Type TestParticleCoords_Type = "TestParticleCoords";

void PICellerator_testParticleCoords( DiscretisationContext* context ) {
	Swarm*     swarm;
	Stream*             stream  = Journal_Register( Info_Type, CURR_MODULE_NAME );
	Name                filename;

	if ( context->timeStep == 0 )
		return;

	Stg_asprintf( &filename, "output%u.dat", context->rank );

	Stream_RedirectFile_WithPrependedPath( stream, context->outputPath, filename );

	Memory_Free( filename );

	swarm = (Swarm*)            LiveComponentRegister_Get( context->CF->LCRegister, "picIntegrationPoints" );
	assert(swarm);

	Swarm_PrintParticleCoords( swarm, stream );
}
	
void _testParticleCoords_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	DiscretisationContext* context;
	context = (DiscretisationContext*)Stg_ComponentFactory_ConstructByName( cf, "context", DiscretisationContext, True, data ); 
	
	ContextEP_Append( context, AbstractContext_EP_Dump, PICellerator_testParticleCoords );
}


void* _testParticleCoords_DefaultNew( Name name ) {
	return _Codelet_New(
			sizeof( Codelet ),
			TestParticleCoords_Type,
			_Codelet_Delete,
			_Codelet_Print,
			_Codelet_Copy,
			_testParticleCoords_DefaultNew,
			_testParticleCoords_Construct,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
}


Index testParticleCoords_Register( PluginsManager* pluginsManager ) {
	Index result;

	result = PluginsManager_Submit( pluginsManager, TestParticleCoords_Type, "0",
		_testParticleCoords_DefaultNew );

	return result;
}
