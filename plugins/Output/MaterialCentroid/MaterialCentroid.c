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
** $Id: MaterialCentroid.c 518 2007-10-11 08:07:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <assert.h>

const Type PICellerator_MaterialCentroid_Type = "PICellerator_MaterialCentroid";

void MaterialCentroid( PICelleratorContext* context ) {
	static IntegrationPointsSwarm*      swarm;
	static Material*            material;
	Coord                       centroid;
	double                      volume;
	static Bool                 beenHere              = False;
	Stream*                     stream                = Journal_Register( Info_Type, "MaterialCentroid" );

	if (!beenHere) {
		Name                 swarmName;
		Name                 materialName;
		Name                 filename;
      Bool                 fileOpened;
      Stream*              errorStream  = Journal_Register( Error_Type, CURR_MODULE_NAME );

		swarmName = Dictionary_GetString_WithDefault( context->dictionary, "MaterialCentroid_Swarm", "picIntegrationPoints" );
		
		/* Get Swarm */
		swarm = (IntegrationPointsSwarm*) LiveComponentRegister_Get( context->CF->LCRegister, swarmName );
		assert( swarm );
	
		/* Get Material */
		materialName = Dictionary_GetString( context->dictionary, "CentroidMaterial" );
		material = Materials_Register_GetByName( swarm->materials_Register, materialName );
		assert( material );
		
		/* Set up stream */
		Stg_asprintf( &filename, "%sCentroid.dat", materialName );
      /* Open File */
      if ( context->rank == 0 ) {
         if ( context->loadFromCheckPoint == False ) {
            /* Always overwrite the file if starting a new run */
            fileOpened = Stream_RedirectFile_WithPrependedPath( stream, context->outputPath, filename );
         } else {
            /* Just append to the file if doing a restart from checkpoint */
            fileOpened = Stream_AppendFile_WithPrependedPath( stream, context->outputPath, filename );
         }
         Journal_Firewall( fileOpened, errorStream, 
               "Could not open file %s/%s. Possibly directory %s does not exist or is not writable.\n"
               "Check 'outputPath' in input file.\n", context->outputPath, filename, context->outputPath );
      }
		Memory_Free( filename );
		Stream_SetAutoFlush( stream, True );

		/* Print header to stream */
		Journal_Printf( stream, 
				"#       Timestep            Time          Volume       CentroidX       CentroidY       CentroidZ\n" );

		beenHere = True;
	}
		
	volume = Material_Volume( material, swarm, centroid );

	Journal_Printf( stream, "    %12.6g    %12.6g    %12.6g    %12.6g    %12.6g    %12.6g\n",
			(double)context->timeStep, context->currentTime, volume, centroid[ I_AXIS ], centroid[ J_AXIS ], centroid[ K_AXIS ] );
}


void _PICellerator_MaterialCentroid_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {

	AbstractContext* context;

	context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data ); 
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput, MaterialCentroid );
}


void* _PICellerator_MaterialCentroid_DefaultNew( Name name ) {
	return _Codelet_New(
			sizeof( Codelet ),
			PICellerator_MaterialCentroid_Type,
			_Codelet_Delete,
			_Codelet_Print,
			_Codelet_Copy,
			_PICellerator_MaterialCentroid_DefaultNew,
			_PICellerator_MaterialCentroid_AssignFromXML,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
}


Index PICellerator_MaterialCentroid_Register( PluginsManager* pluginsManager ) {
	Index result;

	result = PluginsManager_Submit( pluginsManager, PICellerator_MaterialCentroid_Type, "0",
		_PICellerator_MaterialCentroid_DefaultNew );

	return result;
}
