/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**	Computational Infrastructure for Geodynamics - http://www.geodynamics.org
**
** Contributors:
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
**	Kent Humphries, Software Engineer, VPAC. (kenth@vpac.org)
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
** $Id: StandardConditionFunctions.c 532 2006-04-04 00:21:59Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#ifdef HAVE_PYTHON
#include <Python.h>
#endif
#ifdef HAVE_SDL
#include <SDL/SDL.h>
#endif

#include <mpi.h>
/*EP_APPLICATIONS_FINALISE defined in StGermain.h */
#include <StGermain/StGermain.h>
/*Must include StgFEM library, but WON'T call Stg_FEM_Init in this plugin. */
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include "Application.h"

#include <stdio.h>
/*for strcmp */
#include <string.h>

const Type Underworld_Application_Type = "Underworld_Application";

void _Underworld_Application_Construct( void* component, Stg_ComponentFactory* cf, void* data ) 
 {
	AbstractContext* currAbstractContext;
	UnderworldContext* context = NULL;

        AbstractContext* prevContext;
        EntryPoint* applicationsFinalise_EP;
	#if 0
	EntryPoint* applicationsPostConstruct_EP;
	#endif

	/*Get the existing abstract context, as defined by StGermain. */
	/*Get it up here, so the communicator can be used for messages. */
	prevContext = (AbstractContext*)Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data ); 

	/*Only need to initialise a new context, and copy all relevant registers over IF this is the first application */
	/*plugin to be constructed. */
	/*The first application plugin to be constructed is Guaranteed to have the 'largest' context. */
	/*			(ie is an inherited child of ALL other application plugins about to be loaded) */
	if( prevContext->type == AbstractContext_Type )
	{
          /*Set the existing abstract context. */
		currAbstractContext = prevContext;
		/*Create a new, empty UnderworldContext. */
		context = UnderworldContext_New( "context", 0, 0, currAbstractContext->communicator, cf->rootDict );
	
		context->dictionary = cf->rootDict;
	
		/*Initialise Abstract parts of UnderworldContext */
		_AbstractContext_Init((AbstractContext*)context, 0, 0, MPI_COMM_WORLD );
		/*Initialise Discretisation parts of UnderworldContext */
		_DiscretisationContext_Init((DiscretisationContext*)context );
		/*Initialise FiniteElement parts of UnderworldContext */
		_FiniteElementContext_Init( (FiniteElementContext*)context );
		/*Initialise PIcellerator parts of UnderworldContext */
		_PICelleratorContext_Init( (PICelleratorContext*)context );
		/*Initialise Underworld parts of UnderworldContext */
		_UnderworldContext_Init( context );
	
	       	/*Need to get the old CF from currAbstractContext, and use that in my new context. */
	        /*Now I CANNOT delete this currAbstractContext or I'll lose the CF. :( */
	        context->CF = currAbstractContext->CF;

	        /*Need to get the LCRegister componentList, and replace the existing (abstract) context  */
		/*with the newly created (Underworld) context!!! */
	        Stg_ObjectList_Replace( context->CF->LCRegister->componentList,
                                	((Stg_Component*) currAbstractContext)->name, 
					KEEP, 
					(Stg_Component*) context);
	
		/*Recreate the registerRegister link in CF. */
		context->CF->registerRegister = context->register_Register;

		/*Create the EntryPoint for all application plugins' finalise functions to hook into. */
	        applicationsFinalise_EP = EntryPoint_New( EP_APPLICATIONS_FINALISE, EntryPoint_VoidPtr_CastType );
	        EntryPoint_Register_Add(context->entryPoint_Register, (void*)applicationsFinalise_EP);
	} /*close of if(context->type == AbstractContext_Type) */
	else /*prevContext was NOT an abstract context -> that is does NOT need to be replaced */
		context = (UnderworldContext*) prevContext;

        EntryPoint_Append( Context_GetEntryPoint( context, EP_APPLICATIONS_FINALISE ),
                           "Underworld App Finalise",
                           Underworld_Application_Finalise,
                           "Underworld_App_Construct" );
}

void Underworld_Application_Finalise() 
{
	Stream* finaliseStream = 
			Journal_Register(Info_Type, "Underworld_Application");
	Journal_Printf( finaliseStream, 
		      "Finalised: Underworld (Geodynamics Framework).\n");

	Underworld_Finalise();
}

void* _Underworld_Application_DefaultNew( Name name ) {
	return Codelet_New(
			Underworld_Application_Type,
			_Underworld_Application_DefaultNew,
			_Underworld_Application_Construct,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
}

Index Underworld_Application_Register( PluginsManager* pluginsManager ) 
{
  /*Initialise the Underworld context. */
        Underworld_Init( NULL, NULL);
	#ifdef HAVE_PYTHON
	Py_Initialize();
	#endif

	return PluginsManager_Submit( pluginsManager, Underworld_Application_Type, "0", _Underworld_Application_DefaultNew );
}

