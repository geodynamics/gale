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

#include <mpi.h>
/* EP_APPLICATIONS_FINALISE defined in StGermain.h */
#include <StGermain/StGermain.h>
/* Must include StgFEM library, but WON'T call Stg_FEM_Init in this plugin. */
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include "Application.h"

#include <stdio.h>
/* for strcmp */
#include <string.h>

/* All StGermain Application Plugins should use this finalise. */
/* char EP_APPLICATIONS_FINALISE[] = "StGermain_EP_ApplicationsFinalise"; */
const Type PICellerator_Application_Type = "PICellerator_Application";

void _PICellerator_Application_Construct( void* component, Stg_ComponentFactory* cf, void* data ) 
 {
	AbstractContext* currAbstractContext;
	PICelleratorContext* context = NULL;
	int tmp = 0;

        AbstractContext* prevContext;
        EntryPoint* applicationsFinalise_EP;

/* 	Get the existing abstract context, as defined by StGermain. */
/* 	Get it up here, so the communicator can be used for messages. */
	prevContext = (AbstractContext*)Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data ); 

	/* Ensures copyright info always come first in output */
	MPI_Barrier( prevContext->communicator ); 

        /* PICellerator's init message */
        tmp = Stream_GetPrintingRank( Journal_Register( InfoStream_Type, "Context" ) );
        Stream_SetPrintingRank( Journal_Register( InfoStream_Type, "Context" ), 0 );
        Journal_Printf( /* DO NOT CHANGE OR REMOVE */
                Journal_Register( InfoStream_Type, "Context" ),
                "PICellerator revision %s. "
		"Copyright (C) 2005 VPAC & Monash Cluster Computing.\n", 
		    VERSION );
        Stream_Flush( Journal_Register( InfoStream_Type, "Context" ) );
        Stream_SetPrintingRank( Journal_Register( InfoStream_Type, "Context" ), tmp );

	/* Ensures copyright info always come first in output */
        MPI_Barrier( prevContext->communicator ); 

/* 	Only need to initialise a new context, and copy all relevant registers over IF this is the first application */
/* 	plugin to be constructed. */
/* 	The first application plugin to be constructed is Guaranteed to have the 'largest' context. */
/* 				(ie is an inherited child of ALL other application plugins about to be loaded) */
	if( prevContext->type == AbstractContext_Type )
	{
/* 		Set the existing abstract context. */
		currAbstractContext = prevContext;
/* 		Create a new, empty PICelleratorContext. */
		context = PICelleratorContext_New( "context", 0, 0, currAbstractContext->communicator, cf->rootDict );
	
		context->dictionary = cf->rootDict;
	
/* 		Initialise Abstract parts of PICelleratorContext */
		_AbstractContext_Init((AbstractContext*)context, 0, 0, MPI_COMM_WORLD );
/* 		Initialise Discretisation parts of PICelleratorContext */
		_DiscretisationContext_Init((DiscretisationContext*)context );
/* 		Initialise FiniteElement parts of PICelleratorContext */
		_FiniteElementContext_Init( (FiniteElementContext*)context );
/* 		Initialise PIcellerator parts of PICelleratorContext */
		_PICelleratorContext_Init( context );
	
/* 	       	Need to get the old CF from currAbstractContext, and use that in my new context. */
/* 	        Now I CANNOT delete this currAbstractContext or I'll lose the CF. :( */
	        context->CF = currAbstractContext->CF;

/* 	        Need to get the LCRegister componentList, and replace the existing (abstract) context  */
/* 		with the newly created (PICellerator) context!!! */
	        Stg_ObjectList_Replace( context->CF->LCRegister->componentList,
                                	((Stg_Component*) currAbstractContext)->name, 
					KEEP, 
					(Stg_Component*) context);
	
/* 		Recreate the registerRegister link in CF. */
		context->CF->registerRegister = context->register_Register;

/* 		Create the EntryPoint for all application plugins' finalise functions to hook into. */
	        applicationsFinalise_EP = EntryPoint_New( EP_APPLICATIONS_FINALISE, EntryPoint_VoidPtr_CastType );
	        EntryPoint_Register_Add(context->entryPoint_Register, (void*)applicationsFinalise_EP);

	} /* close of if(context->type == AbstractContext_Type) */
	else /* prevContext was NOT an abstract context -> that is does NOT need to be replaced */
		context = (PICelleratorContext*) prevContext;


        EntryPoint_Append( Context_GetEntryPoint( context, EP_APPLICATIONS_FINALISE ),
                           "PICellerator App Finalise",
                           PICellerator_Application_Finalise,
                           "PICellerator_App_Construct" );
}

void PICellerator_Application_Finalise() 
{
	Stream* finaliseStream = 
			Journal_Register(Info_Type, "PICellerator_Application");
	Journal_Printf( finaliseStream, 
		      "Finalised: Particle-In-Cellerator (FEM/PIC Framework).\n");

	PICellerator_Finalise();
}

void* _PICellerator_Application_DefaultNew( Name name ) {
	return Codelet_New(
			PICellerator_Application_Type,
			_PICellerator_Application_DefaultNew,
			_PICellerator_Application_Construct,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
}

Index PICellerator_Application_Register( PluginsManager* pluginsManager ) 
{
/* 	Initialise the PICellerator context. */
        PICellerator_Init( NULL, NULL);
	#ifdef HAVE_PYTHON
	Py_Initialize();
	#endif

	return PluginsManager_Submit( pluginsManager, PICellerator_Application_Type, "0", _PICellerator_Application_DefaultNew );
}

