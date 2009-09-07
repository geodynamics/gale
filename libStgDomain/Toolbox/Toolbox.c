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

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "Toolbox.h"



const Type StgDomain_Toolbox_Type = "StgDomain_Toolbox";

void _StgDomain_Toolbox_Construct( void* component, Stg_ComponentFactory* cf, void* ptrToContext )  {
	DomainContext*			context;
	AbstractContext**		ptrToSelf;
	AbstractContext*		curContext;
 
	/* Get the existing context, as defined by StGermain, or an alternative toolbox. */
	ptrToSelf = (AbstractContext**)ptrToContext;
	curContext = (AbstractContext*)Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, ptrToSelf ); 


	/* Only need to initialise a new context, and copy all relevant registers over IF this is the first toolbox */
        /* to be constructed. */
        /* The first toolbox to be constructed is Guaranteed to have the 'largest' context. */
        /*                       (ie is an inherited child of ALL other toolbox about to be loaded) */
        if( curContext->type == AbstractContext_Type ) {
		/* Create a new, empty FiniteElementContext. */
		context = DomainContext_New( curContext->name, 0, 0, curContext->communicator, curContext->dictionary );
		Memory_CountInc( context );

                context->dictionary = curContext->dictionary;
                /* Need to get the old CF from curContext, and use that in my new context. */
                /* Now I CANNOT delete this curContext or I'll lose the CF. :( */
                context->CF = curContext->CF;

                /* Need to get the LCRegister componentList, and replace the existing (abstract) context */
                /* with the newly created (PICellerator) context!!! */
                Stg_ObjectList_Replace( 
			context->CF->LCRegister->componentList,
			((Stg_Component*) curContext)->name,
			KEEP,
			(Stg_Component*) context);

                /* Recreate the registerRegister link in CF. */
                context->CF->registerRegister = context->register_Register;
        } /* close of if( prevContext->type == AbstractContext_Type) */
	else { /* curContext was NOT an abstract context -> that is does NOT need to be replaced */
		context = (DomainContext*)curContext;
	}
	
	*ptrToSelf = (AbstractContext*)context;
}


void* _StgDomain_Toolbox_DefaultNew( Name name ) {
	return Codelet_New(
			StgDomain_Toolbox_Type,
			_StgDomain_Toolbox_DefaultNew,
			_StgDomain_Toolbox_Construct,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
}

void StgDomain_Toolbox_Initialise( PluginsManager* pluginsManager, int* argc, char*** argv ) {
	StgDomain_Init( argc, argv );
}

void StgDomain_Toolbox_Finalise( PluginsManager* pluginsManager ) {
	StgDomain_Finalise();
	
	Journal_RPrintf( Journal_Register( Info_Type, StgDomain_Toolbox_Type ), "Finalised: StGermain Domain Toolbox.\n" );
}

Index StgDomain_Toolbox_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, StgDomain_Toolbox_Type, "0", _StgDomain_Toolbox_DefaultNew );
}

