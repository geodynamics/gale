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
** $Id: FrequentOutput.c 1219 2008-09-04 23:09:02Z JohnMansour $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <string.h>
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include "FrequentOutput.h"

#ifndef MASTER
	#define MASTER 0
#endif

const Type StgFEM_FrequentOutput_Type = "StgFEM_FrequentOutput";

void _StgFEM_FrequentOutput_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
	StgFEM_FrequentOutput*            self         = (StgFEM_FrequentOutput*) component;
	AbstractContext*                   context;
	Dictionary*                        dictionary;
	Stream*                            stream;
	Name                               frequentOutputFilename;
	Bool                               fileOpened;
	Stream*                            errorStream  = Journal_Register( Error_Type, CURR_MODULE_NAME );
	Dictionary*			   pluginDict	= Codelet_GetPluginDictionary( self, cf->rootDict );

	context = Stg_ComponentFactory_ConstructByName( cf, Dictionary_GetString( pluginDict, "Context" ), AbstractContext, True, data );
	self->context = context;
	dictionary = context->dictionary;
	
	ContextEP_Append( context, AbstractContext_EP_Initialise, StgFEM_FrequentOutput_PrintNewLine );
	ContextEP_Prepend_AlwaysFirst( context, AbstractContext_EP_FrequentOutput, StgFEM_FrequentOutput_PrintTime );
	ContextEP_Append_AlwaysLast(   context, AbstractContext_EP_FrequentOutput, StgFEM_FrequentOutput_PrintNewLine );
	
	/* Create Stream */
	stream = self->stream = Journal_Register( InfoStream_Type, "FrequentOutputFile" );

	/* Set auto flush on stream */
	Stream_SetAutoFlush( stream, True );

	/* Get name of frequent output file */
	frequentOutputFilename = Dictionary_GetString_WithDefault( dictionary, "FrequentOutputFilename", "FrequentOutput.dat" );

	/* Open File */
	if ( context->rank == MASTER ) {
		if ( (context->loadFromCheckPoint == False) && (Dictionary_GetBool_WithDefault( context->dictionary, "visualOnly", False ) == False ) ) {
			/* Always overwrite the file if starting a new run */
			fileOpened = Stream_RedirectFile_WithPrependedPath( stream, context->outputPath, frequentOutputFilename );
		}
		else {
			/* Just append to the file if doing a restart from checkpoint */
			fileOpened = Stream_AppendFile_WithPrependedPath( stream, context->outputPath, frequentOutputFilename );
		}
		Journal_Firewall( fileOpened, errorStream, 
				"Could not open file %s/%s. Possibly directory %s does not exist or is not writable.\n"
				"Check 'outputPath' in input file.\n", context->outputPath, frequentOutputFilename, context->outputPath );
	}
	
	/* Set it so only master processor can print to stream */
	Stream_SetPrintingRank( stream, MASTER );

	/* Read in values from dictionary */
	self->columnWidth   = Dictionary_GetUnsignedInt_WithDefault( dictionary, "FrequentOutputColumnWidth",  12 );
	self->decimalLength = Dictionary_GetUnsignedInt_WithDefault( dictionary, "FrequentOutputDecimalLength", 6 );
	self->borderString  = Dictionary_GetString_WithDefault( dictionary, "FrequentOutputBorderString", "    " );

	StgFEM_FrequentOutput_PrintHeader( context );
}

void* _StgFEM_FrequentOutput_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(StgFEM_FrequentOutput);
	Type                                                      type = StgFEM_FrequentOutput_Type;
	Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
	Stg_Class_PrintFunction*                                _print = _Codelet_Print;
	Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _StgFEM_FrequentOutput_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _StgFEM_FrequentOutput_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _Codelet_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _Codelet_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _Codelet_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _Codelet_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;

	return _Codelet_New(  CODELET_PASSARGS  );
}

Index StgFEM_FrequentOutput_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, StgFEM_FrequentOutput_Type, "0", _StgFEM_FrequentOutput_DefaultNew );
}

void StgFEM_FrequentOutput_PrintString( void* _context, char* string ) {
	AbstractContext*                   context = (AbstractContext*) _context;
	Stream*                            stream;

	StgFEM_FrequentOutput* self = (StgFEM_FrequentOutput*)LiveComponentRegister_Get( 
									context->CF->LCRegister, 
									StgFEM_FrequentOutput_Type );
	stream     = self->stream;

	/* Print some empty space at start */
	Journal_Printf( stream, self->borderString );

	/* Print String - Truncated if nessesary */
	Journal_PrintString_WithLength( stream, string, self->columnWidth );
}

void StgFEM_FrequentOutput_PrintDouble( void* _context, double value ) {
	AbstractContext*                   context = (AbstractContext*) _context;
	char*                              formatString;
	Stream*                            stream;

	StgFEM_FrequentOutput* self = (StgFEM_FrequentOutput*)LiveComponentRegister_Get( 
									context->CF->LCRegister,
									StgFEM_FrequentOutput_Type );

	stream     = self->stream;

	/* Create format String */
	Stg_asprintf( &formatString, "%%%d.%dg", self->columnWidth, self->decimalLength );

	Journal_Printf( self->stream, self->borderString );
	Journal_Printf( self->stream, formatString, value );

	Memory_Free( formatString );
}

void StgFEM_FrequentOutput_PrintHeader( void* _context ) {
	AbstractContext*                   context = (AbstractContext*) _context;
	char*                              firstBorderString;
	Stream*                            stream;

	StgFEM_FrequentOutput* self = (StgFEM_FrequentOutput*)LiveComponentRegister_Get( 
									context->CF->LCRegister, 
									StgFEM_FrequentOutput_Type );
	
	stream     = self->stream;

	/* Print First Boarder with '#' in the front */
	firstBorderString = StG_Strdup( self->borderString );
	firstBorderString[0] = '#';
	Journal_Printf( stream, firstBorderString );
	Memory_Free( firstBorderString );

	Journal_PrintString_WithLength( stream, "Timestep", self->columnWidth );

	StgFEM_FrequentOutput_PrintString( context, "Time" );
}

void StgFEM_FrequentOutput_PrintTime( void* _context ) {
	AbstractContext* context = (AbstractContext*) _context;

	StgFEM_FrequentOutput_PrintValue( context, context->timeStep );
	StgFEM_FrequentOutput_PrintValue( context, context->currentTime );
}

void StgFEM_FrequentOutput_PrintNewLine( void* _context ) {
	AbstractContext*                   context = (AbstractContext*) _context;
	Stream*                            stream;

	StgFEM_FrequentOutput* self = (StgFEM_FrequentOutput*)LiveComponentRegister_Get( 
									context->CF->LCRegister, 
									StgFEM_FrequentOutput_Type );
	
	stream     = self->stream;

	Journal_Printf( stream, "\n" );
}


