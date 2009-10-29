/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**      Melbourne, 3053, Australia.
**
** Primary Contributing Organisations:
**      Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**      Australian Computational Earth Systems Simulator - http://www.access.edu.au
**      Monash Cluster Computing - http://www.mcc.monash.edu.au
**      Computational Infrastructure for Geodynamics - http://www.geodynamics.org
**
** Contributors:
**      Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**      Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**      Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**      David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**      Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**      Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**      Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**      Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**      Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**      Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
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
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <string.h>
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include "FeVariableList.h"

#ifndef MASTER
        #define MASTER 0
#endif

const Type StgFEM_FeVariableList_Type = "StgFEM_FeVariableList";

void _StgFEM_FeVariableList_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
        StgFEM_FeVariableList*        self         = (StgFEM_FeVariableList*) component;
        AbstractContext*                context;
        Dictionary*                     dictionary;
        Stream*                         stream;
        Name                            fevariableListFilename;
        Bool                            fileOpened;
        Bool                            PrintToFile;

        context     = (AbstractContext*)Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data ); 
        dictionary  = context->dictionary;
        
        /* Append to end of build entry point */
        ContextEP_Append( context,  AbstractContext_EP_Build, StgFEM_FeVariableList_PrintVariables );  
        
        /* Create Stream */
        stream = self->stream = Journal_Register( InfoStream_Type, "FeVariableList" );
        
        /* Set auto flush on stream */
        Stream_SetAutoFlush( stream, True );    

        /* Print to screen or to file? */
        PrintToFile = Dictionary_GetBool_WithDefault( dictionary, "FeVariableListPrintToFile", True );
        if(PrintToFile){
                /* Get name of fevariable list file */
                fevariableListFilename = Dictionary_GetString_WithDefault( dictionary, "FeVariableListFilename", "FeVariables.list" );
                /* Open New File */
                if ( context->rank == MASTER ) {
                        Stream* errorStream = Journal_Register( Error_Type, CURR_MODULE_NAME );
                        fileOpened          = Stream_RedirectFile_WithPrependedPath( stream, context->outputPath, fevariableListFilename );
                        Journal_Firewall( fileOpened, errorStream, 
                                        "Could not open file %s/%s. Possibly directory %s does not exist or is not writable.\n"
                                        "Check 'outputPath' in input file.\n", context->outputPath, fevariableListFilename, context->outputPath );
                }
                /* Set it so only master processor can print to stream */
                Stream_SetPrintingRank( stream, MASTER );
        }       
}

void* _StgFEM_FeVariableList_DefaultNew( Name name ) {
        return _Codelet_New(
                        sizeof(StgFEM_FeVariableList),
                        StgFEM_FeVariableList_Type,
                        _Codelet_Delete,
                        _Codelet_Print,
                        _Codelet_Copy,
                        _StgFEM_FeVariableList_DefaultNew,
                        _StgFEM_FeVariableList_AssignFromXML,
                        _Codelet_Build,
                        _Codelet_Initialise,
                        _Codelet_Execute,
                        _Codelet_Destroy,
                        name );
}

Index StgFEM_FeVariableList_Register( PluginsManager* pluginsManager ) {
        return PluginsManager_Submit( pluginsManager, StgFEM_FeVariableList_Type, "0", _StgFEM_FeVariableList_DefaultNew );
}

void StgFEM_FeVariableList_PrintVariables( void* _context ){
        DomainContext*          context         = (DomainContext*) _context;
        Stream*                 stream;
        FieldVariable_Register* fV_Register;
        Index                   variablecount;
        Index                   countindex;
        Index                   columnWidth1     = 70;
        Index                   columnWidth2     = 30;
        
        StgFEM_FeVariableList* self   = (StgFEM_FeVariableList*)LiveComponentRegister_Get( 
                                                                        context->CF->LCRegister, 
                                                                        StgFEM_FeVariableList_Type );
        stream                          = self->stream;

        /* Get FeVariable Register*/
        fV_Register                     = context->fieldVariable_Register; 
        variablecount                   = (Index) fV_Register->objects->count;
        
        /* Print header material */
        Journal_Printf( stream, "\n");
        Journal_PrintString_WithLength( stream, "FeVariable", columnWidth1 );
        Journal_PrintString_WithLength( stream, "FeVariableType", columnWidth2 );
        Journal_Printf( stream, "\n");
        Journal_PrintString_WithLength( stream, "------------------------", columnWidth1 );
        Journal_PrintString_WithLength( stream, "------------------------", columnWidth2 );
        Journal_Printf( stream, "\n");  
        
        /* Print Variables */
        for(countindex = 1; countindex <= variablecount; ++countindex){
                Journal_PrintString_WithLength( stream, fV_Register->objects->data[ countindex - 1 ]->name, columnWidth1 );
                Journal_PrintString_WithLength( stream, fV_Register->objects->data[ countindex - 1 ]->type, columnWidth2 );
                Journal_Printf( stream, "\n");
        }
        Journal_Printf( stream, "\n");
}
