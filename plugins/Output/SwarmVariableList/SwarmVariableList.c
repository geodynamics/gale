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
** $Id: $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <string.h>
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include "SwarmVariableList.h"

#ifndef MASTER
        #define MASTER 0
#endif

const Type StgFEM_SwarmVariableList_Type = "StgFEM_SwarmVariableList";

void _StgFEM_SwarmVariableList_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
        StgFEM_SwarmVariableList* self         = (StgFEM_SwarmVariableList*) component;
        AbstractContext*                context;
        Dictionary*                     dictionary;
        Stream*                         stream;
        Name                            swarmVariableListFilename;
        Bool                            fileOpened;
        Bool                            PrintToFile;
        Swarm_Register*                 swarmRegister;

        self->context       = context = (AbstractContext*)Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data );
        dictionary          = context->dictionary;
        self->swarmRegister = Swarm_Register_GetSwarm_Register();

        /** create stream **/
        stream = self->stream = Journal_Register( InfoStream_Type, "SwarmVariableList" );

        /** set auto flush on stream **/
        Stream_SetAutoFlush( stream, True );

        /** print to screen or to file? **/
        PrintToFile = Dictionary_GetBool_WithDefault( dictionary, "SwarmVariableListPrintToFile", True );
        if(PrintToFile){
                /** get name of SwarmVariable list file **/
                swarmVariableListFilename = Dictionary_GetString_WithDefault( dictionary, "SwarmVariableListFilename", "SwarmVariables.list" );
                /** open new file **/
                if ( context->rank == MASTER ) {
                        Stream* errorStream = Journal_Register( Error_Type, CURR_MODULE_NAME );
                        fileOpened          = Stream_RedirectFile_WithPrependedPath( stream, context->outputPath, swarmVariableListFilename );
                        Journal_Firewall( fileOpened, errorStream,
                                        "Could not open file %s/%s. Possibly directory %s does not exist or is not writable.\n"
                                        "Check 'outputPath' in input file.\n", context->outputPath, swarmVariableListFilename, context->outputPath );
                }
                /** set it so only master processor can print to stream **/
                Stream_SetPrintingRank( stream, MASTER );
        }
}

void* _StgFEM_SwarmVariableList_DefaultNew( Name name ) {
   /* Variables set in this function */
   SizeT                                              _sizeOfSelf = sizeof(StgFEM_SwarmVariableList);
   Type                                                      type = StgFEM_SwarmVariableList_Type;
   Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
   Stg_Class_PrintFunction*                                _print = _Codelet_Print;
   Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
   Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _StgFEM_SwarmVariableList_DefaultNew;
   Stg_Component_ConstructFunction*                    _construct = _StgFEM_SwarmVariableList_AssignFromXML;
   Stg_Component_BuildFunction*                            _build = _Codelet_Build;
   Stg_Component_InitialiseFunction*                  _initialise = StgFEM_SwarmVariableList_PrintVariables;
   Stg_Component_ExecuteFunction*                        _execute = _Codelet_Execute;
   Stg_Component_DestroyFunction*                        _destroy = _Codelet_Destroy;

   /* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
   AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

        return _Codelet_New(  CODELET_PASSARGS  );
}

Index StgFEM_SwarmVariableList_Register( PluginsManager* pluginsManager ) {
        return PluginsManager_Submit( pluginsManager, StgFEM_SwarmVariableList_Type, "0", _StgFEM_SwarmVariableList_DefaultNew );
}

void StgFEM_SwarmVariableList_PrintVariables( void* _self, void* data ){
   StgFEM_SwarmVariableList* self    = (StgFEM_SwarmVariableList*)_self;
   AbstractContext*          context = self->context;
   Stream*                   stream  = self->stream;;
   Index                     swarmCount;
   Index                     variablecount;
   Index                     countindex;
   Index                     swarmcountindex;
   Index                     columnWidth     = 40;
   Swarm*                    currentSwarm;
   
   /** print header material **/
   Journal_Printf( stream, "\n");
   Journal_PrintString_WithLength( stream, "SwarmVariable", columnWidth );
   Journal_PrintString_WithLength( stream, "Swarm", columnWidth );
   Journal_Printf( stream, "\n");
   Journal_PrintString_WithLength( stream, "------------------------", columnWidth );
   Journal_PrintString_WithLength( stream, "------------------------", columnWidth );
   Journal_Printf( stream, "\n");
   
   /** get total number of different swarms **/
   swarmCount  = self->swarmRegister->swarmList->count;
   
   /** print swarm variables **/
   for(swarmcountindex = 0; swarmcountindex < swarmCount; ++swarmcountindex){
          currentSwarm = (Swarm*)Stg_ComponentFactory_ConstructByName( context->CF, self->swarmRegister->swarmList->data[swarmcountindex]->name, Swarm, True, NULL );
          variablecount = currentSwarm->swarmVariable_Register->objects->count;
          for(countindex = 0; countindex < variablecount; ++countindex){
                  Journal_PrintString_WithLength( stream, currentSwarm->swarmVariable_Register->objects->data[ countindex ]->name, columnWidth );
                  Journal_PrintString_WithLength( stream, self->swarmRegister->swarmList->data[swarmcountindex]->name, columnWidth );
                  Journal_Printf( stream, "\n");
          }
   }
   Journal_Printf( stream, "\n");
}


