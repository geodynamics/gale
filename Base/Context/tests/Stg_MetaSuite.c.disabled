/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**   Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**   Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**   Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**   Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**   Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**   Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
** $Id: testJournal-Dictionary.c 2745 2005-03-05 08:12:18Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"
#include "StGermain/Base/Container/Container.h"
#include "StGermain/Base/Automation/Automation.h"
#include "Stg_MetaSuite.h"


typedef struct {
   XML_IO_Handler*         io;
   Dictionary*             metaDict;
   LiveComponentRegister*  lcRegister;
} Stg_MetaSuiteData;


void Stg_MetaSuite_Setup( Stg_MetaSuiteData* data ) {
   data->io = XML_IO_Handler_New();
   data->metaDict = Dictionary_New();
   data->lcRegister = LiveComponentRegister_New( );
}


void Stg_MetaSuite_Teardown( Stg_MetaSuiteData* data ) {
   Stg_Class_Delete( data->io );
   Stg_Class_Delete( data->metaDict );
   Stg_Class_Delete( data->lcRegister );
}


/* Use this function if you need to regenerate the testing .xml file */
void Stg_MetaSuite_RegenerateComponentVC_MetaFile( Stg_MetaSuiteData* data ) {
   Stg_Class_Delete( data->metaDict );
   data->metaDict = Stg_ComponentRegister_GetMetadata( Stg_ComponentRegister_Get_ComponentRegister(), CompositeVC_Type, "0" );
   IO_Handler_WriteAllToFile( data->io, "ComponentVCMetaDict.xml", data->metaDict );
}


void Stg_MetaSuite_TestGetFunctions( Stg_MetaSuiteData* data ) {
   char     metaFilename[PCU_PATH_MAX];
   Index    paramCount=0;
   Index    assocCount=0;

   pcu_filename_input( "ComponentVCMetaDict.xml", metaFilename );
   IO_Handler_ReadAllFromFile( data->io, metaFilename, data->metaDict );

   pcu_check_streq( Stg_Meta_GetType( data->metaDict ), "CompositeVC" );
	pcu_check_streq( Stg_Meta_GetCreator( data->metaDict ), "VPAC" );
	pcu_check_streq( Stg_Meta_GetPublisher( data->metaDict ), "VPAC" );
	pcu_check_streq( Stg_Meta_GetRights( data->metaDict ), "The Gnu Lesser General Public License v2.1 - http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html" );
	pcu_check_streq( Stg_Meta_GetSource( data->metaDict ), "./StGermain/Base/Automation/src/" );
	pcu_check_streq( Stg_Meta_GetSubject( data->metaDict ), "StGermain" );
	pcu_check_streq( Stg_Meta_GetDescription( data->metaDict ), "CompositeVC is used as a container to combine multiple variable conditions." );
	pcu_check_streq( Stg_Meta_GetExampleDocumentation( data->metaDict ), "" );
	pcu_check_streq( Stg_Meta_GetExampleCode( data->metaDict ), "Refer to ./StgFEM/Apps/StokesMomentumUzawa/lidDrivenBCs.xml" );
	pcu_check_streq( Stg_Meta_GetInherits( data->metaDict ), "VariableCondition" );

	pcu_check_streq( Stg_Meta_GetReference( data->metaDict ), "" );
	pcu_check_streq( Stg_Meta_GetEquation( data->metaDict ), "" );

	paramCount = Stg_Meta_GetParameterCount( data->metaDict );
   pcu_check_true( paramCount == 2 );
	pcu_check_streq( Stg_Meta_GetParameterName( data->metaDict, 0 ), "vcName" );
	pcu_check_streq( Stg_Meta_GetParameterType( data->metaDict, 0 ), "xsd:string" );
	pcu_check_streq( Stg_Meta_GetParameterDefault( data->metaDict, 0 ), "self->name" );
	pcu_check_streq( Stg_Meta_GetParameterDocumentation( data->metaDict, 0 ), "Deprecated; should not be used." );
	pcu_check_streq( Stg_Meta_GetParameterName( data->metaDict, 1 ), "vcList" );
	pcu_check_streq( Stg_Meta_GetParameterType( data->metaDict, 1 ), "stg:list" );
	pcu_check_streq( Stg_Meta_GetParameterDefault( data->metaDict, 1 ), "" );
	pcu_check_streq( Stg_Meta_GetParameterDocumentation( data->metaDict, 1 ), "A list of other VariableCondition definitions." );

	assocCount = Stg_Meta_GetAssociationCount( data->metaDict ); 
   pcu_check_true( assocCount == 1 );
	pcu_check_streq( Stg_Meta_GetAssociationName( data->metaDict, 0 ), "Data" );
	pcu_check_streq( Stg_Meta_GetAssociationType( data->metaDict, 0 ), "Stg_Component" );
	pcu_check_streq( Stg_Meta_GetAssociationNillable( data->metaDict, 0 ), "true" );
	pcu_check_streq( Stg_Meta_GetAssociationDocumentation( data->metaDict, 0 ), "User defined data in the form of a Stg_Component." );
} 


void Stg_MetaSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, Stg_MetaSuiteData );
   pcu_suite_setFixtures( suite, Stg_MetaSuite_Setup, Stg_MetaSuite_Teardown );
   pcu_suite_addTest( suite, Stg_MetaSuite_TestGetFunctions );
}
