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
** $Id: ContactVC.c 4153 2007-07-26 02:25:22Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>

#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>

#include "types.h"
#include "WallVC.h"
#include "ContactVC.h"
#include "RegularMeshUtils.h"

#include <string.h>
#include <assert.h>


const Type ContactVC_Type = "ContactVC";
const Name defaultContactVCName = "defaultContactVCName";


/*--------------------------------------------------------------------------------------------------------------------------
** Constructor
*/

VariableCondition* ContactVC_Factory(
	AbstractContext*					context,
	Variable_Register*				variable_Register, 
   ConditionFunction_Register*	conFunc_Register, 
   Dictionary*							dictionary,
   void*									data )
{
   return (VariableCondition*)ContactVC_New( defaultContactVCName, context, NULL, variable_Register, conFunc_Register, dictionary, (Mesh*)data );
}

ContactVC*	ContactVC_New(
   Name									name,
	AbstractContext*					context,
   Name									_dictionaryEntryName, 
   Variable_Register*				variable_Register, 
   ConditionFunction_Register*	conFunc_Register, 
   Dictionary*							dictionary,
   void*									_mesh )
{
   ContactVC* self = _ContactVC_DefaultNew( name );
	
	_VariableCondition_Init( self, context, variable_Register, conFunc_Register, dictionary );
   _WallVC_Init( self, _dictionaryEntryName, _mesh );
	_ContactVC_Init( self, _dictionaryEntryName, _mesh );

	return self;
}

ContactVC* _ContactVC_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                               _sizeOfSelf = sizeof(ContactVC);
	Type                                                       type = ContactVC_Type;
	Stg_Class_DeleteFunction*                               _delete = _ContactVC_Delete;
	Stg_Class_PrintFunction*                                 _print = _WallVC_Print;
	Stg_Class_CopyFunction*                                   _copy = _WallVC_Copy;
	Stg_Component_DefaultConstructorFunction*   _defaultConstructor = (Stg_Component_DefaultConstructorFunction*)_ContactVC_DefaultNew;
	Stg_Component_ConstructFunction*                     _construct = _ContactVC_AssignFromXML;
	Stg_Component_BuildFunction*                             _build = _ContactVC_Build;
	Stg_Component_InitialiseFunction*                   _initialise = _VariableCondition_Initialise;
	Stg_Component_ExecuteFunction*                         _execute = _VariableCondition_Execute;
	Stg_Component_DestroyFunction*                         _destroy = _WallVC_Destroy;
	AllocationType                               nameAllocationType = NON_GLOBAL;
	VariableCondition_BuildSelfFunc*                     _buildSelf = _WallVC_BuildSelf;
	VariableCondition_PrintConciseFunc*               _printConcise = _WallVC_PrintConcise;
	VariableCondition_ReadDictionaryFunc*           _readDictionary = _ContactVC_ReadDictionary;
	VariableCondition_GetSetFunc*                           _getSet = _ContactVC_GetSet;
	VariableCondition_GetVariableCountFunc*       _getVariableCount = _WallVC_GetVariableCount;
	VariableCondition_GetVariableIndexFunc*       _getVariableIndex = _WallVC_GetVariableIndex;
	VariableCondition_GetValueIndexFunc*             _getValueIndex = _WallVC_GetValueIndex;
	VariableCondition_GetValueCountFunc*             _getValueCount = _WallVC_GetValueCount;
	VariableCondition_GetValueFunc*                       _getValue = _WallVC_GetValue;
	VariableCondition_ApplyFunc*                             _apply = _VariableCondition_Apply;

   return _ContactVC_New(  CONTACTVC_PASSARGS  );
}

ContactVC* _ContactVC_New(  CONTACTVC_DEFARGS  ) {
   ContactVC* self;
	
   /* Allocate memory/General info */
   assert( _sizeOfSelf >= sizeof(ContactVC) );
   self = (ContactVC*)_WallVC_New(  WALLVC_PASSARGS  );
	
   /* Virtual info */
	
   /* Stg_Class info */
	
   return self;
}


void _ContactVC_Init( void* wallVC, Name _dictionaryEntryName, void* _mesh ) {
   ContactVC* self = (ContactVC*)wallVC;

   self->deep = False;
}


/*--------------------------------------------------------------------------------------------------------------------------
** General virtual functions
*/

void _ContactVC_ReadDictionary( void* variableCondition, void* dictionary ) {
   ContactVC* self = (ContactVC*)variableCondition;
   Dictionary_Entry_Value *vcDictVal, _vcDictVal, *entryVal;

   _WallVC_ReadDictionary( variableCondition, dictionary );

   /* Find dictionary entry */
   if (self->_dictionaryEntryName)
      vcDictVal = Dictionary_Get( dictionary, (Dictionary_Entry_Key)self->_dictionaryEntryName );
   else
   {
      vcDictVal = &_vcDictVal;
      Dictionary_Entry_Value_InitFromStruct(vcDictVal, dictionary);
   }

   if (vcDictVal) {
      entryVal = Dictionary_Entry_Value_GetMember( vcDictVal, (Dictionary_Entry_Key)"deep" );
      if( entryVal )
         self->deep = Dictionary_Entry_Value_AsBool( entryVal );
   }
   else {
      self->deep = False;
   }
}


void _ContactVC_Delete(void* wallVC) {
   ContactVC*	self = (ContactVC*)wallVC;
	
   /* Stg_Class_Delete parent */
   _WallVC_Delete( self  );
}

void _ContactVC_Destroy(void* wallVC, void* data) {
   ContactVC*	self = (ContactVC*)wallVC;
	
   /* Stg_Class_Delete parent */
   _WallVC_Destroy( self, data);
}

void _ContactVC_Build(  void* wallVC, void* data ) {
   ContactVC*			self = (ContactVC*)wallVC;
	
   _WallVC_Build( self, data );
}
	

/*--------------------------------------------------------------------------------------------------------------------------
** Macros
*/


/*--------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _ContactVC_AssignFromXML( void* wallVC, Stg_ComponentFactory* cf, void* data ) { 
}

IndexSet* _ContactVC_GetSet(void* variableCondition) {
   ContactVC*				self = (ContactVC*)variableCondition;
   IndexSet*				set = NULL;
   Stream*					warningStr = Journal_Register( Error_Type, (Name)self->type );
   unsigned					nDims;
   Grid*						vertGrid;
   CartesianGenerator*	gen;

   nDims = Mesh_GetDimSize( self->_mesh );
   gen = (CartesianGenerator* )self->_mesh->generator;

   if( strcmp( gen->type, CartesianGenerator_Type ) )
      abort();
   vertGrid = *(Grid**)ExtensionManager_Get( self->_mesh->info, self->_mesh, ExtensionManager_GetHandle( self->_mesh->info, (Name)"vertexGrid" ) );

   switch (self->_wall) {
      case WallVC_Wall_Front:
         if ( nDims < 3 || !vertGrid->sizes[2]  ) {
            Journal_Printf( warningStr, "Warning - in %s: Can't build a %s wall VC "
                            "when mesh has no elements in the %s axis. Returning an empty set.\n", __func__,
                            WallVC_WallEnumToStr[self->_wall], "K" );
            set = IndexSet_New( Mesh_GetDomainSize( self->_mesh, MT_VERTEX ) );
         }
         else {
            abort();
            /*set = RegularMeshUtils_CreateContactFrontSet( self->_mesh );*/
         }
         break;
			
      case WallVC_Wall_Back:
         if ( nDims < 3 || !vertGrid->sizes[2] ) {
            Journal_Printf( warningStr, "Warning - in %s: Can't build a %s wall VC "
                            "when mesh has no elements in the %s axis. Returning an empty set.\n", __func__,
                            WallVC_WallEnumToStr[self->_wall], "K" );
            set = IndexSet_New( Mesh_GetDomainSize( self->_mesh, MT_VERTEX ) );
         }
         else {
            abort();
            /*set = RegularMeshUtils_CreateContactBackSet( self->_mesh );*/
         }	
         break;
			
      case WallVC_Wall_Top:
         if ( nDims < 2 || !vertGrid->sizes[1] ) {
            Journal_Printf( warningStr, "Warning - in %s: Can't build a %s wall VC "
                            "when mesh has no elements in the %s axis. Returning an empty set.\n", __func__,
                            WallVC_WallEnumToStr[self->_wall], "J" );
            set = IndexSet_New( Mesh_GetDomainSize( self->_mesh, MT_VERTEX ) );
         }
         else {
            set = RegularMeshUtils_CreateContactTopSet(self->_mesh, gen->contactDepth[0][0], gen->contactDepth[0][1]);
         }	
         break;
			
      case WallVC_Wall_Bottom:
         if ( nDims < 2 || !vertGrid->sizes[1] ) {
            Journal_Printf( warningStr, "Warning - in %s: Can't build a %s wall VC "
                            "when mesh has no elements in the %s axis. Returning an empty set.\n", __func__,
                            WallVC_WallEnumToStr[self->_wall], "J" );
            set = IndexSet_New( Mesh_GetDomainSize( self->_mesh, MT_VERTEX ) );
         }
         else {
            if( self->deep ) {
               set = RegularMeshUtils_CreateContactBottomSet(
                  self->_mesh, 0, 0, gen->contactDepth[1][0] );
            }
            else {
               set = RegularMeshUtils_CreateContactBottomSet(
                  self->_mesh, gen->contactDepth[0][0], gen->contactDepth[0][1], 0 );
            }
         }	
         break;
			
      case WallVC_Wall_Left:
         if ( nDims < 1 ) {
            Journal_Printf( warningStr, "Warning - in %s: Can't build a %s wall VC "
                            "when mesh has no elements in the %s axis. Returning an empty set.\n", __func__,
                            WallVC_WallEnumToStr[self->_wall], "I" );
            set = IndexSet_New( Mesh_GetDomainSize( self->_mesh, MT_VERTEX ) );
         }
         else {
            set = RegularMeshUtils_CreateContactLeftSet(self->_mesh, gen->contactDepth[1][0], gen->contactDepth[1][1]);
         }	
         break;
			
      case WallVC_Wall_Right:
         if( nDims < 1 ) {
            Journal_Printf( warningStr, "Warning - in %s: Can't build a %s wall VC "
                            "when mesh has no elements in the %s axis. Returning an empty set.\n", __func__,
                            WallVC_WallEnumToStr[self->_wall], "I" );
            set = IndexSet_New( Mesh_GetDomainSize( self->_mesh, MT_VERTEX ) );
         }
         else {
            if( self->deep ) {
               set = RegularMeshUtils_CreateContactRightSet(
                  self->_mesh, 0, 0, gen->contactDepth[0][1] );
            }
            else {
               set = RegularMeshUtils_CreateContactRightSet(
                  self->_mesh, gen->contactDepth[1][0], gen->contactDepth[1][1], 0 );
            }
         }
         break;
			
      case WallVC_Wall_Size:
      default:
         assert(0);
         break;
   }
	
   return set;
}

/*--------------------------------------------------------------------------------------------------------------------------
** Build functions
*/


/*--------------------------------------------------------------------------------------------------------------------------
** Functions
*/


