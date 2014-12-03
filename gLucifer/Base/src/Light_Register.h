/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2005, Monash Cluster Computing 
** All rights reserved.
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
** 		* Redistributions of source code must retain the above copyright notice, 
** 			this list of conditions and the following disclaimer.
** 		* Redistributions in binary form must reproduce the above copyright 
**			notice, this list of conditions and the following disclaimer in the 
**			documentation and/or other materials provided with the distribution.
** 		* Neither the name of the Monash University nor the names of its contributors 
**			may be used to endorse or promote products derived from this software 
**			without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
** THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
** PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS 
** BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
** CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
** SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
** HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
** LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
** OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
**
** Contact:
*%		Cecile Duboz - Cecile.Duboz@sci.monash.edu.au
*%
** Contributors:
*+		Cecile Duboz
*+		Robert Turnbull
*+		Alan Lo
*+		Louis Moresi
*+		David Stegman
*+		David May
*+		Stevan Quenette
*+		Patrick Sunter
*+		Greg Watson
*+
** $Id: Arrhenius.c 78 2005-11-29 11:58:21Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
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
**
** This file may be distributed under the terms of the VPAC Public License
** as defined by VPAC of Australia and appearing in the file
** LICENSE.VPL included in the packaging of this file.
**
** This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
** WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
**
*/
/** \file
**  Role:
**	Allows users to access lucLights based on their textual name,
**	or index.
**
** Assumptions:
**
** Comments:
**
**
** $Id: Light_Register.h 419 2005-11-29 12:04:20Z RobertTurnbull $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __lucLight_Register_h__
#define __lucLight_Register_h__
	
	extern const Type lucLight_Register_Type;
	
	#define __lucLight_Register \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__NamedObject_Register \
		\
		/* Virtual functions go here */ \
		\
		/* Class info */ \
		Light_Index currentLightIndex;

	struct lucLight_Register { __lucLight_Register };
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructor
	*/
	
	lucLight_Register*	lucLight_Register_New( void );
	
	/*-----------------------------------------------------------------------------------------------------------------
	** General virtual functions
	*/
	
	/*-----------------------------------------------------------------------------------------------------------------
	** Public functions
	*/
	#define lucLight_Register_Add NamedObject_Register_Add

	#define lucLight_Register_GetIndex NamedObject_Register_GetIndex 

	#define lucLight_Register_GetByName( self, materialName ) \
		( (lucLight*) NamedObject_Register_GetByName( self, materialName ) ) 

	#define lucLight_Register_GetByIndex( self, materialIndex ) \
		( (lucLight*) NamedObject_Register_GetByIndex( self, materialIndex ) )

	#define lucLight_Register_GetCount( self ) \
		(self)->objects->count
	
	#define lucLight_Register_DeleteAllObjects( self ) \
		Stg_ObjectList_DeleteAllObjects( (self)->objects )
	#define lucLight_Register_PrintAllObjects( self, stream ) \
		Stg_ObjectList_PrintAllObjects( (self)->objects, stream )

	/* +++ Public Functions +++ */
        void  lucLight_Register_EnableAll( void * lightRegister );
	Light_Index   lucLight_Register_GetCurrentLightIndex( void * lightRegister );
	void lucLight_Register_SetCurrentLightIndex( void * lightRegister, Light_Index index );
	void lucLight_Register_ChangeCurrentLightIndex( void * lightRegister );


        
#endif 
