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
*%		Louis Moresi - Louis.Moresi@sci.monash.edu.au
*%
** Contributors:
*+		Robert Turnbull
*+		Vincent Lemiale
*+		Louis Moresi
*+		David May
*+		David Stegman
*+		Mirko Velic
*+		Patrick Sunter
*+		Julian Giordani
*+
** $Id: IncompressibleExtensionBC.h 694 2008-04-03 01:06:30Z LouisMoresi $
** 
** Comments:
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __Underworld_IncompressibleExtensionBC_h__
#define __Underworld_IncompressibleExtensionBC_h__

//	double GetLeftWallVelocity( FeVariable* velocityField ) ;
//	double GetRightWallVelocity( FeVariable* velocityField ) ;
	double GetTopWallVelocity( FeVariable* velocityField, void* _context, double y ) ;
	double GetBottomWallVelocity( FeVariable* velocityField, void* _context, double y ) ;
//	void GetVelocity( FeVariable* velocityField, double y, Coord coord, void* _context, double* velocity ) ;
	
	void IncompressibleExtensionBC_TopCondition( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
	void IncompressibleExtensionBC_LeftCondition( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
	void IncompressibleExtensionBC_RightCondition( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;
	void IncompressibleExtensionBC_BottomCondition( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) ;

	Index Underworld_IncompressibleExtensionBC_Register( PluginsManager* pluginsManager );
	void* _Underworld_IncompressibleExtensionBC_DefaultNew( Name name );
	void Underworld_IncompressibleExtensionBC_Construct( void* component, Stg_ComponentFactory* _cf, void* data ) ;

#endif
