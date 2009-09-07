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
** $Id:  $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>

const Type Underworld_ShapeFemIC_Type = "Underworld_ShapeFemIC";
typedef struct {
	__Codelet
} Underworld_ShapeFemIC;


void Underworld_LinearShapeIC( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {

  UnderworldContext*       context            = (UnderworldContext*)_context;
  Dictionary*              theDictionary      = context->dictionary;
  FeMesh*	           theMesh            = NULL;
  double*                  result             = (double*) _result;
  Stg_Shape*               shape;
  Name                     shapeName;
  double*                  coord;

  Dictionary_Entry_Value*  shapeSpecs;	
  Index                    numSpecs = 0;
  Index	                   shapeSpec_I;
  Dictionary_Entry_Value*  shapeSpec = NULL;
  Dictionary*	           shapeSpecDict;
  double                   ox,oy;
  double                   gx,gy,gz;
  double                   xp,yp,zp;    
  double                   dx,dy,dz;  
  double                   x,y,z;
  double                   T0;
  double                   az;
  double                   p1x,p1y,p1z,p2x,p2y,p2z;
  double                   a,b,c,d;
  double                   ux,uy,uz,nu;
  double                   gx2,gy2,gz2,T02;
  double                   alphax,alphay,alphaz;

  theMesh = context->temperatureField->feMesh;

  shapeSpecs = Dictionary_Get( theDictionary, "linearShapeIC" );
  numSpecs = Dictionary_Entry_Value_GetCount( shapeSpecs );

  for( shapeSpec_I = 0; shapeSpec_I < numSpecs; shapeSpec_I++ ) {
    shapeSpec = Dictionary_Entry_Value_GetElement( shapeSpecs, shapeSpec_I );
    shapeSpecDict = Dictionary_Entry_Value_AsDictionary( shapeSpec );
    shapeName = Dictionary_Entry_Value_AsString( Dictionary_Get( shapeSpecDict, "Shape" ) );

    /* Get the shape */
    shape = (Stg_Shape*) LiveComponentRegister_Get( context->CF->LCRegister, shapeName );

    Journal_Firewall( (Bool) shape, 
		      Journal_Register( Error_Type, Underworld_ShapeFemIC_Type ), 
		      "Shape %s not found.\n", shapeName );

    /* Find coordinate of node */
    coord = Mesh_GetVertex( theMesh, node_lI );

    if( Stg_Shape_IsCoordInside( shape, coord ) ) {

      /* -- Basic configuration -- */

      /* rotation angle */
      az = Dictionary_GetDouble_WithDefault( shapeSpecDict, "rotation", 0.0 );
      /* gradient in each direction */
      gx = Dictionary_GetDouble_WithDefault( shapeSpecDict, "gradientx", 0.0 );
      gy = Dictionary_GetDouble_WithDefault( shapeSpecDict, "gradienty", 0.0 );
      gz = Dictionary_GetDouble_WithDefault( shapeSpecDict, "gradientz", 0.0 );
      /* origin (2D case) */
      ox = Dictionary_GetDouble_WithDefault( shapeSpecDict, "originx", 0.0 );
      oy = Dictionary_GetDouble_WithDefault( shapeSpecDict, "originy", 0.0 );
      /* two points defining the rotational axis (3D case) First point is the origin */
      p1x = Dictionary_GetDouble_WithDefault( shapeSpecDict, "p1x", 0.0 );
      p1y = Dictionary_GetDouble_WithDefault( shapeSpecDict, "p1y", 0.0 );
      p1z = Dictionary_GetDouble_WithDefault( shapeSpecDict, "p1z", 0.0 );
      p2x = Dictionary_GetDouble_WithDefault( shapeSpecDict, "p2x", 0.0 );
      p2y = Dictionary_GetDouble_WithDefault( shapeSpecDict, "p2y", 0.0 );
      p2z = Dictionary_GetDouble_WithDefault( shapeSpecDict, "p2z", 1.0 );
      /* value of the field at origin */
      T0 = Dictionary_GetDouble_WithDefault( shapeSpecDict, "valueAtOrigin", 0.0 );

      /*  -- Advanced configuration -- */

      /* second gradient in each direction */
      gx2 = Dictionary_GetDouble_WithDefault( shapeSpecDict, "gradientx2", gx );
      gy2 = Dictionary_GetDouble_WithDefault( shapeSpecDict, "gradienty2", gy );
      gz2 = Dictionary_GetDouble_WithDefault( shapeSpecDict, "gradientz2", gz );
      /* point to impose 2nd gradients (relative to origin) */
      dx = Dictionary_GetDouble_WithDefault( shapeSpecDict, "distancex", 0.0 );
      dy = Dictionary_GetDouble_WithDefault( shapeSpecDict, "distancey", 0.0 );
      dz = Dictionary_GetDouble_WithDefault( shapeSpecDict, "distancez", 0.0 );      
      /* value of the field at origin + [dx,dy,dz] */
      T02 = Dictionary_GetDouble_WithDefault( shapeSpecDict, "valueAtOrigin2", T0 );

      x = coord[I_AXIS];
      y = coord[J_AXIS];
      az = az/180.0*M_PI;

      /* 2D case*/
      if( theMesh->topo->nDims == 2 ) {
	/* rotations and translations */
	xp = (x-ox)*cos(az) - (y-oy)*sin(az);
	yp = (x-ox)*sin(az) + (y-oy)*cos(az);

	/* compute gradients and values at origin*/
	if( dx != 0 ) 
	  alphax = xp/dx;
	else
	  alphax = 0;

	if( dy != 0 ) 
	  alphay = yp/dy;
	else
	  alphay = 0;

	gx = gx*(1.0-alphay) + gx2*alphay;
	gy = gy*(1.0-alphax) + gy2*alphax;
	T0 = T0*(1.0-alphax) + T02*alphax;

	/* compute value at the point */
	*result = T0 + xp*gx + yp*gy; 
      } 

	/* 3D case */
      if( theMesh->topo->nDims == 3 ) {
 	z   = coord[K_AXIS];

	/* unit axis of rotation */
	ux = p1x-p2x;
	uy = p1y-p2y;
	uz = p1z-p2z;
	nu = sqrt(ux*ux + uy*uy + uz*uz);
	ux /= nu;
	uy /= nu;
	uz /= nu;

	/* Rotational magic */
	d = sqrt( uy*uy + uz*uz );
	xp = cos(az)*( -uy*(y-p1y) - uz*(z-p1z) ) - sin(az)*( uz/d*(y-p1y) - uy/ux*(z-p1z) );
	yp = sin(az)*( -uy*(y-p1y) - uz*(z-p1z) ) + cos(az)*( uz/d*(y-p1y) - uy/ux*(z-p1z) );
	zp = d*uy/ux*(y-p1y) + d*uz/ux*(z-p1z);

	/* compute gradients and values at origin*/
	if( dx != 0 ) 
	  alphax = xp/dx;
	else
	  alphax = 0;

	if( dy != 0 ) 
	  alphay = yp/dy;
	else
	  alphay = 0;

	if( dx != 0 ) 
	  alphax = xp/dx;
	else
	  alphax = 0;

	gx = gx*(1.0-alphax) + gx2*alphax;
	gy = gy*(1.0-alphay) + gy2*alphay;
	gz = gz*(1.0-alphaz) + gz2*alphaz;
	T0 = T0*(1.0-alphax) + T02*alphax;

	/* compute value at the point */
	*result = T0 + xp*gx + yp*gy + zp*gz; 
      }
    }
  }
}


void Underworld_SimpleShapeIC( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	UnderworldContext* context    = (UnderworldContext*)_context;
	Dictionary* dictionary = context->dictionary;
	MeshVariable* meshVar    = NULL;
	FeMesh* mesh       = NULL;
	double* result     = (double*) _result;
	Stg_Shape* shape;
	Name shapeName;
	double* coord;

	meshVar = (MeshVariable*)Variable_Register_GetByIndex( context->variable_Register, var_I );
	mesh = (FeMesh*)meshVar->mesh; assert( mesh != NULL );

	shapeName = Dictionary_GetString( dictionary, "ShapeFemIC" );
	shape = (Stg_Shape*) LiveComponentRegister_Get( context->CF->LCRegister, shapeName );
	assert( shape );

	/* Find coordinate of node */
	coord = Mesh_GetVertex( mesh, node_lI );

	if ( Stg_Shape_IsCoordInside( shape, coord ) ) 
		*result = 1.0;
	else 
		*result = 0.0;
}

void Underworld_GaussianIC( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	UnderworldContext*      context            = (UnderworldContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	FeMesh*			mesh               = NULL;
	double*                 result             = (double*) _result;
	Stg_Shape*              shape;
	Name                    shapeName;
	double*                 coord;
	double                  disVec[3];
	double                  amplitude, width;
	double                  rSq;
	
	mesh       = context->temperatureField->feMesh;

	amplitude = Dictionary_GetDouble_WithDefault( dictionary, "GaussianIC-Amplitude", 1.0 );
	width = Dictionary_GetDouble_WithDefault( dictionary, "GaussianIC-Width", 1e-2 );

	shapeName = Dictionary_GetString( dictionary, "ShapeFemIC" );
	shape = (Stg_Shape*) LiveComponentRegister_Get( context->CF->LCRegister, shapeName );
	assert( shape );
	Journal_Firewall( !strcmp(shape->type, "Sphere") || !strcmp(shape->type, "Cylinder"),
			Journal_Register( Error_Type, Underworld_ShapeFemIC_Type ),
			"Error in %s: You're applying the GaussianIC to a shape of type %s, which can't be done."
			" It can only work on Sphere\' or \'Cylinder\' shapes\n", __func__,  shape->type );
	/* Find coordinate of node */
	coord = Mesh_GetVertex( mesh, node_lI );

	if( !strcmp(shape->type, "Sphere") ) {
		_Sphere_DistanceFromCenterAxis( shape, coord, disVec );

		rSq = disVec[0]*disVec[0]+disVec[1]*disVec[1];
		*result = amplitude * exp( -1 * rSq / (2 * width) );
	} 
	else if(  !strcmp(shape->type, "Cylinder") )  {
		_Cylinder_DistanceFromCenterAxis( shape, coord, disVec );

		rSq = disVec[0]*disVec[0]+disVec[1]*disVec[1];
		if( shape->dim == 3 ) rSq += disVec[2]*disVec[2];

		*result = amplitude * exp( -1 * rSq / (2 * width) );
	}

}

void _Underworld_ShapeFemIC_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	ConditionFunction*      condFunc;
	UnderworldContext*      context;

	context = (UnderworldContext*)Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data ); 
	
	condFunc = ConditionFunction_New( Underworld_SimpleShapeIC, "Inside1_Outside0_ShapeIC" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( Underworld_GaussianIC, "GaussianIC" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( Underworld_LinearShapeIC, "linearShapeIC" );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
}

void* _Underworld_ShapeFemIC_DefaultNew( Name name ) {
	return Codelet_New(
		Underworld_ShapeFemIC_Type,
		_Underworld_ShapeFemIC_DefaultNew,
		_Underworld_ShapeFemIC_Construct,
		_Codelet_Build,
		_Codelet_Initialise,
		_Codelet_Execute,
		_Codelet_Destroy,
		name );
}

Index Underworld_ShapeFemIC_Register( PluginsManager* pluginsManager ) {
	Journal_DPrintf( Underworld_Debug, "In: %s( void* )\n", __func__ );

	return PluginsManager_Submit( pluginsManager, Underworld_ShapeFemIC_Type, "0", _Underworld_ShapeFemIC_DefaultNew );
}

