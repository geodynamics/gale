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

#define STD   1
#define PM    2

const Type Underworld_ShapeFemIC_Type = "Underworld_ShapeFemIC";
typedef struct {
	__Codelet
} Underworld_ShapeFemIC;



void Underworld_LinearShapeIC( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {

  UnderworldContext*       context            = (UnderworldContext*)_context;

  Dictionary*              theDictionary      = context->dictionary;
  FeVariable*    tempField   = (FeVariable*)LiveComponentRegister_Get( context->CF->LCRegister, (Name)"TemperatureField" );
  FeMesh*	           theMesh            = NULL;
  double*                  result             = (double* ) _result;
  Stg_Shape*               shape;
  Name                     shapeName;
  double*                  coord;
  Dictionary_Entry_Value*  shapeSpecs;	
  Index                    numSpecs = 0;
  Index	                   shapeSpec_I;
  Dictionary_Entry_Value*  shapeSpec = NULL;
  Dictionary*	           shapeSpecDict;

  int setup;
  /* STD config */
  double gz, p1x, p1y, p1z, p2x, p2y, p2z, z, ux, uy, uz, nu, d, zp;
  double gx, gy, T0, x, y, az, xp, yp;
  double ox, oy;
  /* PM config */
  double ox1,oy1,T1,H1,gy1;
  double ox2,oy2,T2,H2,gy2;
  double T1u, T1l;
  double T2u, T2l;
  double Tmu, Tml;
  double o1ux, o1uy;
  double o1ly;
  double o2ux, o2uy;
  double o2ly;
  double omuy, omly;
  double alpha, beta, W;

  theMesh = tempField->feMesh;

  shapeSpecs = Dictionary_Get( theDictionary, (Dictionary_Entry_Key)"linearShapeIC" );
  numSpecs = Dictionary_Entry_Value_GetCount( shapeSpecs );

  for( shapeSpec_I = 0; shapeSpec_I < numSpecs; shapeSpec_I++  ) {
    shapeSpec = Dictionary_Entry_Value_GetElement( shapeSpecs, shapeSpec_I );
    shapeSpecDict = Dictionary_Entry_Value_AsDictionary( shapeSpec );
    shapeName = Dictionary_Entry_Value_AsString( Dictionary_Get( shapeSpecDict, (Dictionary_Entry_Key)"Shape" ) );

    /* Get the shape */
    shape = (Stg_Shape* ) LiveComponentRegister_Get( context->CF->LCRegister, (Name)shapeName );

    Journal_Firewall( shape!=NULL, 
		      Journal_Register( Error_Type, (Name)Underworld_ShapeFemIC_Type  ), 
		      "Shape %s not found.\n", shapeName );

    /* Find coordinate of node */
    coord = Mesh_GetVertex( theMesh, node_lI );

    if( Stg_Shape_IsCoordInside( shape, coord ) ) {

      setup = Dictionary_GetInt_WithDefault( shapeSpecDict, (Dictionary_Entry_Key)"setup", STD );

      switch( setup  ) {
      case STD:
	/* rotation angle */
	az = Dictionary_GetDouble_WithDefault( shapeSpecDict, (Dictionary_Entry_Key)"rotation", 0.0  );
	/* gradient in each direction */
	gx = Dictionary_GetDouble_WithDefault( shapeSpecDict, (Dictionary_Entry_Key)"gradientx", 0.0  );
	gy = Dictionary_GetDouble_WithDefault( shapeSpecDict, (Dictionary_Entry_Key)"gradienty", 0.0  );
	/* value of the field at origin */
	T0 = Dictionary_GetDouble_WithDefault( shapeSpecDict, (Dictionary_Entry_Key)"valueAtOrigin", 0.0 );

	x = coord[I_AXIS];
	y = coord[J_AXIS];
	az = az/180.0*M_PI;

	switch( theMesh->topo->nDims  ) {
	case 2: 
	  /* origin */
	  ox = Dictionary_GetDouble_WithDefault( shapeSpecDict, (Dictionary_Entry_Key)"originx", 0.0  );
	  oy = Dictionary_GetDouble_WithDefault( shapeSpecDict, (Dictionary_Entry_Key)"originy", 0.0 );
	  /* rotations and translations */
	  xp = (x-ox)*cos(az) - (y-oy)*sin(az);
	  yp = (x-ox)*sin(az) + (y-oy)*cos(az );
	  /* compute value at the point */
	  *result = T0 + xp*gx + yp*gy; 
	  break;
	
	case 3: 
	  gz = Dictionary_GetDouble_WithDefault( shapeSpecDict, (Dictionary_Entry_Key)"gradientz", 0.0 );
	  /* two points defining the rotational axis (3D case ) First point is the origin */
	  p1x = Dictionary_GetDouble_WithDefault( shapeSpecDict, (Dictionary_Entry_Key)"p1x", 0.0  );
	  p1y = Dictionary_GetDouble_WithDefault( shapeSpecDict, (Dictionary_Entry_Key)"p1y", 0.0  );
	  p1z = Dictionary_GetDouble_WithDefault( shapeSpecDict, (Dictionary_Entry_Key)"p1z", 0.0  );
	  p2x = Dictionary_GetDouble_WithDefault( shapeSpecDict, (Dictionary_Entry_Key)"p2x", 0.0  );
	  p2y = Dictionary_GetDouble_WithDefault( shapeSpecDict, (Dictionary_Entry_Key)"p2y", 0.0  );
	  p2z = Dictionary_GetDouble_WithDefault( shapeSpecDict, (Dictionary_Entry_Key)"p2z", 1.0 );

	  z   = coord[K_AXIS];

	  /* unit axis of rotation */
	  ux = p1x-p2x;
	  uy = p1y-p2y;
	  uz = p1z-p2z;
	  nu = sqrt(ux*ux + uy*uy + uz*uz);

	  ux /= nu;
	  uy /= nu;
	  uz /= nu;
	  
	  /* Rotate and translate */
	  d = sqrt( uy*uy + uz*uz );
	  xp = cos(az)*( -uy*(y-p1y) - uz*(z-p1z) ) - sin(az)*( uz/d*(y-p1y) - uy/ux*(z-p1z) );
	  yp = sin(az)*( -uy*(y-p1y) - uz*(z-p1z) ) + cos(az)*( uz/d*(y-p1y) - uy/ux*(z-p1z) );
	  zp = d*uy/ux*(y-p1y) + d*uz/ux*(z-p1z );

	  /* compute value at the point */
	  *result = T0 + xp*gx + yp*gy + zp*gz; 
	  break;	
	}
	break;

      case PM:
	/* Advanced config */
	ox1 = Dictionary_GetDouble_WithDefault( shapeSpecDict, (Dictionary_Entry_Key)"originx1", 0.0  );
	oy1 = Dictionary_GetDouble_WithDefault( shapeSpecDict, (Dictionary_Entry_Key)"originy1", 0.0  );
	gy1 = Dictionary_GetDouble_WithDefault( shapeSpecDict, (Dictionary_Entry_Key)"gradienty1", 0.0  );
	T1 = Dictionary_GetDouble_WithDefault( shapeSpecDict, (Dictionary_Entry_Key)"valueAtOrigin1", 0.0  );
	H1 = Dictionary_GetDouble_WithDefault( shapeSpecDict, (Dictionary_Entry_Key)"height1", 0.0  );

	ox2 = Dictionary_GetDouble_WithDefault( shapeSpecDict, (Dictionary_Entry_Key)"originx2", 0.0  );
	oy2 = Dictionary_GetDouble_WithDefault( shapeSpecDict, (Dictionary_Entry_Key)"originy2", 0.0  );
	gy2 = Dictionary_GetDouble_WithDefault( shapeSpecDict, (Dictionary_Entry_Key)"gradienty2", 0.0  );
	T2 = Dictionary_GetDouble_WithDefault( shapeSpecDict, (Dictionary_Entry_Key)"valueAtOrigin2", 0.0  );
	H2 = Dictionary_GetDouble_WithDefault( shapeSpecDict, (Dictionary_Entry_Key)"height2", 0.0 );

	x = coord[I_AXIS];
	y = coord[J_AXIS];

	T1u = T1;
	T1l = T1 + H1 * gy1;
	T2u = T2;
	T2l = T2 + H2 * gy2;

	o1ux = ox1;
	o1uy = oy1;
	o1ly = oy1 + H1;
            
	o2ux = ox2;
	o2uy = oy2;
	o2ly = oy2 + H2;

	W = o2ux - o1ux;
	alpha = (x - o1ux) / W;
            
	Tmu = (1-alpha)*T1u + alpha*T2u;
	Tml = (1-alpha)*T1l + alpha*T2l;            
	omuy = (1-alpha)*o1uy + alpha*o2uy;
	omly = (1-alpha)*o1ly + alpha*o2ly;
            
	beta = (y - omuy) / (omly - omuy);

	*result = (1-beta )*Tmu + beta*Tml;
	break;
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

	shapeName = Dictionary_GetString( dictionary, (Dictionary_Entry_Key)"ShapeFemIC" );
	shape = (Stg_Shape* ) LiveComponentRegister_Get( context->CF->LCRegister, (Name)shapeName );
	assert( shape  );

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
  FeVariable*    tempField   = (FeVariable*)LiveComponentRegister_Get( context->CF->LCRegister, (Name)"TemperatureField" );
	FeMesh*			mesh               = NULL;
	double*                 result             = (double* ) _result;
	Stg_Shape*              shape;
	Name                    shapeName;
	double*                 coord;
	double                  disVec[3];
	double                  amplitude, width;
	double                  rSq;
	
	mesh       = tempField->feMesh;

	amplitude = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"GaussianIC-Amplitude", 1.0  );
	width = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"GaussianIC-Width", 1e-2  );

	shapeName = Dictionary_GetString( dictionary, (Dictionary_Entry_Key)"ShapeFemIC" );
	shape = (Stg_Shape* ) LiveComponentRegister_Get( context->CF->LCRegister, (Name)shapeName );
	assert( shape  );
	Journal_Firewall( !strcmp(shape->type, "Sphere") || !strcmp(shape->type, "Cylinder"),
			Journal_Register( Error_Type, (Name)Underworld_ShapeFemIC_Type  ),
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

void _Underworld_ShapeFemIC_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
	ConditionFunction*      condFunc;
	UnderworldContext*      context;

	context = (UnderworldContext*)Stg_ComponentFactory_ConstructByName( cf, (Name)"context", UnderworldContext, True, data  ); 
	
	condFunc = ConditionFunction_New( Underworld_SimpleShapeIC, (Name)"Inside1_Outside0_ShapeIC"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( Underworld_GaussianIC, (Name)"GaussianIC"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
	condFunc = ConditionFunction_New( Underworld_LinearShapeIC, (Name)"linearShapeIC"  );
	ConditionFunction_Register_Add( condFunc_Register, condFunc );
}

void* _Underworld_ShapeFemIC_DefaultNew( Name name ) {
	return Codelet_New(
		Underworld_ShapeFemIC_Type,
		_Underworld_ShapeFemIC_DefaultNew,
		_Underworld_ShapeFemIC_AssignFromXML,
		_Codelet_Build,
		_Codelet_Initialise,
		_Codelet_Execute,
		_Codelet_Destroy,
		name );
}

Index Underworld_ShapeFemIC_Register( PluginsManager* pluginsManager ) {
	Journal_DPrintf( Underworld_Debug, "In: %s( void* )\n", __func__ );

	return PluginsManager_Submit( pluginsManager, Underworld_ShapeFemIC_Type, (Name)"0", _Underworld_ShapeFemIC_DefaultNew  );
}



