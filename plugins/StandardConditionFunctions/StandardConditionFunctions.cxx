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
** $Id: StandardConditionFunctions.c 1196 2008-08-04 16:29:30Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <string.h>
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <assert.h>
#include "StandardConditionFunctions.h"
#include <list>
#include "boost/math/special_functions/erf.hpp"

const Type StgFEM_StandardConditionFunctions_Type = "StgFEM_StandardConditionFunctions";

void _StgFEM_StandardConditionFunctions_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
  Codelet*		self		= (Codelet*)component;
  AbstractContext*        context;
  ConditionFunction*      condFunc;
  Dictionary*		pluginDict	= Codelet_GetPluginDictionary( component, cf->rootDict );

  context = (AbstractContext*)Stg_ComponentFactory_ConstructByName( cf, Dictionary_GetString( pluginDict, (Dictionary_Entry_Key)"Context"  ), AbstractContext, True, data );
  self->context = context;
	
  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SolidBodyRotation, "Velocity_SolidBodyRotation"  );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
	
  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_PartialRotationX, "Velocity_PartialRotationX"  );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
		
  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_PartialRotationY, "Velocity_PartialRotationY"  );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
	
  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_TaperedRotationX, "TaperedRotationX" );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
		
  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_TaperedRotationY, "TaperedRotationY" );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
	
  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SimpleShear, "Velocity_SimpleShear" );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SimpleShearInverted, "Velocity_SimpleShearInverted" );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_ShearZ, "ShearZ" );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_Extension, "Velocity_Extension" );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_PartialLid_TopLayer, "Velocity_PartialLid_TopLayer"  );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_Trigonometry, "Temperature_Trigonometry"  );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_LinearInterpolationLid, "Velocity_LinearInterpolationLid"  );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_Lid_RampWithCentralMax, "Velocity_Lid_RampWithCentralMax"  );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
	
  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_LinearVelocityLeftWall, "LinearVelocityLeftWall"  );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
	
  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_LinearVelocityRightWall, "LinearVelocityRightWall"  );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SinusoidalLid, "Velocity_SinusoidalLid"  );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_TemperatureCosineHill, "Temperature_CosineHill"  );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
	
  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_ConvectionBenchmark, "Temperature_ConvectionBenchmark"  );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
	
  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_LinearWithSinusoidalPerturbation, "LinearWithSinusoidalPerturbation"  );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_EdgeDriveConvectionIC, "EdgeDriveConvectionIC"  );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_ThermalEdgeDriveConvectionIC, "ThermalEdgeDriveConvectionIC"  );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
	
  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_AnalyticalTemperatureIC, "AnalyticalTemperatureIC"  );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New( Stg_FEM_VelicTemperatureIC, "VelicTemperatureIC" );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New( Stg_FEM_VelicTemperatureIC_SolB, "VelicTemperatureIC_SolB" );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
	
  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SinusoidalExtension, "SinusoidalExtension" );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_StepFunction, "StepFunction" );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StG_FEM_StandardConditionFunctions_StepFunctionProduct1, "StepFunctionProduct1");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StG_FEM_StandardConditionFunctions_StepFunctionProduct2, "StepFunctionProduct2");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StG_FEM_StandardConditionFunctions_StepFunctionProduct3, "StepFunctionProduct3");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StG_FEM_StandardConditionFunctions_StepFunctionProduct4, "StepFunctionProduct4");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_TemperatureProfile, "TemperatureProfile");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StG_FEM_StandardConditionFunctions_Gaussian, "Gaussian");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_ERF,
                                   (char*)"ERF");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_ERFC,
                                   (char*)"ERFC");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_RubberSheet,
                                   (char*)"RubberSheet");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_MovingStepFunction, "MovingStepFunction");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_SpecRidge3D, "SpecRidge3D" );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_ErrorFunc, "ErrorFunc" );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_ConstantVector, "ConstantVector" );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_GaussianDistribution, "GaussianDistribution" );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_1DGaussianDistribution, "1DGaussianDistribution" );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_HalfContainer, "HalfContainer" );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_ConstantValue, "ConstantValue" );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_DiagonalLine, "DiagonalLine" );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_DeltaFunction, "DeltaFunction" );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_InflowBottom, "InflowBottom" );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_GaussianTube, "GaussianTube" );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New( StgFEM_StandardConditionFunctions_GravitationalPotential, "GravitationalPotential" );
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_WarsTemperature,
                                   "WarsTemperature");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_Quadratic,
                                   "Quadratic");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );

  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File0,
                                   "File0");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File1,
                                   "File1");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File2,
                                   "File2");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File3,
                                   "File3");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File4,
                                   "File4");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File5,
                                   "File5");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File6,
                                   "File6");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File7,
                                   "File7");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File8,
                                   "File8");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File9,
                                   "File9");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File10,
                                   "File10");
  ConditionFunction_Register_Add(condFunc_Register,condFunc);
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File11,
                                   "File11");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File12,
                                   "File12");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File13,
                                   "File13");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File14,
                                   "File14");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File15,
                                   "File15");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File16,
                                   "File16");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File17,
                                   "File17");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File18,
                                   "File18");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File19,
                                   "File19");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File20,
                                   "File20");
  ConditionFunction_Register_Add(condFunc_Register,condFunc);
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File21,
                                   "File21");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File22,
                                   "File22");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File23,
                                   "File23");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File24,
                                   "File24");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File25,
                                   "File25");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File26,
                                   "File26");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File27,
                                   "File27");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File28,
                                   "File28");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
  condFunc = ConditionFunction_New(StgFEM_StandardConditionFunctions_File29,
                                   "File29");
  ConditionFunction_Register_Add( condFunc_Register, condFunc );
}

void _StgFEM_StandardConditionFunctions_Destroy( void* _self, void* data ) {
   /* This function will totally clean the condFunc_Register
    *
    * This could be trouble some if other code adds or deletes condition functions on this register
    */

   unsigned *refCount = &(condFunc_Register->count);

   /* first check if there are things still on the condFunc_Register, if so .... */
   if( *refCount != 0 ) {
      while( *refCount != 0 ) {

         _ConditionFunction_Delete( condFunc_Register->_cf[ *refCount-1 ] );
         condFunc_Register->_cf[ *refCount-1 ] = NULL;

         *refCount = *refCount - 1;
      }
   }
   _Codelet_Destroy( _self, data );
}
void* _StgFEM_StandardConditionFunctions_DefaultNew( Name name ) {
	return Codelet_New(
		StgFEM_StandardConditionFunctions_Type,
		_StgFEM_StandardConditionFunctions_DefaultNew,
		_StgFEM_StandardConditionFunctions_AssignFromXML,
		_Codelet_Build,
		_Codelet_Initialise,
		_Codelet_Execute,
		_StgFEM_StandardConditionFunctions_Destroy,
		name );
}

Index StgFEM_StandardConditionFunctions_Register( PluginsManager* pluginsManager ) {
	Journal_DPrintf( StgFEM_Debug, "In: %s( void* )\n", __func__ );

	return PluginsManager_Submit( pluginsManager, StgFEM_StandardConditionFunctions_Type, (Name)"0", _StgFEM_StandardConditionFunctions_DefaultNew  );
}

Bool StgFEM_StandardConditionFunctions_Init( int* argc, char** argv[] ) {
  Stg_ComponentRegister* componentsRegister = Stg_ComponentRegister_Get_ComponentRegister();
  Stg_ComponentRegister_Add(componentsRegister,
                            StgFEM_StandardConditionFunctions_Type, (Name)"0",
                            _StgFEM_StandardConditionFunctions_DefaultNew );
  RegisterParent( StgFEM_StandardConditionFunctions_Type, Stg_Component_Type );
  return True;
}

void StgFEM_StandardConditionFunctions_SolidBodyRotation(const double *coord, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	Coord                   centre;
	Coord                   vector;
	double                  omega;
	
	/* Find Centre of Solid Body Rotation */
	centre[ I_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationCentreX", 0.0  );
	centre[ J_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationCentreY", 0.0  );
	centre[ K_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationCentreZ", 0.0  );
	omega            = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationOmega", 1.0  );

	/* Find vector from centre to node */
	StGermain_VectorSubtraction( vector, coord, centre, 2 );

	result[ I_AXIS ] = -omega * vector[ J_AXIS ];
	result[ J_AXIS ] =  omega * vector[ I_AXIS ];
}


void StgFEM_StandardConditionFunctions_PartialRotationX(const double *coord, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	Coord                   centre;
	Coord                   vector;
	double                  omega;
	double			size;

	/* Find Centre of Solid Body Rotation */
	centre[ I_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationCentreX", 0.0  );
	centre[ J_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationCentreY", 0.0  );
	centre[ K_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationCentreZ", 0.0  );
	size             = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"RadiusCylinder", 0.0  );
	omega            = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationOmega", 1.0  );

	/* Find vector from centre to node */
	StGermain_VectorSubtraction( vector, coord, centre, 2 );

	/*if (context->currentTime > 1.33e-6)
	  omega=0.0;*/
	
	if ((vector[ I_AXIS ]*vector[ I_AXIS ]+vector[ J_AXIS ]*vector[ J_AXIS ])<=size*size)
		*result = -omega * vector[ J_AXIS ];
	else
		*result = 0.0;
}

void StgFEM_StandardConditionFunctions_PartialRotationY(const double *coord, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	Coord                   centre;
	Coord                   vector;
	double                  omega;
	double			size;

	/* Find Centre of Solid Body Rotation */
	centre[ I_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationCentreX", 0.0  );
	centre[ J_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationCentreY", 0.0  );
	centre[ K_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationCentreZ", 0.0  );
	size             = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"RadiusCylinder", 0.0  );
	omega            = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationOmega", 1.0  );

	/* Find vector from centre to node */
	StGermain_VectorSubtraction( vector, coord, centre, 2 );
	
	if ((vector[ I_AXIS ]*vector[ I_AXIS ]+vector[ J_AXIS ]*vector[ J_AXIS ])<=size*size)
		*result =  omega * vector[ I_AXIS ];
	else 
		*result = 0.0;
}


void StgFEM_StandardConditionFunctions_TaperedRotationX(const double *coord, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	Coord                   centre;
	Coord                   vector;
	double                  omega;
	double			size, r, taper;

	/* Find Centre of Solid Body Rotation */
	centre[ I_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationCentreX", 0.0 );
	centre[ J_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationCentreY", 0.0 );
	centre[ K_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationCentreZ", 0.0 );
	size             = Dictionary_GetDouble_WithDefault( dictionary, "RadiusCylinder", 0.0 );
	omega            = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationOmega",   1.0 );

	taper            = Dictionary_GetDouble_WithDefault( dictionary, "TaperedRadius",   0.0 );

	/* Find vector from centre to node */
	StGermain_VectorSubtraction( vector, coord, centre, 2 );

        r=sqrt(vector[ I_AXIS ]*vector[ I_AXIS ]
               +vector[ J_AXIS ]*vector[ J_AXIS ]);
	if (r<=size)
          *result = -omega * vector[ J_AXIS ];
	else if(r<=taper)
          *result = -omega * vector[ J_AXIS ]*(taper-r)/(taper-size);
        else
          *result = 0;
}

void StgFEM_StandardConditionFunctions_TaperedRotationY(const double *coord, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	Coord                   centre;
	Coord                   vector;
	double                  omega;
	double			size, r, taper;

	/* Find Centre of Solid Body Rotation */
	centre[ I_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationCentreX", 0.0 );
	centre[ J_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationCentreY", 0.0 );
	centre[ K_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationCentreZ", 0.0 );
	size             = Dictionary_GetDouble_WithDefault( dictionary, "RadiusCylinder", 0.0 );
	omega            = Dictionary_GetDouble_WithDefault( dictionary, "SolidBodyRotationOmega",   1.0 );

	taper            = Dictionary_GetDouble_WithDefault( dictionary, "TaperedRadius",   0.0 );

	/* Find vector from centre to node */
	StGermain_VectorSubtraction( vector, coord, centre, 2 );

        r=sqrt(vector[ I_AXIS ]*vector[ I_AXIS ]
               +vector[ J_AXIS ]*vector[ J_AXIS ]);
	if (r<=size)
          *result = omega * vector[ I_AXIS ];
	else if(r<=taper)
          *result = omega * vector[ I_AXIS ]*(taper-r)/(taper-size);
        else
          *result = 0;
}




void StgFEM_StandardConditionFunctions_SimpleShear(const double *coord, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double                  centre;
	double                  factor;
	
	/* Find Centre of Solid Body Rotation */
	centre = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SimpleShearCentreY", 0.0  );
	factor = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SimpleShearFactor", 1.0  );

	*result = factor * (coord[ J_AXIS ] - centre);
}

void StgFEM_StandardConditionFunctions_ShearZ(const double *coord, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double                  centre;
	double                  factor;
	
	/* Find Centre of Solid Body Rotation */
	centre = Dictionary_GetDouble_WithDefault( dictionary, "ShearZCentre", 0.0 );
	factor = Dictionary_GetDouble_WithDefault( dictionary, "ShearZFactor", 1.0 );

	*result = factor * (coord[ K_AXIS ] - centre);
}

void StgFEM_StandardConditionFunctions_SimpleShearInverted(const double *coord, void* _context, void* _result ) {
        DomainContext*  context            = (DomainContext*)_context;
        Dictionary*             dictionary         = context->dictionary;
        double*                 result             = (double*) _result;
        double                  factor;

        /* Find Centre of Solid Body Rotation */
        factor = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SimpleShearFactor", 1.0  );

        *result = factor * ( 1.0 - coord[ J_AXIS ] ) ;
}


void StgFEM_StandardConditionFunctions_Extension(const double *coord, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double                  centre;
	double                  factor;
	
	/* Find Centre of Solid Body Rotation */
	centre = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"ExtensionCentreX", 0.0  );
	factor = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"ExtensionFactor", 1.0  );

	*result = factor * (coord[ I_AXIS ] - centre);
}


void StgFEM_StandardConditionFunctions_PartialLid_TopLayer(const double *coord, void* _context, void* result ) {
	DomainContext*	context = (DomainContext*)_context;
	double*			velResult = (double*)result;
	double                  margin = 0;
	double			min[3], max[3];
	FeVariable*             velVar = NULL;
	FeMesh*			mesh = NULL;
	
	velVar = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh = velVar->feMesh;
	Mesh_GetMinimumSeparation( mesh, &margin, NULL );
	Mesh_GetGlobalCoordRange( mesh, min, max );
	margin *= 1.1;
	if((coord[I_AXIS] < (max[I_AXIS] - margin )) && 
           (coord[I_AXIS] > (min[I_AXIS] + margin )))
	{
		(*velResult) = 1;
	}
	else {
		(*velResult) = 0;
	}
}

void StgFEM_StandardConditionFunctions_LinearInterpolationLid(const double *coord, void* _context, void* result ) {
	DomainContext*	context = (DomainContext*)_context;
	FeVariable*             velVar = NULL;
	FeMesh*			mesh = NULL;
	double*			velResult = (double*)result;
	double			boxLength = 0;
	double			leftHandSideValue = 0;
	double			rightHandSideValue = 0;
	double			gradient = 0;
	double			min[3], max[3];
	
	velVar = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh = velVar->feMesh;

	Mesh_GetGlobalCoordRange( mesh, min, max );
	boxLength = max[I_AXIS] - min[I_AXIS];
	leftHandSideValue = Dictionary_GetDouble_WithDefault( context->dictionary, (Dictionary_Entry_Key)"bcLeftHandSideValue", 0.0  );
	rightHandSideValue = Dictionary_GetDouble_WithDefault( context->dictionary, (Dictionary_Entry_Key)"bcRightHandSideValue", 1.0 );
	gradient = (rightHandSideValue - leftHandSideValue) / boxLength;
	(*velResult ) = leftHandSideValue + gradient * (coord[I_AXIS] - min[I_AXIS] );
}


void StgFEM_StandardConditionFunctions_Lid_RampWithCentralMax(const double *coord, void* _context, void* result ) {
	DomainContext*	context = (DomainContext*)_context;
	FeVariable*             velVar = NULL;
	FeMesh*			mesh = NULL;
	double*			velResult = (double*)result;
	double			boxLength = 0;
	double			xPosRelativeToTopLeft = 0;
	double			min[3], max[3];
	
	velVar = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh = velVar->feMesh;

	Mesh_GetGlobalCoordRange( mesh, min, max );
	xPosRelativeToTopLeft = coord[I_AXIS] - min[I_AXIS];
	boxLength = max[I_AXIS] - min[I_AXIS];
	if ( xPosRelativeToTopLeft < boxLength / 2 ) {
		(*velResult) =  2 * xPosRelativeToTopLeft / boxLength;
	}
	else {
		(*velResult) = 1 - 2 * ( xPosRelativeToTopLeft - (boxLength/2) );
	}
}
void StgFEM_StandardConditionFunctions_LinearVelocityLeftWall(const double *coord, void* _context, void* result ) {
	DomainContext*	context = (DomainContext*)_context;
	FeVariable*             velVar = NULL;
	FeMesh*			mesh = NULL;
	double*			velResult = (double*)result;
	Dictionary*             dictionary         = context->dictionary;
	double			min[3], max[3];
	double			gradient, maxvel;
	velVar = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh = velVar->feMesh;

	maxvel = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"MaximumVelocity_Left", 0.0  );
	Mesh_GetGlobalCoordRange( mesh, min, max );
	gradient = maxvel/(min[1] - max[1]);
	
	(*velResult)   = gradient*coord[J_AXIS];
}
void StgFEM_StandardConditionFunctions_LinearVelocityRightWall(const double *coord, void* _context, void* result ) {
	DomainContext*	context = (DomainContext*)_context;
	FeVariable*             velVar = NULL;
	FeMesh*			mesh = NULL;
	double*			velResult = (double*)result;
	Dictionary*             dictionary         = context->dictionary;
	double			min[3], max[3];
	double			gradient, maxvel;
	velVar = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh = velVar->feMesh;

	maxvel = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"MaximumVelocity_Right", 0.0  );
	Mesh_GetGlobalCoordRange( mesh, min, max );
	gradient = maxvel/(max[1] - min[1]);
	 
	(*velResult)   = maxvel - gradient*coord[J_AXIS];
}


void StgFEM_StandardConditionFunctions_SinusoidalLid(const double *coord, void* _context, void* result ) {
	DomainContext*	context = (DomainContext*)_context;
	FeVariable*             velVar = NULL;
	FeMesh*			mesh = NULL;
	double*			velResult = (double*)result;
	double			boxLength = 0;
	double			linearInterp = 0;
	double          	wavenumber;
	double			min[3], max[3];

	wavenumber = Dictionary_GetDouble_WithDefault( context->dictionary, (Dictionary_Entry_Key)"sinusoidalLidWavenumber", 1 );
	
	velVar = (FeVariable* )FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh = velVar->feMesh;

	Mesh_GetGlobalCoordRange( mesh, min, max );
	boxLength = max[I_AXIS] - min[I_AXIS];
	linearInterp = (coord[I_AXIS] - min[I_AXIS] ) / boxLength;
	(*velResult) = sin( linearInterp * M_PI * wavenumber );
}

double StGermain_CosineHillValue( const double* centre, const double* position, const double height, const double diameterAtBase, Dimension_Index dim ) {
	double distanceFromCentre = StGermain_DistanceBetweenPoints( centre, position, dim );
	
	if (distanceFromCentre < diameterAtBase * 0.5 ) 
		return height * (0.5 + 0.5 * cos( 2.0 * M_PI/diameterAtBase * distanceFromCentre ) );
	else
		return 0.0;
}

void StgFEM_StandardConditionFunctions_TemperatureCosineHill(const double *coord, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	Coord                   centre;
	Coord                   rotationCentre;
	double                  omega;
	double                  hillHeight;
	double                  hillDiameter;
	
	/* Read values from dictionary */
	hillHeight       = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"CosineHillHeight"  , 1.0  );
	hillDiameter     = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"CosineHillDiameter", 1.0  );
	centre[ I_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"CosineHillCentreX" , 0.0  );
	centre[ J_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"CosineHillCentreY" , 0.0  );
	centre[ K_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"CosineHillCentreZ" , 0.0  );

	if ( Dictionary_GetBool( dictionary, "RotateCosineHill" ) ) {
		/* Assume solid body rotation */
		rotationCentre[ I_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationCentreX", 0.0  );
		rotationCentre[ J_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationCentreY", 0.0  );
		rotationCentre[ K_AXIS ] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationCentreZ", 0.0  );
		omega                    = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SolidBodyRotationOmega", 1.0  );

		StGermain_VectorSubtraction( centre, rotationCentre, centre, context->dim );
		StGermain_RotateCoordinateAxis( centre, centre, K_AXIS, omega * context->currentTime );
		StGermain_VectorAddition( centre, centre, rotationCentre, context->dim );
	}

	*result = StGermain_CosineHillValue( centre, coord, hillHeight, hillDiameter, context->dim );
}


void StgFEM_StandardConditionFunctions_LinearWithSinusoidalPerturbation(const double *coord, void* _context, void* _result ) {
	DomainContext*	context = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	FeVariable*             feVariable = NULL;
	FeMesh*			feMesh = NULL;
	double*                 result = (double*) _result;
	double                  topLayerBC;
	double                  bottomLayerBC;
	double                  perturbationAmplitude;
	double                  horizontalWaveNumber;
	double                  verticalWaveNumber;
	double                  scaleFactor;
	Coord                   relScaledCoord; 
	double			min[3], max[3], topLayerCoord, bottomLayerCoord;
	
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	feMesh       = feVariable->feMesh;

	Mesh_GetGlobalCoordRange( feMesh, min, max );

	topLayerCoord = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalTempIC_TopLayerCoord", max[J_AXIS]  );
	bottomLayerCoord = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalTempIC_BottomLayerCoord", min[J_AXIS]  );

	topLayerBC = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalTempIC_TopLayerBC", 0.0  );
	bottomLayerBC = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalTempIC_BottomLayerBC", 1.0  );
	scaleFactor = bottomLayerBC - topLayerBC;
	perturbationAmplitude = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalTempIC_PerturbationAmplitude", 0.1  );
	/* Note: these are both multiplied by pi, so wavenumber = 1 means the perturbation goes from 0 to pi, which is
	 * half a full sin or cos cycle. Wavenumber = 3 means the range is 0 -> 3pi, or 1 and a half full cycles. */
	horizontalWaveNumber = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalTempIC_HorizontalWaveNumber", 1.0  );
	verticalWaveNumber = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalTempIC_VerticalWaveNumber", 1.0  );

	/* if node is outside IC shape set to 0 temperature */
	if( coord[J_AXIS] > topLayerCoord || coord[J_AXIS] < bottomLayerCoord ) {
		*result = 0; return ;
	}

	/* make coord relative to box bottom left corner, then scale from 0 to 1 between box min & max */
	relScaledCoord[I_AXIS] = (coord[0] - min[0]) / (max[0] - min[0]);
	relScaledCoord[J_AXIS] = (coord[1] - bottomLayerCoord) / (topLayerCoord - bottomLayerCoord);


	/* Note: ok to use the 1.0 below since we've already scaled the coord to somewhere between 0 to 1 */
	*result = topLayerBC + scaleFactor * ( 1.0 - relScaledCoord[ J_AXIS ] )
		+ perturbationAmplitude * ( cos( horizontalWaveNumber * M_PI * coord[ I_AXIS ] )
					    * sin( verticalWaveNumber * M_PI * relScaledCoord[ J_AXIS ] ) );
}

void StgFEM_StandardConditionFunctions_Trigonometry(const double *coord, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	FeVariable*             feVariable         = NULL;
	FeMesh*			feMesh               = NULL;
	double*                 result             = (double*) _result;
	double                  width;
	double			min[3], max[3];

	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	feMesh       = feVariable->feMesh;

	Mesh_GetGlobalCoordRange( feMesh, min, max );

	/* Get Aspect Ratio */
	width  = max[ I_AXIS ] - min[ I_AXIS ];
	
	*result = 1.0 - 0.5 * M_PI * coord[ J_AXIS ] * sin( M_PI * coord[ I_AXIS ]/width );
}

#define SMALL 1.0e-5
void Stg_FEM_VelicTemperatureIC(const double *coord, void* _context, void* _result ) {
	DomainContext*  context            = (DomainContext*)_context;
	FeVariable*             temperatureField   = (FeVariable*) FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	FeMesh*			feMesh               = temperatureField->feMesh;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double                  x; 
	double                  y;
	double                  kx;
	double                  ky;
	int                     wavenumberX;
	double                  wavenumberY;
	double                  sigma;
	double                  Lx;
	double			min[3], max[3];
	
	Mesh_GetGlobalCoordRange( feMesh, min, max );

	/* Make sure that the box has right dimensions */
	assert( ( max[ J_AXIS ] - min[ J_AXIS ] - 1.0 ) < SMALL );
	Lx = max[ I_AXIS ] - min[ I_AXIS ];

	x = coord[ I_AXIS ] - min[ I_AXIS ];
	y = coord[ J_AXIS ] - min[ J_AXIS ];

	wavenumberX = Dictionary_GetInt_WithDefault( dictionary, (Dictionary_Entry_Key)"wavenumberX", 1  );
	wavenumberY = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"wavenumberY", 1.0  );
	sigma = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"sigma", 1.0 );
	
	assert( sigma > 0.0 );
	assert( wavenumberY > 0.0 );
	assert( wavenumberX > 0.0 );
	
	kx = (double)wavenumberX * M_PI / Lx;
	ky = (double)wavenumberY * M_PI;

	*result = sigma * sin( ky * y ) * cos( kx * x  );
}

/* IC from Mirko Velic. This is the IC temperature for his solB, from his Analytic Suite. Added 22-May-2006 */
void Stg_FEM_VelicTemperatureIC_SolB(const double *coord, void* _context, void* _result ) {
	DomainContext*  context            = (DomainContext*)_context;
	FeVariable*             temperatureField   = (FeVariable*) FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	FeMesh*			feMesh               = temperatureField->feMesh;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double                  x; 
	double                  y;
	double                  km; /*  for y-direction */
	double                  kn; /*  for x-direction */
	double                  wavenumberX;
	double                  wavenumberY;
	double                  L;
	double                  sigma;
	double			min[3], max[3];
	
	Mesh_GetGlobalCoordRange( feMesh, min, max );

	/* Make sure that the box has right dimensions */
	assert( (max[ J_AXIS ] - min[ J_AXIS ] - 1.0 ) < SMALL );
	L = max[ I_AXIS ] - min[ I_AXIS ];

	x = coord[ I_AXIS ] - min[ I_AXIS ];
	y = coord[ J_AXIS ] - min[ J_AXIS ];

	wavenumberX = Dictionary_GetInt_WithDefault( dictionary, (Dictionary_Entry_Key)"wavenumberX", 1  );
	wavenumberY = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"wavenumberY", 2.0 );
	assert( wavenumberX != wavenumberY  );
	sigma = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"sigma", 1.0 );

	kn = wavenumberX * M_PI / L;
/* 	 TODO: Re-write Mirko's code and/or Documentation so the input parameters for these ICs are less confusing */
	km = wavenumberY / L;

	*result = sigma * sinh( km * y ) * cos( kn * x  );
}


/* Initial Condition derived from Boundary Layer theory -
   taken from P. E. van Keken, S. D. King, U. R. Schmeling, U. R. Christensen, D. Neumeister, and M.-P. Doin. A comparison of methods for the modeling of thermochemical convection. Journal of Geophysical Research, 102(B10):22477-22496, october 1997. */
void StgFEM_StandardConditionFunctions_AnalyticalTemperatureIC(const double *coord, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	FeVariable*             feVariable         = NULL;
	FeMesh*			feMesh               = NULL;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double                  u0, v0, Q;
	double                  x, y;
	double                  RaT;
	double                  lambda, height, width;
	double                  Tu, Tl, Tr, Ts;
	double			min[3], max[3];

	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	feMesh       = feVariable->feMesh;

	Mesh_GetGlobalCoordRange( feMesh, min, max );

	/* Get Aspect Ratio */
	height = max[ J_AXIS ] - min[ J_AXIS ];
	width  = max[ I_AXIS ] - min[ I_AXIS ];
	lambda = width/height;
	
	x = coord[ I_AXIS ] - min[ I_AXIS ];
	y = coord[ J_AXIS ] - min[ J_AXIS ];
	
	/* Get thermal Rayleigh Number from Dictionary */
	RaT = Dictionary_GetDouble( dictionary, "RaT" );
	
	/* Horizontal fluid velocity at upper boundary & lower boundary - Equation A3 */
	u0 = pow( lambda , 7.0/3.0 )/ pow(1 + lambda*lambda*lambda*lambda, 2.0/3.0) * pow(0.5*RaT/sqrt(M_PI) , 2.0/3.0);

	/* Vertical velocity of the upwelling and downwellings - Modified from Van Keken to match Turcotte and Shubert */
	v0 = u0; /*lambda; */
	
	/* Total rate of heat flow out of the top of the cell per unit distance along the axis of the roll - Equation A3 */
	Q = 2.0 * sqrt(M_1_PI * lambda/u0);
	Tu = 0.5 * boost::math::erf( 0.5 * ( 1 - y ) * sqrt(u0/x) );                                                      /* Equation A2a */
	Tl = 1.0 - 0.5 * boost::math::erf(0.5 * y * sqrt(u0/(lambda-x)));                                                 /* Equation A2b */
	Tr = 0.5 + 0.5*Q/sqrt(M_PI) * sqrt(v0/(y+1)) * exp( -x*x*v0/(4*y+4) );                               /* Equation A2c */
	Ts = 0.5 - 0.5*Q/sqrt(M_PI) * sqrt(v0/(2-y)) * exp( -(lambda - x) * (lambda - x) * v0 / (8 - 4*y) ); /* Equation A2d */

	/* Equation A1 */
	*result = Tu + Tl + Tr + Ts - 1.5;

	/* Crop result */
	if ( *result > 1.0 ) 
		*result = 1.0;
	else if ( *result < 0.0 ) 
		*result = 0.0;
	
}

void StgFEM_StandardConditionFunctions_EdgeDriveConvectionIC(const double *coord, void* _context, void* _result ) 
{        
	DomainContext*  context = (DomainContext*)_context;        
	Dictionary*             dictionary         = context->dictionary;        
	double*                 result = (double*) _result;        
	double                  perturbationAmplitude;        
	double                  thermalAnomalyOffset;        
	
	perturbationAmplitude = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalTempIC_PerturbationAmplitude", 0.1  );        
	thermalAnomalyOffset = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"thermalAnomalyOffset", 0.0  );        
	/* eqn 1 from S.D.King & D.L. Anderson, "Edge-drive convection", EPSL 160 (1998) 289-296 */        
	
	*result = 1.0 + perturbationAmplitude * sin( M_PI * coord[ J_AXIS ] ) * cos( 0.5 * M_PI * ( coord[ I_AXIS ] + thermalAnomalyOffset ) );
}

void StgFEM_StandardConditionFunctions_ThermalEdgeDriveConvectionIC(const double *coord, void* _context, void* _result )
{
        DomainContext*  context = (DomainContext*)_context;
        Dictionary*             dictionary         = context->dictionary;
        double*                 result = (double*) _result;
        int                     dim;
        double                  contStartX, contEndX;
        double                  contStartY, contEndY;
        double                  contStartZ, contEndZ;
        double                  interiorTemp;

	dim = Dictionary_GetInt_WithDefault( dictionary, (Dictionary_Entry_Key)"dim", 0.0  );
        contStartX = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"contStartX", 0.0  );
        contEndX = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"contEndX", 0.0  );
        contStartY = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"contStartY", 0.0  );
        contEndY = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"contEndY", 0.0  );
	interiorTemp = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"interiorTemp", 1.0 );
        if ( dim == 3  ) {
                contStartZ = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"contStartZ", 0.0  );
                contEndZ = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"contEndZ", 0.0 );
        }

        if(( coord[I_AXIS] >= contStartX && coord[ I_AXIS ] <= contEndX ) && ( coord[J_AXIS] >= contStartY && coord[ J_AXIS ] <= contEndY )) {
                if ( dim == 3 ) {
                        if ( coord[K_AXIS] >= contStartZ && coord[ K_AXIS ] <= contEndZ  )
                                        *result = 0.0;
                        else
                                        *result = interiorTemp;
                }
        }
        else
                        *result = interiorTemp;
}

void StgFEM_StandardConditionFunctions_SinusoidalExtension(const double *coord, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double                  frequency;
	double                  vel0;
	double                  amplitude;
	double                  phaseShift;

	frequency  = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalExtensionFrequency", 1.0  );
	vel0       = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalExtensionVelocity", 0.0  );
	amplitude  = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalExtensionAmplitude", 0.0  );
	phaseShift = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalExtensionPhaseShift", 0.0 );


	*result = vel0 + amplitude * cos( 2.0 * M_PI * frequency * (context->currentTime + context->dt - phaseShift )  );
}


void StgFEM_StandardConditionFunctions_StepFunction(const double *coord, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double                  lower_offset, upper_offset;
	double                  value, lower_value, upper_value;
	unsigned		dim;

	lower_offset = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionLowerOffset", 0.0 );
	upper_offset = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionUpperOffset", lower_offset );
	value = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionValue", 0.0 );
	dim = Dictionary_GetUnsignedInt_WithDefault( dictionary, "StepFunctionDim", 0 );

        lower_value = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionLowerValue", 0.0 );
        upper_value = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionUpperValue", value );

        if(dim==3)
          {
            dim=0;
            coord=&(context->currentTime);
          }

        if(coord[dim] < lower_offset) {
          *result=lower_value;
        } else if(coord[dim] < upper_offset) {
          *result=lower_value + 
            (upper_value-lower_value)
            *(coord[dim] - lower_offset)/(upper_offset-lower_offset);
        } else {
          *result=upper_value;
        }
}


void StG_FEM_StandardConditionFunctions_StepFunctionProduct1(const double *coord, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double                  start, end;
	double                  value;
	unsigned		dim;

	start = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionProduct1Start", 0.0 );
	end = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionProduct1End", 0.0 );
	value = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionProduct1Value", 0.0 );
	dim = Dictionary_GetUnsignedInt_WithDefault( dictionary, "StepFunctionProduct1Dim", 0 );

        if( coord[dim] > start && coord[dim] < end ) {
          *result = value;
        }
        else {
          *result = 0;
        }
}

void StG_FEM_StandardConditionFunctions_StepFunctionProduct2(const double *coord, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double                  start, end;
	double                  value;
	unsigned		dim;

	start = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionProduct2Start", 0.0 );
	end = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionProduct2End", 0.0 );
	value = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionProduct2Value", 0.0 );
	dim = Dictionary_GetUnsignedInt_WithDefault( dictionary, "StepFunctionProduct2Dim", 0 );

        if( coord[dim] > start && coord[dim] < end ) {
          *result = value;
        }
        else {
          *result = 0;
        }
}


void StG_FEM_StandardConditionFunctions_StepFunctionProduct3(const double *coord, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double                  start, end;
	double                  value;
	unsigned		dim;

	start = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionProduct3Start", 0.0 );
	end = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionProduct3End", 0.0 );
	value = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionProduct3Value", 0.0 );
	dim = Dictionary_GetUnsignedInt_WithDefault( dictionary, "StepFunctionProduct3Dim", 1 );

        if( coord[dim] > start && coord[dim] < end ) {
          *result = value;
        }
        else {
          *result = 0;
        }
}

void StG_FEM_StandardConditionFunctions_StepFunctionProduct4(const double *coord, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double                  start, end;
	double                  value;
	unsigned		dim;

	start = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionProduct4Start", 0.0 );
	end = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionProduct4End", 0.0 );
	value = Dictionary_GetDouble_WithDefault( dictionary, "StepFunctionProduct4Value", 0.0 );
	dim = Dictionary_GetUnsignedInt_WithDefault( dictionary, "StepFunctionProduct4Dim", 1 );

        if( coord[dim] > start && coord[dim] < end ) {
          *result = value;
        }
        else {
          *result = 0;
        }
}

/* A Gaussian GaussianHeight*exp(-((GaussianCenter-x)/GaussianWidth)^2) */

void StG_FEM_StandardConditionFunctions_Gaussian
(const double *coord, void* _context,
  void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double                  center, width, height;
	unsigned		dim;

        center = Dictionary_GetDouble_WithDefault( dictionary,
                                                   "GaussianCenter", 0.0 );
	width = Dictionary_GetDouble_WithDefault( dictionary,
                                                  "GaussianWidth", 1.0 );
	height = Dictionary_GetDouble_WithDefault( dictionary,
                                                   "GaussianHeight", 1.0 );
	dim = Dictionary_GetUnsignedInt_WithDefault( dictionary,
                                                     "GaussianDim", 0 );

        *result=height*exp(-(center-coord[dim])*(center-coord[dim])
                           /(width*width));
}

void StgFEM_StandardConditionFunctions_MovingStepFunction
(const double *coord, void* _ctx, void* _result ) {
   FiniteElementContext* ctx = (FiniteElementContext*)_ctx;
   FeVariable* velField;
   FeMesh* mesh;
   Dictionary* dict = ctx->dictionary;
   double* result = (double*)_result;
   double offsetLower, offsetUpper, left, right;
   double *wallCrd, pos;
   int dim, wallDepth;
   unsigned ijk[3];
   char* movingWall;
   Grid* grid;

   /*
   ** Get the velocity field. */
   velField = (FeVariable*)FieldVariable_Register_GetByName(
      ctx->fieldVariable_Register, "VelocityField" );

   /*
   ** Get the mesh. */
   mesh = velField->feMesh;

   /*
   ** Extract all the parameters we need from the dictionary. */
   offsetLower = Dictionary_GetDouble_WithDefault( dict, (Dictionary_Entry_Key)"MovingStepFunctionOffsetLower", 0.0  );
   offsetUpper = Dictionary_GetDouble_WithDefault( dict, (Dictionary_Entry_Key)"MovingStepFunctionOffsetUpper", 0.0  );
   dim = Dictionary_GetUnsignedInt_WithDefault( dict, "MovingStepFunctionDim", 0 );
   left = Dictionary_GetDouble_WithDefault( dict, (Dictionary_Entry_Key)"MovingStepFunctionLeftSide", 0.0  );
   right = Dictionary_GetDouble_WithDefault( dict, (Dictionary_Entry_Key)"MovingStepFunctionRightSide", 0.0  );
   movingWall = Dictionary_GetString_WithDefault( dict, "MovingStepFunctionMovingWall", "lower" );
   wallDepth = Dictionary_GetInt_WithDefault( dict, (Dictionary_Entry_Key)"MovingStepFunctionWallDepth", 0  );

   /*
   ** Because we're dealing with a moving step function, we need to calculate
   ** from where the offset should be applied. */
   grid = *(Grid**)Mesh_GetExtension( mesh, Grid**, "vertexGrid" );
   assert( grid );
   memset( ijk, 0, 3 * sizeof(unsigned) );
   if( !strcmp( movingWall, "lower" ) ) {
      ijk[dim] = wallDepth;
      wallCrd = Mesh_GetVertex( mesh, Grid_Project( grid, ijk ) );
      offsetLower += wallCrd[dim];
      offsetUpper += wallCrd[dim];
   }
   else {
      ijk[dim] = grid->sizes[dim] - wallDepth - 1;
      wallCrd = Mesh_GetVertex( mesh, Grid_Project( grid, ijk ) );
      offsetLower += wallCrd[dim];
      offsetUpper += wallCrd[dim];
   }

   /*
   ** Apply the set of parameters to this node. */
   pos = coord[dim];
   if( pos <= offsetLower )
      *result = left;
   else if( pos >= offsetUpper )
      *result = right;
   else {
      *result = left + ((pos - offsetLower) / (offsetUpper - offsetLower)) * (right - left);
   }
}

void StgFEM_StandardConditionFunctions_ConvectionBenchmark(const double *coord, void* _context, void* _result ) {
	/* This IC is for the 2D ConvectionBenchmark defined in
	 * http://www.mcc.monash.edu.au/twiki/view/Research/ConvectionBenchmarks
	 */
	
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	FeVariable*             feVariable         = NULL;
	FeMesh*			mesh;
	double*                 result             = (double*) _result;
	double			min[3], max[3];
	double                  x,y;
	double                  Lx, Ly;

	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	mesh       = (FeMesh*)feVariable->feMesh;

	Mesh_GetGlobalCoordRange( mesh, min, max );
	
	Lx = max[ I_AXIS ] - min[ I_AXIS ];
	Ly = max[ J_AXIS ] - min[ J_AXIS ];
	
	x = ( coord[0] - min[ I_AXIS ] ) / Lx;
	y = ( coord[1] - min[ J_AXIS ] ) / Ly;


	*result = ( 1 - y ) + ( cos( M_PI * x ) * sin( M_PI * y ) ) / 100 ;
}

void StgFEM_StandardConditionFunctions_ConstantVector(const double *coord, void* _context, void* _result ) {
	DomainContext*		context            = (DomainContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	
	result[0] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"ConstantValueX", 0.0  );
	result[1] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"ConstantValueY", 0.0 );
  if (context->dim == 3  ) 
    result[2] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"ConstantValueZ", 0.0 );
}

/* 3D spec ridge top BC (for milestone 1 of magma project ) 
 * to be applied to the top x-z plane of the domain */
void StgFEM_StandardConditionFunctions_SpecRidge3D(const double *coord, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;

	double			leftVal;
	double			rightVal;
	double			xOffset1;
	double			xOffset2;
	double			yOffset1, yOffset2;
	double			xBegin, xEnd;
	double			zBegin, zEnd;


	leftVal = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SpecRidge3DLeftSide", 0.0  );
	rightVal = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SpecRidge3DRightSide", 0.0  );
	xOffset1 = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SpecRidge3DXOffset1", 0.0  );
	xOffset2 = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SpecRidge3DXOffset2", 0.0  );
	yOffset1 = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SpecRidge3DZOffset1", 0.0  );
	yOffset2 = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SpecRidge3DZOffset2", 0.0  );
	xBegin = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SpecRidge3DXBegin", 0.0  );
	xEnd = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SpecRidge3DXEnd", 0.0  );
	zBegin = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SpecRidge3DZBegin", 0.0  );
	zEnd = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"SpecRidge3DZEnd", 0.0 );

	if( coord[0] < xBegin || coord[0] > xEnd ||
	    coord[2] < zBegin || coord[2] > zEnd )
	{
		*result = 0.0;
	}
	else if( coord[0] < xOffset1 )
		*result = leftVal;
	else if( coord[0] < xOffset2 && coord[2] > yOffset1 && coord[2] < yOffset2  )
		*result = leftVal;
	else
		*result = rightVal;
}

void StgFEM_StandardConditionFunctions_TemperatureProfile(const double *coord, void* _context, void* _result ) {
  FiniteElementContext *	context            = (FiniteElementContext*)_context;
  Dictionary*             dictionary         = context->dictionary;
  double*                 result             = (double*) _result;
  double                  T_0, H_0, dH, H, H_m, A, B, C, x_min, x_max, y_max, T_m, xc, dum;
  /* G.Ito 10/08 added variables x_min, x_max, T_m, Xc, to do variation in x
     and limit maximum T */
  
  T_0 = Dictionary_GetDouble_WithDefault( dictionary, "TemperatureProfileTop", 0.0 );
  T_m = Dictionary_GetDouble_WithDefault( dictionary, "TemperatureProfileMax", 10000.0 );
  H_0 = Dictionary_GetDouble_WithDefault( dictionary, "TemperatureProfileH0", -1.0 );
  H_m = Dictionary_GetDouble_WithDefault( dictionary, "TemperatureProfileHm", 1.0e+8 );
  dH = Dictionary_GetDouble_WithDefault( dictionary, "TemperatureProfiledH", 0.0 );     
  A = Dictionary_GetDouble_WithDefault( dictionary, "TemperatureProfileLinearCoefficient", 0.0 );
  B = Dictionary_GetDouble_WithDefault( dictionary, "TemperatureProfileExponentialCoefficient1", 0.0 );
  C = Dictionary_GetDouble_WithDefault( dictionary, "TemperatureProfileExponentialCoefficient2", 0.0 );
  y_max = Dictionary_GetDouble_WithDefault( dictionary, "maxY", 0.0 );
  x_max = Dictionary_GetDouble_WithDefault( dictionary, "maxX", 0.0 );
  x_min = Dictionary_GetDouble_WithDefault( dictionary, "minX", 0.0 );
  xc = Dictionary_GetDouble_WithDefault( dictionary, "ExtensionCentreX", 0.0 );
  
  if (H_0<0.0)
    {
      if(coord[1]>y_max)
        {
          *result=T_0;
        }
      else
        {
          *result=T_0 + A*(y_max-coord[1]) + B*(1-exp(-C*(y_max-coord[1])));
        }
    }
  else
    {
      if(coord[1]>=y_max)
        {
          *result=T_0;
        }
      else
        {
          H=H_0 + 2*fabs(coord[0]-xc)/(x_max-x_min)*dH;
          if (H>H_m) H=H_m;
          
          dum=T_0 + ((T_m-T_0)/H)*(y_max-coord[1])
            + B*(1-exp(-C*(y_max-coord[1])));
          if (dum>T_m) dum=T_m;
          *result=dum;
        }
    }

}

void StgFEM_StandardConditionFunctions_ERF(const double *coord, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double                  width, scale, dilate, offset, constant;
	unsigned		dim;

	width = Dictionary_GetDouble_WithDefault( dictionary, "ERFWidth", 0.0 );
	offset= Dictionary_GetDouble_WithDefault(dictionary, "ERFOffset",0.0 );
	constant=Dictionary_GetDouble_WithDefault(dictionary,"ERFConstant",0.0);
        scale = Dictionary_GetDouble_WithDefault( dictionary, "ERFScale", 1.0 );
	dilate = Dictionary_GetDouble_WithDefault( dictionary,"ERFDilate",1.0 );
	dim = Dictionary_GetUnsignedInt_WithDefault( dictionary, "ERFDim", 0 );

        if(dim==3)
          {
            dim=0;
            coord=&(context->currentTime);
          }

        if(coord[dim]+offset < -width && width!=0)
          *result=constant-scale;
        else if(coord[dim]+offset > width && width!=0)
          *result=constant+scale;
        else
          *result=constant+scale*boost::math::erf((coord[dim]+offset)/dilate);
}

void StgFEM_StandardConditionFunctions_ERFC(const double *coord,
                                            void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double                  width, scale, dilate, offset, constant;
	unsigned		dim;

	width = Dictionary_GetDouble_WithDefault(dictionary, "ERFCWidth", 0.0 );
	offset= Dictionary_GetDouble_WithDefault(dictionary, "ERFCOffset",0.0 );
	constant=Dictionary_GetDouble_WithDefault(dictionary,"ERFCConstant",0.0);
        scale = Dictionary_GetDouble_WithDefault(dictionary, "ERFCScale", 1.0 );
	dilate = Dictionary_GetDouble_WithDefault(dictionary,"ERFCDilate",1.0 );
	dim = Dictionary_GetUnsignedInt_WithDefault(dictionary, "ERFCDim", 0 );

        if(dim==3)
          {
            dim=0;
            coord=&(context->currentTime);
          }

        if(coord[dim]+offset < -width && width!=0)
          *result=constant-scale;
        else if(coord[dim]+offset > width && width!=0)
          *result=constant+scale;
        else
          *result=constant+scale*boost::math::erfc((coord[dim]+offset)/dilate);
}

void StgFEM_StandardConditionFunctions_RubberSheet(const double *coord,
                                                   void* _context, void* _result)
{
  FiniteElementContext *	context            = (FiniteElementContext*)_context;
  FeVariable*             feVariable         = NULL;
  FeMesh*			feMesh               = NULL;
  Dictionary*             dictionary         = context->dictionary;
  double*                 result             = (double*) _result;
  double                  lower_offset, upper_offset;
  double                  lower_value, upper_value, time;
  unsigned		dim;

  feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
  feMesh       = feVariable->feMesh;

  lower_offset = Dictionary_GetDouble_WithDefault( dictionary,
                                                   "RubberSheetLowerOffset",
                                                   0.0 );
  upper_offset = Dictionary_GetDouble_WithDefault( dictionary,
                                                   "RubberSheetUpperOffset",
                                                   lower_offset );
  dim = Dictionary_GetUnsignedInt_WithDefault( dictionary,
                                               "RubberSheetDim", 0 );

  lower_value = Dictionary_GetDouble_WithDefault( dictionary,
                                                  "RubberSheetLowerValue",
                                                  0.0 );
  upper_value = Dictionary_GetDouble_WithDefault( dictionary,
                                                  "RubberSheetUpperValue",
                                                  0.0 );

  time=context->currentTime;

  if(coord[dim] < lower_offset + lower_value*time)
    {
      *result=lower_value;
    }
  else if(coord[dim] < upper_offset + upper_value*time)
    {
      double min[3], max[3];
      Mesh_GetGlobalCoordRange( feMesh, min, max );
      *result=lower_value + 
        (upper_value-lower_value)
        *(coord[dim] - min[dim])/(max[dim]-min[dim]);
    }
  else
    {
      *result=upper_value;
    }
}

/* error function for use in 3D spec ridge top BC */
double errorFunction( double z, int n ) {
	double		pi	= 3.1415926535;
	double 		a;
	double 		erf 	= 0.0;
	int		denom;
	int		i, j;

	a = 2.0/sqrt( pi );

	for( i=0 ; i<n ; i++ ) {
		denom = 1;
		for( j=1 ; j<=2*i+1 ; j+=2 ) 
			denom *= j; 
		
		erf += pow( 2.0, i )*pow( z, 2*i+1 )/denom;
	}

	return erf *= a*exp( -1.0*z*z );
}



void StgFEM_StandardConditionFunctions_ErrorFunc(const double *coord, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double			dilate;
	double			width; 

	dilate      = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"ErrorFuncDilate", 0.0  );
	width       = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"ErrorFuncWidth", 0.0 );

	if( coord[0] < -1.0*width ) {
		*result = -1.0;
	}
	else if( coord[0] > width  ) {
		*result = 1.0;
	}
	else {
		*result     = errorFunction( coord[0]/dilate, 5 );	
	}
}

void StgFEM_StandardConditionFunctions_GaussianDistribution(const double *coord, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	unsigned		nDims              = context->dim;
	unsigned		dim_I;
	double			orig[3];
	double			sigma              = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"sigma", 1.0  );
	double			gaussianScale      = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"GaussianScale", 1.0  );
	double			background         = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"backgroundValue", 1.0  );
	double			distsq             = 0.0;

	orig[0] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"x0", 0.0  );
	orig[1] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"y0", 0.0  );
	orig[2] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"z0", 0.0 );

	for( dim_I = 0; dim_I < nDims; dim_I++ )
		distsq += ( coord[dim_I] - orig[dim_I] ) * ( coord[dim_I] - orig[dim_I] );

	*result = gaussianScale * exp( -distsq / ( 2.0 * sigma * sigma )  ) + background;
}

void StgFEM_StandardConditionFunctions_GravitationalPotential(const double *coord, void* _context, void* _result ) {
	double* result=(double*) _result;

	*result = -1.0 * coord[J_AXIS];
}

void StgFEM_StandardConditionFunctions_1DGaussianDistribution(const double *coord, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double			orig[3];
	double			sigma              = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"sigma", 1.0  );
	double			gaussianScale      = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"GaussianScale", 1.0  );
	double			background         = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"backgroundValue", 1.0  );
	double			distsq             = 0.0;

	orig[0] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"x0", 0.0  );
	orig[1] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"y0", 0.0  );
	orig[2] = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"z0", 0.0 );

	distsq = ( coord[J_AXIS] - orig[J_AXIS] ) * ( coord[J_AXIS] - orig[J_AXIS] );

	*result = gaussianScale * exp( -distsq / ( 2.0 * sigma * sigma )  ) + background;
}

void StgFEM_StandardConditionFunctions_HalfContainer(const double *coord, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double			halfPoint          = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"halfPoint", 0.0  );

	if( coord[1] < halfPoint )
		*result = 1;
	else
		*result = 0;	
}

void StgFEM_StandardConditionFunctions_ConstantValue(const double *coord, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double			value              = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"constantValue", 1.0  );

	*result = value;
}

void StgFEM_StandardConditionFunctions_DiagonalLine(const double *coord, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double			width              = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"lineWidth", 1.0  );

	if( fabs( coord[0] - coord[1] ) < width )
		*result = 1.0;
	else
		*result = 0.0;
}

void StgFEM_StandardConditionFunctions_DeltaFunction(const double *coord, void* _context, void* _result ) {
	FiniteElementContext *	context            = (FiniteElementContext*)_context;
	Dictionary*             dictionary         = context->dictionary;
	double*                 result             = (double*) _result;
	double			epsilon		   = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"deltaFunctionEpsilon", 0.001  );
	unsigned		dim		   = Dictionary_GetUnsignedInt_WithDefault( dictionary, "deltaFunctionDim", 0 );
	double			centre		   = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"deltaFunctionCentre", 0.5  );
	double			value		   = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"deltaFunctionValue", 1.0  );
	
	*result = (fabs( coord[dim] - centre ) < epsilon) ? value : 0.0;
}

void StgFEM_StandardConditionFunctions_InflowBottom(const double *coord, void* _context, void* _result ) {
	DomainContext*	context            = (DomainContext*)_context;
	FeVariable* feVariable;
	Dictionary*             dictionary         = context->dictionary;
	FeMesh*			mesh               = NULL;
	double*                 result             = (double*) _result;
	double sideLength, wallLength, sideV;
	double min[3], max[3];
	
	feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	mesh       = feVariable->feMesh;
	Mesh_GetGlobalCoordRange( mesh, min, max );
	sideLength = max[1] - min[1];
	wallLength = max[0] - min[0];
	sideV = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"InflowSideVelocity", 1.0  );

	*result = 2.0 * sideV * sideLength / wallLength;
}

void StgFEM_StandardConditionFunctions_GaussianTube(const double *coord, void* _context, void* _result ) 
{
        DomainContext*  context = (DomainContext*)_context;
        Dictionary*             dictionary         = context->dictionary;
        FeVariable*             feVariable = NULL;
        FeMesh*                 feMesh = NULL;
        unsigned                nDims;
        double*                 result = (double*) _result;
        double                  a1,b1,c1, a2,b2,c2, x,y,z,r_y,r_yz;
        double                  min[3], max[3];
	double                  y_shift, z_shift;

        feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
        feMesh       = feVariable->feMesh;

        nDims = Mesh_GetDimSize( feMesh );
        Mesh_GetGlobalCoordRange( feMesh, min, max );

 	a1 = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"GaussianTube_a1", 1.0  ); /* Scales the magnitude of the perturbation. */
	c1 = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"GaussianTube_c1", 0.1  ); /* Controls the smoothing length. Smaller values produce less smoothing. */

        a2 = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"GaussianTube_a2", 0.05  ); /* Controls ampltude of oscillations */
        b2 = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"GaussianTube_b2", 6.28318530718  ); /* Controls frequency of oscillations */
	c2 = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"GaussianTube_c2", 1.570796326795  ); /* Shifts oscillations */

	y_shift = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"GaussianTube_y_origin", 0.0  );
	z_shift = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"GaussianTube_z_origin", 0.0  );


	x = coord[ I_AXIS ];
	y = coord[ J_AXIS ];

	y = y - y_shift;
	if (nDims==2) {
		b1 = a2 * sin( b2*x - c2 );
		r_y = sqrt( (y-b1)*(y-b1) );
		*result = a1 * exp( -(r_y * r_y) / (2.0*c1*c1) );
	}
	if (nDims==3) {
		z = coord[ K_AXIS ];
		z = z - z_shift;

		b1 = a2 * sin( b2*x - c2 );
		r_yz = sqrt( (y-b1)*(y-b1) + z*z );
		*result = a1 * exp( -(r_yz * r_yz) / (2.0*c1*c1) );
	}


}


void StgFEM_StandardConditionFunctions_WarsTemperature(const double *coord, void* _context, void* _result ) 
{
  FiniteElementContext *	context            = (FiniteElementContext*)_context;
  Dictionary*             dictionary         = context->dictionary;
  double*                 result             = (double*) _result;
  double                  EAEnd, WarsStart, WarsHeight, WarsTTop,
    WarsTBottom, h, maxY;
  
  EAEnd = Dictionary_GetDouble( dictionary, "EAEnd");
  WarsStart = Dictionary_GetDouble( dictionary, "WarsStart");
  WarsHeight = Dictionary_GetDouble( dictionary, "WarsHeight");
  WarsTTop = Dictionary_GetDouble( dictionary, "WarsTTop");
  WarsTBottom = Dictionary_GetDouble( dictionary, "WarsTBottom");     
  maxY=Dictionary_GetDouble( dictionary, "maxY");


  h=WarsHeight*(coord[0]-EAEnd)/(WarsStart-EAEnd);
  if(coord[0]<EAEnd)
    h=0;
  if(coord[0]>WarsStart)
    h=WarsHeight;
  *result=WarsTBottom + ((coord[1]-h)/(maxY-h))*(WarsTTop-WarsTBottom);
}

void StgFEM_StandardConditionFunctions_Quadratic(const double *coord, void* _context, void* _result ) 
{
  FiniteElementContext *	context            = (FiniteElementContext*)_context;
  Dictionary*             dictionary         = context->dictionary;
  double*                 result             = (double*) _result;
  int                     dim;
  double                  a, b, c;
  
  dim = Dictionary_GetInt( dictionary, "Quadratic_Dim");
  a = Dictionary_GetDouble( dictionary, "Quadratic_Constant");
  b = Dictionary_GetDouble( dictionary, "Quadratic_Linear");
  c = Dictionary_GetDouble( dictionary, "Quadratic_Quadratic");

  *result= a + coord[dim]*(b + c*coord[dim]);
}

int Binary_Search(double *data, int s, int e, double value);

void StgFEM_StandardConditionFunctions_FileN(const double *coord, void* _context, void* _result, int file_num, double **coords, double **data);

void StgFEM_StandardConditionFunctions_File0(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,0,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File1(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,1,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File2(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,2,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File3(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,3,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File4(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,4,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File5(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,5,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File6(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,6,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File7(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,7,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File8(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,8,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File9(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,9,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File10(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,10,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File11(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,11,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File12(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,12,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File13(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,13,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File14(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,14,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File15(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,15,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File16(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,16,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File17(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,17,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File18(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,18,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File19(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,19,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File20(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,10,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File21(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,11,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File22(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,12,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File23(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,13,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File24(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,14,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File25(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,15,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File26(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,16,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File27(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,17,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File28(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,18,&coords,&data);
}

void StgFEM_StandardConditionFunctions_File29(const double *coord, void* _context, void* _result ) 
{
  static double *coords=NULL;
  static double *data=NULL;
  StgFEM_StandardConditionFunctions_FileN(coord,_context,_result,19,&coords,&data);
}

void StgFEM_StandardConditionFunctions_FileN(const double *coord, void* _context, void* _result, int file_num, double **coords, double **data)
{
  FiniteElementContext *	context            = (FiniteElementContext*)_context;
  Dictionary*             dictionary         = context->dictionary;
  double*                 result             = (double*) _result;
  int                     dim, dim2, dim3(-1), i, j, k;
  char *filename;
  int N, N2, N3, ndims;
  double factor, factor2, factor3;
  
  char fileN_number[10], fileN_dim[15], fileN_dim2[15], fileN_dim3[15],
    fileN_name[128], fileN_N[15], fileN_N2[15], fileN_N3[15];
  sprintf(fileN_number,"File%d",file_num);
  sprintf(fileN_dim,"File%d_Dim",file_num);
  sprintf(fileN_dim2,"File%d_Dim2",file_num);
  sprintf(fileN_dim3,"File%d_Dim3",file_num);
  sprintf(fileN_name,"File%d_Name",file_num);
  sprintf(fileN_N,"File%d_N",file_num);
  sprintf(fileN_N2,"File%d_N2",file_num);
  sprintf(fileN_N3,"File%d_N3",file_num);

  filename = Dictionary_GetString( dictionary, fileN_name);
  N = Dictionary_GetInt( dictionary, fileN_N);
  dim = Dictionary_GetInt( dictionary, fileN_dim);
  N2 = Dictionary_GetInt_WithDefault( dictionary, fileN_N2,-1);
  N3 = Dictionary_GetInt_WithDefault( dictionary, fileN_N3,-1);

  if(N2==-1)
    ndims=1;
  else if(N3==-1)
    ndims=2;
  else
    ndims=3;

  if(ndims>1)
    {
      dim2 = Dictionary_GetInt( dictionary, fileN_dim2);
      if(ndims>2)
        dim3 = Dictionary_GetInt( dictionary, fileN_dim3);
    }

  Journal_Firewall(dim>=0 && dim<3,
                   Journal_Register( Error_Type,"StgFEM_StandardConditionFunctions_FileN"),
                   "%s must be either 0, 1, or 2, but was set to %d\n",
                   fileN_dim,dim);
  Journal_Firewall(N>0,
                   Journal_Register( Error_Type,"StgFEM_StandardConditionFunctions_FileN"),
                   "%s must be greater than zero, but was set to %d.\n",
                   fileN_N,N);
  if(*data==NULL)
    {
      FILE *fp=fopen(filename,"r");
      Journal_Firewall(fp!=NULL,
                       Journal_Register( Error_Type,"StgFEM_StandardConditionFunctions_FileN"),
                       "Bad filename for %s.  Could not open %s\n",
                       fileN_name,filename);

      /* In 1D, data and coords are simple 1D arrays.  In 2D, data is
         a 2D arrays, and coord is still a 1D array, with the first N
         elements being the coordinates in the Dim direction, and the
         next N2 elements being the coordinates in the Dim2
         dirction.  Similarly in 3D.  */
      if(ndims==1)
        {
          *data=(double *)malloc(N*sizeof(double));
          *coords=(double *)malloc(N*sizeof(double));
        }
      else if(ndims==2)
        {
          *data=(double *)malloc(N*N2*sizeof(double));
          *coords=(double *)malloc((N+N2)*sizeof(double));
        }
      else
        {
          *data=(double *)malloc(N*N2*N3*sizeof(double));
          *coords=(double *)malloc((N+N2+N3)*sizeof(double));
        }

      Journal_Firewall(*data!=NULL && *coords!=NULL,
                       Journal_Register( Error_Type,"StgFEM_StandardConditionFunctions_FileN"),
                       "Could not allocate enough memory for %s\n",file_num);
      if(ndims==1)
        {
          for(i=0;i<N;++i)
            fscanf(fp,"%lf %lf",*coords+i,*data+i);
        }
      else if(ndims==2)
        {
          for(i=0;i<N;++i)
            for(j=0;j<N2;++j)
              {
                fscanf(fp,"%lf %lf %lf",*coords+i,*coords+N+j,*data+i+N*j);
              }
        }
      else if(ndims==3)
        {
          for(i=0;i<N;++i)
            for(j=0;j<N2;++j)
              for(k=0;k<N3;++k)
                fscanf(fp,"%lf %lf %lf %lf",*coords+i,*coords+N+j,*coords+N+N2+k,
                       *data+i+N*(j+N2*k));
        }
      fclose(fp);
    }

  Journal_Firewall(!(coord[dim]<(*coords)[0] || coord[dim]>(*coords)[N-1]),
                   Journal_Register( Error_Type,"StgFEM_StandardConditionFunctions_FileN"),
                   "The range in the file '%s' does not cover this value %g\nIt only covers %g to %g in the %d direction.\n",
                   filename,coord[dim],(*coords)[0],(*coords)[N-1],dim);
  if(ndims>1)
    Journal_Firewall(!(coord[dim2]<(*coords)[N] || coord[dim2]>(*coords)[N+N2-1]),
                     Journal_Register( Error_Type,"StgFEM_StandardConditionFunctions_FileN"),
                     "The range in the file '%s' does not cover this value %g\nIt only covers %g to %g in the %d direction.\n",
                     filename,coord[dim2],(*coords)[N],(*coords)[N+N2-1],dim2);
  if(ndims>2)
    Journal_Firewall(!(coord[dim3]<(*coords)[N+N2] || coord[dim3]>(*coords)[N+N2+N3-1]),
                     Journal_Register( Error_Type,"StgFEM_StandardConditionFunctions_FileN"),
                     "The range in the file '%s' does not cover this value %g\nIt only covers %g to %g in the %d direction.\n",
                     filename,coord[dim3],(*coords)[N+N2],(*coords)[N+N2+N3-1],dim2);

  // i=std::lower_bound(coords,coords+N,coord[dim])-coords;
  i=Binary_Search(*coords,0,N-1,coord[dim]);
  factor=((*coords)[i+1]-coord[dim])/((*coords)[i+1]-(*coords)[i]);
    
  if(ndims>1)
    {
      j=Binary_Search(*coords,N,N+N2-1,coord[dim2]);
      factor2=((*coords)[j+1]-coord[dim2])/((*coords)[j+1]-(*coords)[j]);
      j=j-N;
      if(ndims>2)
        {
          k=Binary_Search(*coords,N,N+N2+N3-1,coord[dim3]);
          factor3=((*coords)[k+1]-coord[dim3])/((*coords)[k+1]-(*coords)[k]);
          k=k-N-N2;
        }
    }

  switch(ndims)
    {
    case 1:
      *result=(*data)[i]*factor + (*data)[i+1]*(1-factor);
      break;
    case 2:
      *result=factor*(factor2*(*data)[i+N*j] + (1-factor2)*(*data)[i+N*(j+1)])
        + (1-factor)*(factor2*(*data)[i+1+N*j] + (1-factor2)*(*data)[i+1+N*(j+1)]);
      break;
    case 3:
      *result=factor*(factor2*(factor3*(*data)[i+N*(j+N2*k)]
                               + (1-factor3)*(*data)[i+N*(j+N2*(k+1))])
                      + (1-factor2)*(factor3*(*data)[i+N*((j+1)+N2*k)]
                                     + (1-factor3)*(*data)[i+N*((j+1)+N2*(k+1))]))
        + (1-factor)*(factor2*(factor3*(*data)[i+1+N*(j+N2*k)]
                               + (1-factor3)*(*data)[i+1+N*(j+N2*(k+1))])
                      + (1-factor2)*(factor3*(*data)[i+1+N*((j+1)+N2*k)]
                                     + (1-factor3)*(*data)[i+1+N*((j+1)+N2*(k+1))]));
      break;
    }
}
    

int Binary_Search(double *data, int s, int e, const double value)
{
  int start, end, midpoint;

  start=s;
  end=e;
  midpoint=(end-start)/2 + start;
  while(start!=midpoint)
    {
      if(data[midpoint]>=value)
        end=midpoint;
      else
        start=midpoint;
      midpoint=(end-start)/2 + start;
    }
  return start;
}
