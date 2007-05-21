/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
** Copyright (c) 2005-2006, Monash Cluster Computing, Building 28, Monash University Clayton Campus,
**	Victoria, 3800, Australia
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**
** Contributors:
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Patrick D. Sunter, Software Engineer, VPAC. (patrick@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	David Stegman, Postdoctoral Fellow, Monash University. (david.stegman@sci.monash.edu.au)
**	Wendy Sharples, PhD Student, Monash University (wendy.sharples@sci.monash.edu.au)
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
** $Id: testElementIntegral.c 464 2007-05-21 04:11:40Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/Voronoi/Voronoi.h>
#include <PICellerator/Weights/Weights.h>

#include <math.h>
#include <string.h>
#include <assert.h>

/* forward declaration. Implementation exists in MaterialPoints module */
extern void _IntegrationPointsSwarm_UpdateHook( void* timeIntegrator, void* swarm );


const Type TestElementIntegral_Type = "TestElementIntegral";

void ConstantFunction( void* feVariable, Element_DomainIndex dElement_I, Coord xi, double* value ) {
	*value = -3.0;
}
void LinearFunction( void* feVariable, Element_DomainIndex dElement_I, Coord xi, double* value ) {
	double x = xi[0];
	double y = xi[1];

	*value = 2.0 + 2.2 * x - y;
}
void ShapeFunction( void* feVariable, Element_DomainIndex dElement_I, Coord xi, double* value ) {
	double x = xi[0];
	double y = xi[1];

	*value = 1 + x + y + x * y;
}
void PolynomialFunction( void* feVariable, Element_DomainIndex dElement_I, Coord xi, double* value ) {
	double x = xi[0];
	double y = xi[1];

	*value = 11 + 2*x*x + 3*x*x*x*y + y + x*x*x + 2*y*y;
}
void QuadraticFunction( void* feVariable, Element_DomainIndex dElement_I, Coord xi, double* value ) {
	double x = xi[0];

	*value = 1 + x + x * x;
}
void ExponentialFunction( void* feVariable, Element_DomainIndex dElement_I, Coord xi, double* value ) {
	double x = xi[0];
	double y = xi[1];

	*value = 5*exp(2*x*x*x + 2*y*y*y) * (1-x) * (1+y);
}
void ExponentialInterface( void* feVariable, Element_DomainIndex dElement_I, Coord xi, double* value ) {
	double x = xi[0];
	double y = xi[1];

	*value = (double) (y <= 0.1 * exp( 2*x ));
}
void CircleInterface( void* feVariable, Element_DomainIndex dElement_I, Coord xi, double* value ) {
	double x = xi[0];
	double y = xi[1];

	*value = (double) (x*x + y*y <= 1.0);
}


void testElementIntegral_CompareAgainstReferenceSolution(DiscretisationContext* context, Stream* stream, double mean, double standardDeviation) {
	double meanTolerance, stdDevTolerance;
	double expectedMean, expectedStdDev;
	double differenceMean, differenceStdDev;	
	/* Get the tolerance */
	meanTolerance = Dictionary_GetDouble_WithDefault( context->dictionary, "mean-tolerance", 0.005 );
	stdDevTolerance = Dictionary_GetDouble_WithDefault( context->dictionary, "standardDeviation-tolerance", 0.005 );
	
	/* Get the expected values */
	expectedMean = Dictionary_GetDouble_WithDefault(context->dictionary, "mean-expectedValue", 0.5);
	expectedStdDev = Dictionary_GetDouble_WithDefault(context->dictionary, "standardDeviation-expectedValue", 0.5);
	
	/* compare the values */
	differenceMean = fabs(mean - expectedMean);
	differenceStdDev = fabs(standardDeviation - expectedStdDev);
	
	Stream_RedirectFile_WithPrependedPath( stream, context->outputPath, "output.dat" );
	
	if (differenceMean <= meanTolerance ) {
		Journal_Printf(stream, "Mean is within tolerance, %f of %f\n",
			meanTolerance, expectedMean);
	}
	else {
		Journal_Printf(stream, "Mean is not within tolerance, %f of %f\n",
			meanTolerance, expectedMean);
		Journal_Printf(stream, "	value = ( %f ) \n", mean );
	}
	if ( differenceStdDev <= stdDevTolerance ) {
		Journal_Printf(stream, "Standard Deviation is within tolerance, %f of %f\n",
			stdDevTolerance, expectedStdDev);
	}
	else {
		Journal_Printf(stream, "Standard Deviation is not within tolerance, %f of %f\n",
			stdDevTolerance, expectedStdDev);
		Journal_Printf(stream, "	value = ( %f ) \n", standardDeviation);
	}
	
	
}


void PICellerator_testElementIntegral( DiscretisationContext* context ) {
	Swarm*              integrationSwarm = (Swarm*)LiveComponentRegister_Get( context->CF->LCRegister, "swarm" );
	Swarm*              materialSwarm    = (Swarm*)LiveComponentRegister_Get( context->CF->LCRegister, "picIntegrationPoints" );
	FeMesh* mesh             = (FeMesh*) LiveComponentRegister_Get( context->CF->LCRegister, "mesh-linear" );
	WeightsCalculator*  weights          = (WeightsCalculator*) LiveComponentRegister_Get( context->CF->LCRegister, "weights" );
	FeVariable*         feVariable;
	Element_LocalIndex  lElement_I       = 0;
	double              analyticValue    = 0.0;
	double              integral         = 0.0;
	double              error;
	double              errorSquaredSum  = 0.0;
	double              errorSum         = 0.0;
	/*Particle_Index      lParticle_I;*/
	/*IntegrationPoint*   particle;*/
	double              mean;
	double              standardDeviation;
	Index               loop_I;
	Index               count            = Dictionary_GetUnsignedInt_WithDefault( context->dictionary, "SampleSize", 5000 );
	Name                funcName;
	Stream*             stream           = Journal_Register( Info_Type, CURR_MODULE_NAME );
	/*ElementLayout*      elementLayout    = mesh->layout->elementLayout;*/
	void*               data;

	assert( integrationSwarm );
	assert( materialSwarm );
	assert( mesh );
	assert( weights );

	/* Create FeVariable */
	feVariable = FeVariable_New_Full( "feVariable", mesh, NULL, NULL, NULL, NULL, NULL, NULL, 
			1, context->dim, False, StgFEM_Native_ImportExportType, StgFEM_Native_ImportExportType,
			MPI_COMM_WORLD, context->fieldVariable_Register );

	funcName = Dictionary_GetString( context->dictionary, "FunctionName" );
	if ( strcasecmp( funcName, "ShapeFunction" ) == 0 ) {
		feVariable->_interpolateWithinElement = ShapeFunction;
		analyticValue = 4.0;
	}
	else if ( strcasecmp( funcName, "ConstantFunction" ) == 0 ) {
		feVariable->_interpolateWithinElement = ConstantFunction;
		analyticValue = -12.0;
	}
	else if ( strcasecmp( funcName, "LinearFunction" ) == 0 ) {
		feVariable->_interpolateWithinElement = LinearFunction;
		analyticValue = 8.0;
	}
	else if ( strcasecmp( funcName, "QuadraticFunction" ) == 0 ) {
		feVariable->_interpolateWithinElement = QuadraticFunction;
		analyticValue = 16.0/3.0;
	}
	else if ( strcasecmp( funcName, "PolynomialFunction" ) == 0 ) {
		feVariable->_interpolateWithinElement = PolynomialFunction;
		analyticValue = 148.0/3.0;
	}
	else if ( strcasecmp( funcName, "ExponentialFunction" ) == 0 ) {
		feVariable->_interpolateWithinElement = ExponentialFunction;
		analyticValue = 0.0 /*TODO*/;
		abort();
	}
	else if ( strcasecmp( funcName, "ExponentialInterface" ) == 0 ) {
		feVariable->_interpolateWithinElement = ExponentialInterface;
		analyticValue = 0.05 * (exp(2) - exp(-2)) + 2.0;
	}
	else if ( strcasecmp( funcName, "CircleInterface" ) == 0 ) {
		feVariable->_interpolateWithinElement = CircleInterface;
		analyticValue = M_PI;
	}
	else 
		Journal_Firewall( False,
				Journal_Register( Error_Type, CURR_MODULE_NAME ),
				"Cannot understand function name '%s'\n", funcName );

	for ( loop_I = 0 ; loop_I < count ; loop_I++ ) {
		/* Layout Particles */
		Swarm_Random_Seed( (long) loop_I );
		_Swarm_InitialiseParticles( materialSwarm );

		_IntegrationPointsSwarm_UpdateHook( NULL, integrationSwarm );
		
		WeightsCalculator_CalculateCell( weights, integrationSwarm, lElement_I );

		/* Evaluate Integral */
		integral = FeVariable_IntegrateElement( feVariable, integrationSwarm, lElement_I );

		/* Calculate Error */
		error = fabs( integral - analyticValue )/fabs( analyticValue );
		errorSum += error;
		errorSquaredSum += error*error;
	}

	/* Calculate Mean and Standard Deviation */
	mean = errorSum / (double) count;
	standardDeviation = sqrt( errorSquaredSum / (double) count - mean * mean );

	/* compare mean and standardDeviation against ref. solution
		to a tolerance taken from the xml files. */
	testElementIntegral_CompareAgainstReferenceSolution(context, stream, mean, standardDeviation );

	
	/* Write to file */
/*	Stream_RedirectFile_WithPrependedPath( stream, context->outputPath, "output.dat" );
	/Journal_Printf( stream, "%u \t %.5g %.5g\n", integrationSwarm->cellParticleCountTbl[ lElement_I ], mean, standardDeviation );
*/	
}
	
	
void _testElementIntegral_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	DiscretisationContext* context;
	context = (DiscretisationContext*)Stg_ComponentFactory_ConstructByName( cf, "context", DiscretisationContext, True, data ); 

	ContextEP_ReplaceAll( context, AbstractContext_EP_Execute, PICellerator_testElementIntegral );
}


void* _testElementIntegral_DefaultNew( Name name ) {
	return _Codelet_New(
			sizeof( Codelet ),
			TestElementIntegral_Type,
			_Codelet_Delete,
			_Codelet_Print,
			_Codelet_Copy,
			_testElementIntegral_DefaultNew,
			_testElementIntegral_Construct,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
}


Index testElementIntegral_Register( PluginsManager* pluginsManager ) {
	Index result;

	result = PluginsManager_Submit( pluginsManager, TestElementIntegral_Type, "0",
		_testElementIntegral_DefaultNew );

	return result;
}
