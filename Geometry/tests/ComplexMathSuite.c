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
** Role:
**   Tests the ComplexMathSuite
**
** $Id: testTemplate.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h> 
#include "StgDomain/Geometry/Geometry.h"
#include "StgDomain/Shape/Shape.h"
#include "StgDomain/Mesh/Mesh.h" 
#include "StgDomain/Utils/Utils.h"
#include "StgDomain/Swarm/Swarm.h"

#include "ComplexMathSuite.h"

#undef Journal_PrintCmplx
/** Print a complex number.  Will use %.7f formatting */
#define Journal_PrintCmplx( stream, self ) \
                Journal_Printf( stream, #self " = %.7f %c %.7f i\n", (self)[ REAL_PART ], (self)[ IMAG_PART ] >= 0.0 ? '+' : '-', fabs( (self)[ IMAG_PART ] ) )

#define TESTCOMPLEXMATH_TOL 1e-15

typedef struct {
	MPI_Comm	comm;
	unsigned rank;
	unsigned nProcs;
} ComplexMathSuiteData;

void ComplexMathSuite_Setup( ComplexMathSuiteData* data ) {
	/* MPI Initializations */
	data->comm = MPI_COMM_WORLD;
	MPI_Comm_rank( data->comm, &data->rank );
	MPI_Comm_size( data->comm, &data->nProcs );
}

void ComplexMathSuite_Teardown( ComplexMathSuiteData* data ) {
}

void ComplexMathSuite_TestComplexJournalPrintingMacro( ComplexMathSuiteData* data ) {
	unsigned	procToWatch = data->nProcs >=2 ? 1 : 0;
	
	if (data->rank == procToWatch) {
		Cmplx		x = {1, 2};
		Cmplx		y = {-1.5, 3};
		Cmplx		u = {-1.5, -3};
		Cmplx		v = {1.5, -3};
		Cmplx		i = {0, 1};
		char		expected_file[PCU_PATH_MAX];
		Stream*	stream = Journal_Register( InfoStream_Type, "ComplexMathPrintingMacro" );
		Stream_RedirectFile( stream, "testComplexMathPrintingMacro.dat" );

		Journal_Printf(stream, "----------------- Testing Complex Journal Printing Macro -----------------\n" );
		Journal_PrintCmplx( stream, x );	
		Journal_PrintCmplx( stream, y );	
		Journal_PrintCmplx( stream, u );	
		Journal_PrintCmplx( stream, v );	
		Journal_PrintCmplx( stream, i );	

		pcu_filename_expected( "testComplexJournalPrintingMacro.expected", expected_file );
		pcu_check_fileEq( "testComplexMathPrintingMacro.dat", expected_file );
		remove( "testComplexMathPrintingMacro.dat" );

		Stream_CloseAndFreeFile( stream );
	}
}

void ComplexMathSuite_TestAddition( ComplexMathSuiteData* data ) {
	unsigned	procToWatch = data->nProcs >=2 ? 1 : 0;
	
	if (data->rank == procToWatch) {
		Cmplx		x = {1, 2};
		Cmplx		y = {-1.5, 3};
		Cmplx		dest = {0, 0};
		char		expected_file[PCU_PATH_MAX];
		Stream*	stream = Journal_Register( InfoStream_Type, "ComplexMathAddition" );
		Stream_RedirectFile( stream, "testComplexMathAddition.dat" );

		Journal_Printf(stream, "----------------- Testing Addition -----------------\n" );
		Cmplx_Add( x, y, dest );
		Cmplx_Add( x, y, x );
		Journal_PrintCmplx( stream, dest );	
		Journal_PrintCmplx( stream, x );	

		pcu_filename_expected( "testComplexMathAddition.expected", expected_file );
		pcu_check_fileEq( "testComplexMathAddition.dat", expected_file );
		remove( "testComplexMathAddition.dat" );

		Stream_CloseAndFreeFile( stream );
	}
}

void ComplexMathSuite_TestSubtraction( ComplexMathSuiteData* data ) {
	unsigned	procToWatch = data->nProcs >=2 ? 1 : 0;
	
	if (data->rank == procToWatch) {
		Cmplx		x = {1, 2};
		Cmplx		y = {-1.5, 3};
		Cmplx		dest = {0, 0};
		char		expected_file[PCU_PATH_MAX];
		Stream*	stream = Journal_Register( InfoStream_Type, "ComplexMathSubtraction" );
		Stream_RedirectFile( stream, "testComplexMathSubtraction.dat" );

		Journal_Printf(stream, "----------------- Testing Subtraction -----------------\n" );
		Cmplx_Subtract( x, y, x );
		Cmplx_Subtract( x, y, dest );
		Journal_PrintCmplx( stream, x );
		Journal_PrintCmplx( stream, dest );

		pcu_filename_expected( "testComplexMathSubtraction.expected", expected_file );
		pcu_check_fileEq( "testComplexMathSubtraction.dat", expected_file );
		remove( "testComplexMathSubtraction.dat" );

		Stream_CloseAndFreeFile( stream );
	}
}

void ComplexMathSuite_TestMultiplication( ComplexMathSuiteData* data ) {
	unsigned	procToWatch = data->nProcs >=2 ? 1 : 0;
	
	if (data->rank == procToWatch) {
		Cmplx		x = {1, 2};
		Cmplx		y = {-1.5, 3};
		Cmplx		dest = {0, 0};
		char		expected_file[PCU_PATH_MAX];
		Stream*	stream = Journal_Register( InfoStream_Type, "ComplexMathMultiplication" );
		Stream_RedirectFile( stream, "testComplexMathMultiplication.dat" );

		Journal_Printf(stream, "----------------- Testing Multiplication -----------------\n" );
		Cmplx_Multiply( x, y, dest );
		Cmplx_Multiply( x, y, y );
		Journal_PrintCmplx( stream, dest );	
		Journal_PrintCmplx( stream, y );	

		pcu_filename_expected( "testComplexMathMultiplication.expected", expected_file );
		pcu_check_fileEq( "testComplexMathMultiplication.dat", expected_file );
		remove( "testComplexMathMultiplication.dat" );

		Stream_CloseAndFreeFile( stream );
	}
}

void ComplexMathSuite_TestDivision( ComplexMathSuiteData* data ) {
	unsigned	procToWatch = data->nProcs >=2 ? 1 : 0;
	
	if (data->rank == procToWatch) {
		Cmplx		x = {1, 2};
		Cmplx		y = {-1.5, 3};
		Cmplx		dest = {0, 0};
		char		expected_file[PCU_PATH_MAX];
		Stream*	stream = Journal_Register( InfoStream_Type, "ComplexMathDivision" );
		Stream_RedirectFile( stream, "testComplexMathDivision.dat" );

		Journal_Printf(stream, "----------------- Testing Division -----------------\n" );
		Cmplx_Division( y, x, y );
		Cmplx_Division( x, y, dest );
		Journal_PrintCmplx( stream, y );	
		Journal_PrintCmplx( stream, dest );	

		pcu_filename_expected( "testComplexMathDivision.expected", expected_file );
		pcu_check_fileEq( "testComplexMathDivision.dat", expected_file );
		remove( "testComplexMathDivision.dat" );

		Stream_CloseAndFreeFile( stream );
	}
}

void ComplexMathSuite_TestRealNumber( ComplexMathSuiteData* data ) {
	unsigned	procToWatch = data->nProcs >=2 ? 1 : 0;
	
	if (data->rank == procToWatch) {
		Cmplx		x = {1, 2};
		Cmplx		y = {-1.5, 3};
		Cmplx		dest = {0, 0};
		char		expected_file[PCU_PATH_MAX];
		Stream*	stream = Journal_Register( InfoStream_Type, "ComplexMathRealNumber" );
		Stream_RedirectFile( stream, "testComplexMathRealNumber.dat" );

		Journal_Printf(stream, "----------------- Testing Real Number Math stuff -----------------\n" );
		Cmplx_AddReal( y, 2, dest );
		Journal_PrintCmplx( stream, dest );	
		Cmplx_RealMinusCmplx( y, 4, dest );
		Journal_PrintCmplx( stream, dest );	
		Cmplx_RealMultiply( y, 0.1, dest );
		Journal_PrintCmplx( stream, dest );	
		Cmplx_RealDivideByCmplx( y, 1.0, dest );
		Journal_PrintCmplx( stream, dest );	
		Cmplx_Multiply( y, dest, dest );

		if (fabs(dest[IMAG_PART]) <= TESTCOMPLEXMATH_TOL ) {
			Journal_Printf(stream, "Answer within tolerance %g of %.5f + i %.5f\n", TESTCOMPLEXMATH_TOL,
				dest[REAL_PART], fabs(0));
		}
		else{
			Journal_Printf(stream, "Answer not within tolerance %g of %.5f + i %.5f\n", 
				TESTCOMPLEXMATH_TOL, 1, 0);
		}

		pcu_filename_expected( "testComplexMathRealNumber.expected", expected_file );
		pcu_check_fileEq( "testComplexMathRealNumber.dat", expected_file );
		remove( "testComplexMathRealNumber.dat" );

		Stream_CloseAndFreeFile( stream );
	}
}

void ComplexMathSuite_TestConjugate( ComplexMathSuiteData* data ) {
	unsigned	procToWatch = data->nProcs >=2 ? 1 : 0;
	
	if (data->rank == procToWatch) {
		Cmplx		x = {1, 2};
		Cmplx		y = {-1.5, 3};
		Cmplx		dest = {0, 0};
		char		expected_file[PCU_PATH_MAX];
		Stream*	stream = Journal_Register( InfoStream_Type, "ComplexMathConjugate" );
		Stream_RedirectFile( stream, "testComplexMathConjugate.dat" );

		Journal_Printf(stream, "----------------- Testing Conjugate -----------------\n" );
		Cmplx_Conjugate( x, x );
		Journal_PrintCmplx( stream, x );	
		Cmplx_Conjugate( x, x );
		Journal_PrintCmplx( stream, x );	

		pcu_filename_expected( "testComplexMathConjugate.expected", expected_file );
		pcu_check_fileEq( "testComplexMathConjugate.dat", expected_file );
		remove( "testComplexMathConjugate.dat" );

		Stream_CloseAndFreeFile( stream );
	}
}

void ComplexMathSuite_TestPolar( ComplexMathSuiteData* data ) {
	unsigned	procToWatch = data->nProcs >=2 ? 1 : 0;
	
	if (data->rank == procToWatch) {
		Cmplx		x = {1, 2};
		Cmplx		y = {-1.5, 3};
		Cmplx		u = {-1.5, -3};
		Cmplx		v = {1.5, -3};
		Cmplx		i = {0, 1};
		Cmplx		minus_i = {0, -1};
		Cmplx		dest = {0, 0};
		double	mod, theta;
		char		expected_file[PCU_PATH_MAX];
		Stream*	stream = Journal_Register( InfoStream_Type, "ComplexMathPolar" );
		Stream_RedirectFile( stream, "testComplexMathPolar.dat" );

		Journal_Printf(stream, "----------------- Testing Complex Polar Stuff -----------------\n" );
		mod = Cmplx_Modulus(x);
		theta = Cmplx_Argument( x );
		Journal_Printf(stream, "x = %2.4lf e^{i %2.4lf} = %2.4lf + %2.4lf\n", mod, theta, mod * cos(theta), mod *sin(theta));
		mod = Cmplx_Modulus(y);
		theta = Cmplx_Argument( y );
		Journal_Printf(stream, "y = %2.4lf e^{i %2.4lf} = %2.4lf + %2.4lf\n", mod, theta, mod * cos(theta), mod *sin(theta));
		mod = Cmplx_Modulus(u);
		theta = Cmplx_Argument( u );
		Journal_Printf(stream, "u = %2.4lf e^{i %2.4lf} = %2.4lf + %2.4lf\n", mod, theta, mod * cos(theta), mod *sin(theta));
		mod = Cmplx_Modulus(v);
		theta = Cmplx_Argument( v );
		Journal_Printf(stream, "v = %2.4lf e^{i %2.4lf} = %2.4lf + %2.4lf\n", mod, theta, mod * cos(theta), mod *sin(theta));
		mod = Cmplx_Modulus(i);
		theta = Cmplx_Argument( i );
		Journal_Printf(stream, "i = %2.4lf e^{i %2.4lf} = %2.4lf + %2.4lf\n", mod, theta, mod * cos(theta), mod *sin(theta));
		mod = Cmplx_Modulus(minus_i);
		theta = Cmplx_Argument( minus_i );
		Journal_Printf(stream, "-i = %2.4lf e^{i %2.4lf} = %2.4lf + %2.4lf\n", mod, theta, mod * cos(theta), mod *sin(theta));

		pcu_filename_expected( "testComplexMathPolar.expected", expected_file );
		pcu_check_fileEq( "testComplexMathPolar.dat", expected_file );
		remove( "testComplexMathPolar.dat" );

		Stream_CloseAndFreeFile( stream );
	}
}

void ComplexMathSuite_TestPower( ComplexMathSuiteData* data ) {
	unsigned	procToWatch = data->nProcs >=2 ? 1 : 0;
	
	if (data->rank == procToWatch) {
		Cmplx		x = {1, 2};
		Cmplx		y = {-1.5, 3};
		Cmplx		dest = {0, 0};
		char		expected_file[PCU_PATH_MAX];
		Stream*	stream = Journal_Register( InfoStream_Type, "ComplexMathPower" );
		Stream_RedirectFile( stream, "testComplexMathPower.dat" );

		Journal_Printf(stream, "----------------- Testing Complex to real Power Stuff -----------------\n" );
		Cmplx_RealPower( x, 2.0, dest );
		Journal_PrintCmplx( stream, dest );	
		Cmplx_Multiply( x, x, dest );
		Journal_PrintCmplx( stream, dest );	
		Cmplx_Multiply( dest, x, dest );
		Journal_PrintCmplx( stream, dest );	
		Cmplx_RealPower( x, 3.0, dest );
		Journal_PrintCmplx( stream, dest );	
		Cmplx_Sqrt( y, dest );
		Journal_PrintCmplx( stream, dest );	
		Cmplx_RealPower( dest, 2.0, dest );
		Journal_PrintCmplx( stream, dest );	

		Journal_Printf(stream, "\n----------------- Testing Complex to complex Power Stuff -----------------\n" );
		Cmplx_CmplxPower( x, y, dest );
		Journal_PrintCmplx( stream, dest );	

		pcu_filename_expected( "testComplexMathPower.expected", expected_file );
		pcu_check_fileEq( "testComplexMathPower.dat", expected_file );
		remove( "testComplexMathPower.dat" );

		Stream_CloseAndFreeFile( stream );
	}
}

void ComplexMathSuite_TestBeautifulEquation( ComplexMathSuiteData* data ) {
	unsigned	procToWatch = data->nProcs >=2 ? 1 : 0;
	
	if (data->rank == procToWatch) {
		Cmplx		x = {1, 2};
		Cmplx		y = {-1.5, 3};
		Cmplx		i = {0, 1};
		Cmplx		e = {M_E, 0};
		Cmplx		ipi = {0, M_PI};
		Cmplx		dest = {0, 0};
		char		expected_file[PCU_PATH_MAX];
		Stream*	stream = Journal_Register( InfoStream_Type, "ComplexMathBeautifulEquation" );
		Stream_RedirectFile( stream, "testComplexMathBeautifulEquation.dat" );

		Journal_Printf(stream, "----------------- Testing The Most Beautiful Equation in Mathematics e^{i \\pi} + 1 = 0 -----------------\n" );
		Journal_Printf( stream, "e^{i \\pi} = ");
		Cmplx_CmplxPower( e, ipi, dest );
		Journal_PrintCmplx( stream, dest );	

		Journal_Printf(stream, "\n----------------- Another Beautiful Equation i^i = e^{-\\pi/2} -----------------\n" );
		Cmplx_CmplxPower( i, i, dest );
		Journal_Printf( stream, "e^{-\\pi/2} = %2.5lf = ", exp( -M_PI * 0.5 ) );
		Journal_PrintCmplx( stream, dest );	

		pcu_filename_expected( "testComplexMathBeautifulEquation.expected", expected_file );
		pcu_check_fileEq( "testComplexMathBeautifulEquation.dat", expected_file );
		remove( "testComplexMathBeautifulEquation.dat" );

		Stream_CloseAndFreeFile( stream );
	}
}

void ComplexMathSuite_TestExponential( ComplexMathSuiteData* data ) {
	unsigned	procToWatch = data->nProcs >=2 ? 1 : 0;
	
	if (data->rank == procToWatch) {
		Cmplx		x = {1, 2};
		Cmplx		y = {-1.5, 3};
		Cmplx		e = {M_E, 0};
		Cmplx		ipi = {0, M_PI};
		Cmplx		dest = {0, 0};
		char		expected_file[PCU_PATH_MAX];
		Stream*	stream = Journal_Register( InfoStream_Type, "ComplexMathExponential" );
		Stream_RedirectFile( stream, "testComplexMathExponential.dat" );

		Journal_Printf(stream, "----------------- Testing Exponential -----------------\n" );
		Journal_Printf( stream, "e^{2 + 3i} = ");
		x[ REAL_PART ] = 2.0;
		x[ IMAG_PART ] = 3.0;
		Cmplx_Exp( x, dest );
		Journal_PrintCmplx( stream, dest );	
		Cmplx_CmplxPower( e, x, dest );
		Journal_PrintCmplx( stream, dest );	

		Journal_Printf( stream, "e^{-5 + 7i} = ");
		x[ REAL_PART ] = -5.0;
		x[ IMAG_PART ] = 7.0;
		Cmplx_Exp( x, dest );
		Journal_PrintCmplx( stream, dest );	
		Cmplx_CmplxPower( e, x, dest );
		Journal_PrintCmplx( stream, dest );	

		Journal_Printf( stream, "e^{i \\pi} = ");
		Cmplx_Exp( ipi, dest );
		Journal_PrintCmplx( stream, dest );	

		pcu_filename_expected( "testComplexMathExponential.expected", expected_file );
		pcu_check_fileEq( "testComplexMathExponential.dat", expected_file );
		remove( "testComplexMathExponential.dat" );

		Stream_CloseAndFreeFile( stream );
	}
}

void ComplexMathSuite_TestCopyAndZero( ComplexMathSuiteData* data ) {
	unsigned	procToWatch = data->nProcs >=2 ? 1 : 0;
	
	if (data->rank == procToWatch) {
		Cmplx		x = {1, 2};
		Cmplx		y = {-1.5, 3};
		Cmplx		i = {0, 1};
		Cmplx		e = {M_E, 0};
		Cmplx		ipi = {0, M_PI};
		Cmplx		dest = {0, 0};
		char		expected_file[PCU_PATH_MAX];
		Stream*	stream = Journal_Register( InfoStream_Type, "ComplexMathCopyAndZero" );
		Stream_RedirectFile( stream, "testComplexMathCopyAndZero.dat" );

		Journal_Printf(stream, "----------------- Testing Copy and Zero -----------------\n" );
		Journal_PrintCmplx( stream, x );	
		Journal_PrintCmplx( stream, dest );	
		Cmplx_Copy( x, dest );
		Journal_PrintCmplx( stream, dest );	
		Cmplx_Zero( dest );
		Journal_PrintCmplx( stream, dest );	

		pcu_filename_expected( "testComplexMathCopyAndZero.expected", expected_file );
		pcu_check_fileEq( "testComplexMathCopyAndZero.dat", expected_file );
		remove( "testComplexMathCopyAndZero.dat" );

		Stream_CloseAndFreeFile( stream );
	}
}

void ComplexMathSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, ComplexMathSuiteData );
   pcu_suite_setFixtures( suite, ComplexMathSuite_Setup, ComplexMathSuite_Teardown );
	pcu_suite_addTest( suite, ComplexMathSuite_TestComplexJournalPrintingMacro );
   pcu_suite_addTest( suite, ComplexMathSuite_TestAddition );
   pcu_suite_addTest( suite, ComplexMathSuite_TestSubtraction );
   pcu_suite_addTest( suite, ComplexMathSuite_TestMultiplication );
   pcu_suite_addTest( suite, ComplexMathSuite_TestDivision );
   pcu_suite_addTest( suite, ComplexMathSuite_TestRealNumber );
   pcu_suite_addTest( suite, ComplexMathSuite_TestConjugate );
   pcu_suite_addTest( suite, ComplexMathSuite_TestPolar );
   pcu_suite_addTest( suite, ComplexMathSuite_TestPower );
   pcu_suite_addTest( suite, ComplexMathSuite_TestBeautifulEquation );
   pcu_suite_addTest( suite, ComplexMathSuite_TestExponential );
   pcu_suite_addTest( suite, ComplexMathSuite_TestCopyAndZero );
}
