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
** $Id: LidDrivenStokesAnalytic.c 922 2007-07-25 03:01:19Z DavidLee $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include <assert.h>
#include <math.h>

const Type StgFEM_LidDrivenStokesAnalytic_Type = "StgFEM_LidDrivenStokesAnalytic";

typedef struct {
	__Codelet
	unsigned int n;
	//double A, B, C, D;
	FieldTest* fieldTest;
} StgFEM_LidDrivenStokesAnalytic;

void StgFEM_LidDrivenStokesAnalytic_CalculateConstants( FieldTest *fieldTest ) {
	double                  E;
	double                  e_nPI;
	double                  e_2nPI;
	double                  e_4nPI;
	double                  n;

	unsigned int*		intArray;
	double*			dblArray;

	intArray = (unsigned int*)fieldTest->data[0];
	dblArray = (double*)fieldTest->data[1];
	n = (double) *intArray;

	e_nPI = exp( n * M_PI );
	e_2nPI = e_nPI * e_nPI;
	e_4nPI = e_2nPI * e_2nPI;

	E = (4.0 * n * n * M_PI * M_PI + 2.0 ) * e_2nPI - e_4nPI - 1.0;

	dblArray[0] = ( e_2nPI - 1.0 )* e_nPI / E;
	dblArray[1] = - dblArray[0];

	dblArray[2] =   ( 2.0 * n * M_PI - e_2nPI + 1.0 ) * e_nPI / E;
	dblArray[3] = - ( 2.0 * n * M_PI * e_2nPI - e_2nPI + 1.0 ) * e_nPI / E;
}

void StgFEM_LidDrivenStokesAnalytic_VelocityFunction( void* data, double* coord, double* velocity ) {
	FieldTest*	fieldTest = (FieldTest*) data;
	double x,y;
	double n;
	double A, B, C, D;
	unsigned int*	intArray = (unsigned int*)fieldTest->data[0];
	double*		dblArray = (double*)fieldTest->data[1];

	/* Get local copy of constants */
	n = (double) *intArray;
	A = dblArray[0];
	B = dblArray[1];
	C = dblArray[2];
	D = dblArray[3];

	/* get copy of coords */
	x = coord[I_AXIS];
	y = coord[J_AXIS];
	
	velocity[ I_AXIS ] = sin( n * M_PI * x ) * 
		( ( A * n * M_PI + C + C * n * M_PI * y) *exp( n * M_PI * y ) 
		- ( B * n * M_PI - D + D * n * M_PI * y ) * exp( - n * M_PI * y ) );
	velocity[ J_AXIS ] = - n * M_PI * cos( n * M_PI * x ) * 
		( ( A + C * y ) * exp( n * M_PI * y ) 
		+ ( B + D * y ) * exp( - n * M_PI * y ) );
}


void StgFEM_LidDrivenStokesAnalytic_PressureFunction( void* data, double* coord, double* pressure ) {
	FieldTest*	fieldTest = (FieldTest*) data;
	double x,y;
	double n;
	double A, B, C, D;
	unsigned int*	intArray = (unsigned int*)fieldTest->data[0];
	double*		dblArray = (double*)fieldTest->data[1];

	/* Get local copy of constants */
	n = (double) *intArray;
	A = dblArray[0];
	B = dblArray[1];
	C = dblArray[2];
	D = dblArray[3];

	/* get copy of coords */
	x = coord[I_AXIS];
	y = coord[J_AXIS];
	
	*pressure = - 2.0 * n * M_PI * cos( n * M_PI * x ) * ( C * exp( n * M_PI * y ) + D * exp( - n * M_PI * y ) );
}

void _StgFEM_LidDrivenStokesAnalytic_Construct( void* codelet, Stg_ComponentFactory* cf, void* data ) {
	StgFEM_LidDrivenStokesAnalytic *self = (StgFEM_LidDrivenStokesAnalytic*)codelet;
	
	unsigned int* waveSpeed;

	self->fieldTest = Stg_ComponentFactory_ConstructByName( cf, "ReferenceFields", FieldTest, True, data );

	/* Set constants */
	self->fieldTest->data = Memory_Alloc_Array_Unnamed( void*, 2 );
	self->fieldTest->data[0] = (unsigned int*) Memory_Alloc_Array_Unnamed( unsigned int, 1 );
	self->fieldTest->data[1] = (double*) Memory_Alloc_Array_Unnamed( double, 4 );

	waveSpeed = Memory_Alloc_Array_Unnamed( unsigned int, 1 );
	*waveSpeed = Stg_ComponentFactory_GetRootDictUnsignedInt( cf, "sinusoidalLidWavenumber", 1 );
	self->fieldTest->data[0] = waveSpeed;
	StgFEM_LidDrivenStokesAnalytic_CalculateConstants( self->fieldTest );
}

void _StgFEM_LidDrivenStokesAnalytic_Build( void* codelet, void* data ) {
	StgFEM_LidDrivenStokesAnalytic *self = (StgFEM_LidDrivenStokesAnalytic*)codelet;

	_Codelet_Build( self, data );

	/* build the field test component so that the analytic function array is initialised */
	Stg_Component_Build( self->fieldTest, data, False );

	/* set the analytic solution functions for the feVariables - must be in the same order as the XML */
	FieldTest_AddAnalyticSolutionFuncToListAtIndex( self->fieldTest, 0, StgFEM_LidDrivenStokesAnalytic_VelocityFunction, 0 );
	FieldTest_AddAnalyticSolutionFuncToListAtIndex( self->fieldTest, 1, StgFEM_LidDrivenStokesAnalytic_PressureFunction, 1 );
}

void* _StgFEM_LidDrivenStokesAnalytic_DefaultNew( Name name ) {
	return Codelet_New(
		StgFEM_LidDrivenStokesAnalytic_Type,
		_StgFEM_LidDrivenStokesAnalytic_DefaultNew,
		_StgFEM_LidDrivenStokesAnalytic_Construct,
		_StgFEM_LidDrivenStokesAnalytic_Build,
		_Codelet_Initialise,
		_Codelet_Execute,
		_Codelet_Destroy,
		name );
	/*return (void*) _FieldTest_New( 
			sizeof(StgFEM_LidDrivenStokesAnalytic),
			StgFEM_LidDrivenStokesAnalytic_Type,
			_FieldTest_Delete,
			_FieldTest_Print,
			_FieldTest_Copy,
			_StgFEM_LidDrivenStokesAnalytic_DefaultNew,
			_StgFEM_LidDrivenStokesAnalytic_Construct,
			_StgFEM_LidDrivenStokesAnalytic_Build, 
			_FieldTest_Initialise,
			_FieldTest_Execute,
			_FieldTest_Destroy,
			name );*/
}

Index StgFEM_LidDrivenStokesAnalytic_Register( PluginsManager* pluginsManager ) {
	Journal_DPrintf( StgFEM_Debug, "In: %s( void* )\n", __func__ ); 

	return PluginsManager_Submit( pluginsManager, StgFEM_LidDrivenStokesAnalytic_Type, "0", _StgFEM_LidDrivenStokesAnalytic_DefaultNew );
}


