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
** $Id: testStiffnessMatrixCompare.c 964 2007-10-11 08:03:06Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/Discretisation/Discretisation.h>
#include <StgFEM/SLE/SLE.h>
#include <assert.h>

const Type testStiffnessMatrixCompare_Type = "testStiffnessMatrixCompare";
typedef struct {
	__Codelet
} testStiffnessMatrixCompare;

void testStiffnessMatrixCompareFunction( FiniteElementContext* context ) {
	StiffnessMatrix* matrix;
	Matrix*          savedMatrix;
	Name             matrixName;
	Dictionary*      dictionary     = context->dictionary;
	double           tolerance;
	double           errorNorm;
	double           matrixNorm;
	Stream*          stream;
	Name             filename;
	Stream*          infoStream     = Journal_Register( Info_Type, CURR_MODULE_NAME );

	/* Get Matrix */
	matrixName = Dictionary_GetString( dictionary, "CompareStiffnessMatrix" );
	Journal_Printf( infoStream, "Comparing stiffness matrix '%s'\n", matrixName );
	matrix = (StiffnessMatrix*) LiveComponentRegister_Get( context->CF->LCRegister, matrixName );
	assert( matrix );

	/* Get Stored Matrix from file */
	filename = Dictionary_GetString( dictionary, "StiffnessMatrixCompareFilename" );
	Journal_Printf( infoStream, "Checking with file '%s'\n", filename );
	if ( Dictionary_GetBool( dictionary, "patchtests" ) ) {
		Matrix_Dump( matrix->matrix, filename );
		Journal_Printf( infoStream, "Patching file '%s'\n", filename );
	}

#ifndef HAVE_PETSC
#error No PETSc!
#endif
	savedMatrix = (Matrix*)PETScMatrix_New( "" );
	Matrix_Load( savedMatrix, filename );

	/* Get Error */
	matrixNorm = Matrix_L2Norm( savedMatrix );
	Matrix_AddScaled( savedMatrix, -1.0, matrix->matrix );
	errorNorm = Matrix_L2Norm( savedMatrix );

	Journal_PrintValue( infoStream, matrixNorm );
	Journal_PrintValue( infoStream, errorNorm );

	/* Check tolerance */
	stream = Journal_Register( Info_Type, "StiffnessMatrixComparison" );
	Stream_RedirectFile_WithPrependedPath( stream, context->outputPath, "StiffnessMatrixCompare.dat" );
	tolerance = Dictionary_GetDouble( dictionary, "StiffnessMatrixCompareTolerance" );
	Journal_PrintValue( infoStream, tolerance );
	Journal_Printf( stream, "Comparison between stiffness matrix '%s' %s with tolerance %4g.\n", 
			matrixName, 
			( errorNorm/matrixNorm < tolerance ? "passed" : "failed" ),
			tolerance );

	/* Clean Up */
	FreeObject( savedMatrix );
	Stream_CloseFile( stream );
}

void _testStiffnessMatrixCompare_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	AbstractContext* context = (AbstractContext*)Stg_ComponentFactory_ConstructByName( 
		cf, 
		"context",
		AbstractContext, 
		True,
		data );

	ContextEP_Append( context, AbstractContext_EP_Solve, testStiffnessMatrixCompareFunction );
}

void* _testStiffnessMatrixCompare_DefaultNew( Name name ) {
	return Codelet_New(
			testStiffnessMatrixCompare_Type,
			_testStiffnessMatrixCompare_DefaultNew,
			_testStiffnessMatrixCompare_Construct,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
}

Index testStiffnessMatrixCompare_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, testStiffnessMatrixCompare_Type, "0", _testStiffnessMatrixCompare_DefaultNew );
}





