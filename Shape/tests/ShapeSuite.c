#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h>
#include "StgDomain/StgDomain.h"

typedef struct {
	Stg_ComponentFactory* cf;
} ShapeSuiteData;

void ShapeSuite_Setup( ShapeSuiteData* data ) { 
}

void ShapeSuite_Teardown( ShapeSuiteData* data ) {
   stgMainDestroy( data->cf );
}

void ShapeSuite_GeneratePoints( ShapeSuiteData* data, Dimension_Index dim, char* inputFileName ) {
/** Test Definition: */
	Stg_ComponentFactory* cf;
	DomainContext*   context = NULL;
	Dictionary*      dictionary;
	Stg_Shape*       shape;
	unsigned         testCoordCount, index;
	Name             outputPath;
	Coord            coord;
	Stream*          stream = Journal_Register( Info_Type, inputFileName );
	char xml_input[PCU_PATH_MAX];

	/* read in the xml input file */
	pcu_filename_input( inputFileName, xml_input );
	data->cf = cf = stgMainInitFromXML( xml_input, MPI_COMM_WORLD, NULL );
	stgMainBuildAndInitialise( cf );
	context = (DomainContext*)LiveComponentRegister_Get( cf->LCRegister, "context" ); 


	dictionary = context->dictionary;
	outputPath = Dictionary_GetString( dictionary, "outputPath" );
	Stream_RedirectFile_WithPrependedPath( stream, outputPath, "test.dat" );
	shape = (Stg_Shape*) LiveComponentRegister_Get( context->CF->LCRegister, "shape" );
	assert( shape );

	testCoordCount = Dictionary_GetUnsignedInt_WithDefault( dictionary, "testCoordCount", 10000 );

	/* Test to see if random points are in shape */
	srand48(0);
	for (index = 0 ; index < testCoordCount ; index++ ) {
		coord[ I_AXIS ] = drand48() - 1.0;
		coord[ J_AXIS ] = drand48() - 1.0;
		if ( dim == 3 ) 
			coord[ K_AXIS ] = drand48() - 1.0;

		if ( Stg_Shape_IsCoordInside( shape, coord ) ) 
			Journal_Printf( stream, "%u\n", index );
	}

	Stream_CloseAndFreeFile( stream );
}

void ShapeSuite_TestBox2D( ShapeSuiteData* data ) {
	Dimension_Index  dim = 2;
	char expected_file[PCU_PATH_MAX];
	ShapeSuite_GeneratePoints( data, dim, "testBox2D.xml" );
	
	pcu_filename_expected( "testBox2D.expected", expected_file );
	pcu_check_fileEq( "output/test.dat", expected_file );
	remove("output/test.dat");
}
void ShapeSuite_TestBox3D( ShapeSuiteData* data ) {
	Dimension_Index  dim = 3;
	char expected_file[PCU_PATH_MAX];
	ShapeSuite_GeneratePoints( data, dim, "testBox3D.xml" );

	pcu_filename_expected( "testBox3D.expected", expected_file );
	pcu_check_fileEq( "output/test.dat", expected_file );
	remove("output/test.dat");
}

void ShapeSuite_TestSphere2D( ShapeSuiteData* data ) {
	Dimension_Index  dim = 2;
	char expected_file[PCU_PATH_MAX];
	ShapeSuite_GeneratePoints( data, dim, "testSphere2D.xml" );

	pcu_filename_expected( "testSphere2D.expected", expected_file );
	pcu_check_fileEq( "output/test.dat", expected_file );
	remove("output/test.dat");
}

void ShapeSuite_TestSphere2D_Invert( ShapeSuiteData* data ) {
	Dimension_Index  dim = 2;
	char expected_file[PCU_PATH_MAX];
	ShapeSuite_GeneratePoints( data, dim, "testSphere-invert.xml" );

	pcu_filename_expected( "testSphere2D-invert.expected", expected_file );
	pcu_check_fileEq( "output/test.dat", expected_file );
	remove("output/test.dat");
}

void ShapeSuite_TestSphere3D( ShapeSuiteData* data ) {
	Dimension_Index  dim = 3;
	char expected_file[PCU_PATH_MAX];
	ShapeSuite_GeneratePoints( data, dim, "testSphere3D.xml" );

	pcu_filename_expected( "testSphere3D.expected", expected_file );
	pcu_check_fileEq( "output/test.dat", expected_file );
	remove("output/test.dat");
}

void ShapeSuite_TestConvexHull2D( ShapeSuiteData* data ) {
	Dimension_Index  dim = 2;
	char expected_file[PCU_PATH_MAX];
	ShapeSuite_GeneratePoints( data, dim, "testConvexHull2D.xml" );

	pcu_filename_expected( "testConvexHull2D.expected", expected_file );
	pcu_check_fileEq( "output/test.dat", expected_file );
	remove("output/test.dat");
}

void ShapeSuite_TestConvexHull3D( ShapeSuiteData* data ) {
	Dimension_Index  dim = 3;
	char expected_file[PCU_PATH_MAX];
	ShapeSuite_GeneratePoints( data, dim, "testConvexHull3D.xml" );

	pcu_filename_expected( "testConvexHull3D.expected", expected_file );
	pcu_check_fileEq( "output/test.dat", expected_file );
	remove("output/test.dat");
}

void ShapeSuite_TestUnion2D( ShapeSuiteData* data ) {
	Dimension_Index  dim = 2;
	char expected_file[PCU_PATH_MAX];
	ShapeSuite_GeneratePoints( data, dim, "testUnion2D.xml" );

	pcu_filename_expected( "testUnion2D.expected", expected_file );
	pcu_check_fileEq( "output/test.dat", expected_file );
	remove("output/test.dat");
}

void ShapeSuite_TestUnion3D( ShapeSuiteData* data ) {
	Dimension_Index  dim = 3;
	char expected_file[PCU_PATH_MAX];
	ShapeSuite_GeneratePoints( data, dim, "testUnion3D.xml" );

	pcu_filename_expected( "testUnion3D.expected", expected_file );
	pcu_check_fileEq( "output/test.dat", expected_file );
	remove("output/test.dat");
}

void ShapeSuite_TestUnion2DSingleNOT( ShapeSuiteData* data ) {
	Dimension_Index  dim = 2;
	char expected_file[PCU_PATH_MAX];
	ShapeSuite_GeneratePoints( data, dim, "testUnion2DSingleNot.xml" );

	pcu_filename_expected( "testUnion2DSingleNot.expected", expected_file );
	pcu_check_fileEq( "output/test.dat", expected_file );
	remove("output/test.dat");
}

void ShapeSuite_TestUnion3DSingleNOT( ShapeSuiteData* data ) {
	Dimension_Index  dim = 3;
	char expected_file[PCU_PATH_MAX];
	ShapeSuite_GeneratePoints( data, dim, "testUnion3DSingleNot.xml" );

	pcu_filename_expected( "testUnion3DSingleNot.expected", expected_file );
	pcu_check_fileEq( "output/test.dat", expected_file );
	remove("output/test.dat");
}

void ShapeSuite_TestIntersection2D( ShapeSuiteData* data ) {
	Dimension_Index  dim = 2;
	char expected_file[PCU_PATH_MAX];
	ShapeSuite_GeneratePoints( data, dim, "testIntersection2D.xml" );

	pcu_filename_expected( "testIntersection2D.expected", expected_file );
	pcu_check_fileEq( "output/test.dat", expected_file );
	remove("output/test.dat");
}

void ShapeSuite_TestIntersection3D( ShapeSuiteData* data ) {
	Dimension_Index  dim = 3;
	char expected_file[PCU_PATH_MAX];
	ShapeSuite_GeneratePoints( data, dim, "testIntersection3D.xml" );

	pcu_filename_expected( "testIntersection3D.expected", expected_file );
	pcu_check_fileEq( "output/test.dat", expected_file );
	remove("output/test.dat");
}

void ShapeSuite_TestCylinder( ShapeSuiteData* data ) {
	Dimension_Index  dim = 2;
	char expected_file[PCU_PATH_MAX];
	ShapeSuite_GeneratePoints( data, dim, "testCylinder.xml" );

	pcu_filename_expected( "testCylinder.expected", expected_file );
	pcu_check_fileEq( "output/test.dat", expected_file );
	remove("output/test.dat");
}

void ShapeSuite_TestPolygonShape2D( ShapeSuiteData* data ) {
	Dimension_Index  dim = 2;
	char expected_file[PCU_PATH_MAX];
	ShapeSuite_GeneratePoints( data, dim, "testPolygonShape.xml" );

	pcu_filename_expected( "testPolygonShape2D.expected", expected_file );
	pcu_check_fileEq( "output/test.dat", expected_file );
	remove("output/test.dat");
}

void ShapeSuite_TestSuperellipsoid2D( ShapeSuiteData* data ) {
	Dimension_Index  dim = 2;
	char expected_file[PCU_PATH_MAX];
	ShapeSuite_GeneratePoints( data, dim, "testSuperellipsoid2D.xml" );

	pcu_filename_expected( "testSuperellipsoid2D.expected", expected_file );
	pcu_check_fileEq( "output/test.dat", expected_file );
	remove("output/test.dat");
}

void ShapeSuite_TestSuperellipsoid3D( ShapeSuiteData* data ) {
	Dimension_Index  dim = 3;
	char expected_file[PCU_PATH_MAX];
	ShapeSuite_GeneratePoints( data, dim, "testSuperellipsoid3D.xml" );

	pcu_filename_expected( "testSuperellipsoid3D.expected", expected_file );
	pcu_check_fileEq( "output/test.dat", expected_file );
	remove("output/test.dat");
}

void ShapeSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, ShapeSuiteData );
   pcu_suite_setFixtures( suite, ShapeSuite_Setup, ShapeSuite_Teardown );
   pcu_suite_addTest( suite, ShapeSuite_TestBox2D );
   pcu_suite_addTest( suite, ShapeSuite_TestBox3D );
   pcu_suite_addTest( suite, ShapeSuite_TestSphere2D );
   pcu_suite_addTest( suite, ShapeSuite_TestSphere2D_Invert );
   pcu_suite_addTest( suite, ShapeSuite_TestSphere3D );
   pcu_suite_addTest( suite, ShapeSuite_TestConvexHull2D );
   pcu_suite_addTest( suite, ShapeSuite_TestConvexHull3D );
   pcu_suite_addTest( suite, ShapeSuite_TestUnion2D );
   pcu_suite_addTest( suite, ShapeSuite_TestUnion3D );
   pcu_suite_addTest( suite, ShapeSuite_TestUnion2DSingleNOT );
   pcu_suite_addTest( suite, ShapeSuite_TestUnion3DSingleNOT );
   pcu_suite_addTest( suite, ShapeSuite_TestIntersection2D );
   pcu_suite_addTest( suite, ShapeSuite_TestIntersection3D );
   pcu_suite_addTest( suite, ShapeSuite_TestCylinder );
   pcu_suite_addTest( suite, ShapeSuite_TestPolygonShape2D );
   pcu_suite_addTest( suite, ShapeSuite_TestSuperellipsoid2D );
   pcu_suite_addTest( suite, ShapeSuite_TestSuperellipsoid3D );
}
