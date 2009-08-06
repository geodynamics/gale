#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h>
#include "StgDomain/StgDomain.h"

typedef struct {
	DomainContext* context;
} ShapeSuiteData;

void ShapeSuite_Setup( ShapeSuiteData* data ) { }

void ShapeSuite_Teardown( ShapeSuiteData* data ) {
	//Stg_Component_Destroy( data->context, 0 /* dummy */, False );
	//Stg_Class_Delete( data->context );
	/* remove generated dat file */
	//if( remove("output/test.dat") != 0 ) { pcu_assert(0); }
	//if( remove("output/input.xml") != 0 ) { pcu_assert(0); }
	//remove("output/input.xml");
}


void ShapeSuite_GenerateBox( ShapeSuiteData* data, Dimension_Index dim, char* testFileName ) {
/** Test Definition: */
	DomainContext*   context = NULL;
	Dictionary*      dictionary;
	Stg_Shape*       shape;
	unsigned         testCoordCount, index;
	Name             outputPath;
	Coord            coord;
	Stream*          stream = Journal_Register( Info_Type, testFileName );
	char xml_input[PCU_PATH_MAX];

	/* read in the xml input file */
	pcu_filename_input( "testBox.xml", xml_input );
	context = (DomainContext*)stgMainInitFromXML( xml_input, MPI_COMM_WORLD );
	dictionary = context->dictionary;
	outputPath = Dictionary_GetString( dictionary, "outputPath" );
	Stream_RedirectFile_WithPrependedPath( stream, outputPath, testFileName );
	shape   = (Stg_Shape*) LiveComponentRegister_Get( context->CF->LCRegister, "shape" );
	assert( shape );

	testCoordCount = Dictionary_GetUnsignedInt_WithDefault( dictionary, "testCoordCount", 10000 );
	//dim = Dictionary_GetUnsignedInt( dictionary, "dim" );

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
	ShapeSuite_GenerateBox( data, dim, "test2D.dat" );
	
	pcu_filename_expected( "testBox2D.expected", expected_file );
	pcu_check_fileEq( "output/test2D.dat", expected_file );
	//remove("output/test2D.dat");
}
void ShapeSuite_TestBox3D( ShapeSuiteData* data ) {
	Dimension_Index  dim = 3;
	char expected_file[PCU_PATH_MAX];
	ShapeSuite_GenerateBox( data, dim, "test3D.dat" );

	pcu_filename_expected( "testBox3D.expected", expected_file );
	pcu_check_fileEq( "output/test3D.dat", expected_file );
	//remove("output/test3D.dat");
}

void ShapeSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, ShapeSuiteData );
   pcu_suite_setFixtures( suite, ShapeSuite_Setup, ShapeSuite_Teardown );
   pcu_suite_addTest( suite, ShapeSuite_TestBox2D );
   pcu_suite_addTest( suite, ShapeSuite_TestBox3D );
}
