#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include<string.h>
#include<assert.h>
#include<stdlib.h>


const Type Underworld_OutputForSPModel_Type = "Underworld_OutputForSPModel";

FILE *iFile, *oFile;
unsigned outputTimestep=0;
Index Underworld_OutputForSPModel_Register( PluginsManager* pluginsManager );
void _Underworld_OutputForSPModel_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data );
void* _Underworld_OutputForSPModel_DefaultNew( Name name );
void Underworld_OutputForSPModelDo( UnderworldContext* context );


void Underworld_DoHeightInterpolation( FeMesh* mesh, ElementType* elType, unsigned *nodes, double *Nx, unsigned nodesPerFace, double *height ) {
	double weight, nodeHeight;
	int node_I;
	/* initialise height */
	*height = 0;
	for( node_I = 0 ; node_I < nodesPerFace ; node_I++ ) {
		/* interpolate each surface node's height to the 
		 * wantedCoord using the shapeFunctions */
		weight = Nx[node_I];
		nodeHeight = Mesh_GetVertex(mesh, nodes[node_I])[1];
		*height += (weight*nodeHeight);
	}
}

void myConvertGlobalCoordToElLocal( double *tmpMax, double *tmpMin, double *wantCoord, double *lProjCoord ) {
	/* converts the a global coordinate to local coordinate defined between [-1,1] */
	lProjCoord[0] = ( (2*wantCoord[0]-tmpMin[0]-tmpMax[0])/(tmpMax[0]-tmpMin[0]) );
	lProjCoord[1] = ( (2*wantCoord[2]-tmpMin[2]-tmpMax[2])/(tmpMax[2]-tmpMin[2]) );
}

Bool OutputForSPModel_Recursive( FeMesh* mesh, double* wantCoord, unsigned *minIJK, unsigned *maxIJK, unsigned* set, int axis ) {
	unsigned IJK[3];
	double maxX, halfX;

	/* finish condition */
	if ( (maxIJK[axis] - minIJK[axis]) == 1 ) {
		set[0] = minIJK[axis];
		set[1] = maxIJK[axis];
		return True;
	}

	/* set IJK to be the same as minIJK but in the axis of interest */
	memcpy( IJK, minIJK, 3*sizeof(unsigned) );
	/* new index = x0 + (x1-x0)/2 */
	IJK[axis] =  (double)( minIJK[axis] + ((maxIJK[axis] - minIJK[axis])/2) );

	halfX = Mesh_GetVertex( mesh, RegularMeshUtils_Node_3DTo1D( mesh, IJK ) )[axis];
	maxX = Mesh_GetVertex( mesh, RegularMeshUtils_Node_3DTo1D( mesh, maxIJK ) )[axis];

	if( Num_InRange( wantCoord[axis], halfX, maxX ) ) {
		/* halfX < wantCoord <= maxX */
		return OutputForSPModel_Recursive( mesh, wantCoord, IJK, maxIJK, set, axis );
	} else {
		/* minX < wantCoord <= halfX */
		return OutputForSPModel_Recursive( mesh, wantCoord, minIJK, IJK, set, axis );
	}
	/* report error */
	return 0;
}

int Underworld_OutputForSPModel_InterpolateHeightInXZ( FeMesh* mesh, double* wantCoord ) {
	Grid* vertGrid;
	IArray* inc;
	ElementType* elType;
	double *coord, tmpMin[3], tmpMax[3], lProjCoord[2], Nx[4], tmpGMin[3], tmpGMax[3];
	unsigned setI[2], setK[2], nodes[4], node_I, minIJK[3], maxIJK[3], nodesPerFace, gMinVertOnProc, gMaxVertOnProc;
	int ii, jj, kk, ijk[3];

	inc = IArray_New();

	assert( mesh );

	/* first check if location is on this proc */
	Mesh_GetLocalCoordRange( mesh, tmpMin, tmpMax );
	Mesh_GetGlobalCoordRange( mesh, tmpGMin, tmpGMax );
	if( !Num_InRange( wantCoord[0], tmpMin[0], tmpMax[0] ) || !Num_InRange( wantCoord[2], tmpMin[2], tmpMax[2] ) ) {
		printf("Error ... the coord (x,z) = (%g, %g), is not a local coord!!!\n", wantCoord[0], wantCoord[2] );
		assert(0);
	}

	vertGrid = *(Grid**)ExtensionManager_Get( mesh->info, mesh,
  	ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

	nodesPerFace = 4;
	/* set up IJK parametrisation of nodes */
	minIJK[0] = minIJK[1] = minIJK[2] = 0;

	maxIJK[0] = vertGrid->sizes[0];
	maxIJK[1] = vertGrid->sizes[1];
	maxIJK[2] = vertGrid->sizes[2];

	jj = maxIJK[1]-1;
	/* get global indicies of min max nodes
	 * ASSUME they are the left bottom closest corner and the right heighest farthest corner ??? */
	gMinVertOnProc = Mesh_DomainToGlobal( mesh, MT_VERTEX, 0 ); 
	gMaxVertOnProc = Mesh_DomainToGlobal( mesh, MT_VERTEX, Mesh_GetLocalSize(mesh, MT_VERTEX)-1 ); 
	
	/* convert to ijk definitions */
	RegularMeshUtils_Node_1DTo3D( mesh, gMinVertOnProc, minIJK );
	RegularMeshUtils_Node_1DTo3D( mesh, gMaxVertOnProc, maxIJK );

	/* get the 2 nodes indecies in the I-axis */
	OutputForSPModel_Recursive( mesh, wantCoord, minIJK, maxIJK, setI, 0 );

	/* get the 2 nodes indecies in the K-axis */
	OutputForSPModel_Recursive( mesh, wantCoord, minIJK, maxIJK, setK, 2 );

	/* set y to maximum */
	minIJK[1] = jj;
	/* move in x */
	minIJK[0] = setI[0]; minIJK[2] = setK[0];
	nodes[0] = RegularMeshUtils_Node_3DTo1D( mesh, minIJK );
	minIJK[0] = setI[1]; minIJK[2] = setK[0];
	nodes[1] = RegularMeshUtils_Node_3DTo1D( mesh, minIJK );

	/* move in z */
	minIJK[0] = setI[0]; minIJK[2] = setK[1];
	nodes[2] = RegularMeshUtils_Node_3DTo1D( mesh, minIJK );
	minIJK[0] = setI[1]; minIJK[2] = setK[1];
	nodes[3] = RegularMeshUtils_Node_3DTo1D( mesh, minIJK );

	tmpMax[0] = tmpMin[0] = Mesh_GetVertex(mesh, nodes[0])[0];
	tmpMax[2] = tmpMin[2] = Mesh_GetVertex(mesh, nodes[0])[2];

	for( node_I = 1 ; node_I < nodesPerFace ; node_I++ ) { 
		coord = Mesh_GetVertex( mesh, nodes[node_I] );
		/* x-coord test */
		if( coord[0] > tmpMax[0] )
			tmpMax[0] = coord[0];
		else if (coord[0] < tmpMin[0] ) 
			tmpMin[0] = coord[0];
		/* z-coord test */
		if( coord[2] > tmpMax[2] )
			tmpMax[2] = coord[2];
		else if (coord[2] < tmpMin[2] ) 
			tmpMin[2] = coord[2];
	}

	if( Num_InRange( wantCoord[0], tmpMin[0], tmpMax[0] ) )  {
		/* lets get interested */
		if( Num_InRange( wantCoord[2], tmpMin[2], tmpMax[2] ) ) {
			/* wantCoord is in this face, now take the shapeFunctions of this element
			 * to aid interpolation */
			myConvertGlobalCoordToElLocal( tmpMax, tmpMin, wantCoord, lProjCoord );
			/* evaluate the shape functions at that local spot */
			_BilinearElementType_SF_allNodes( elType, lProjCoord, Nx ); 

			/* now we have the nodes, next interpolate heights */
			Underworld_DoHeightInterpolation( mesh, elType, nodes, Nx, nodesPerFace, &wantCoord[1] );
			return 1;
		}
	}
	/* report failure */
	return 0;
}

void Underworld_OutputForSPModelDo( UnderworldContext* context ) {
	FeMesh* mesh = context->velocityField->feMesh;
	double wantCoord[3], minInputCrd[2], maxInputCrd[2], uw_Coord[3], Trans[3][3];
	double originalHeight, bc;
	int coordCount, line_I = 0;
	char buffer[512];
	
	MPI_Comm comm;
	int      myRank;
	int      nProcs;

	/* don't do anything if this isn't the outputTimestep */
	if( context->timeStep != outputTimestep ) return; 

	double startTime;

	oFile = NULL;
	iFile = NULL; 
	
	char *iFileName = Dictionary_GetString_WithDefault( context->dictionary, "OutputForSPModel_InputFile", "spm.input" );
	char *oFileName = Memory_Alloc_Array_Unnamed( char, 
			/*		outputPath/SPModelMeshOutput.time.dat */ 
			strlen(context->outputPath)+1+strlen("spmMeshOutput")+2+5+1+3 );

	comm = Comm_GetMPIComm( Mesh_GetCommTopology( mesh, MT_VERTEX ) );

	MPI_Comm_size( comm, (int*)&nProcs );
	MPI_Comm_rank( comm, (int*)&myRank );

	sprintf( oFileName, "%s/spmMeshOutput.%.5u.dat", context->outputPath, context->timeStep );

	startTime = MPI_Wtime();

	/* open the input file and get the number of test coord */
	if( (iFile=fopen(iFileName, "r" )) == NULL ) {
		printf("ERROR in %s: couldn't open file %s\n", __func__, iFileName );
		exit(0);
	}

	maxInputCrd[0] = maxInputCrd[1] = 0;
	minInputCrd[0] = minInputCrd[1] = 0;

	while( fgets( buffer, 512, iFile ) != NULL ) { 
		line_I++; 
		/* test for the max and min coord in x and z */
		sscanf( buffer, "%lf %lf", &wantCoord[0], &wantCoord[2] );
		if( wantCoord[0] > maxInputCrd[0] ) maxInputCrd[0] = wantCoord[0];
		else if (wantCoord[0] < minInputCrd[0] ) minInputCrd[0] = wantCoord[0];

		if( wantCoord[2] > maxInputCrd[1] ) maxInputCrd[1] = wantCoord[2];
		else if (wantCoord[2] < minInputCrd[1] ) minInputCrd[1] = wantCoord[2];
	}
	/* allocate the right amount of memory for the coords */
	coordCount = line_I;

	/* close the input file */
	fclose( iFile );


	/* setup affine transformation.
	 * |x'|   | A 0 xt | |x|
	 * |z'| = | 0 B yt | |z|
	 * |1 |   | 0 0 1  | |1|
	 *
	 * where:
	 *  x' and z' are the uw_Coord
	 *  x and z are the spmodel coords
	 *  xt and yt and the scaled offsets in each axis
	 */

	/* setup coord transform matrix */
	Trans[0][0] = Trans[0][1] = Trans[0][2] = 0;
	Trans[1][0] = Trans[1][1] = Trans[1][2] = 0;
	Trans[2][0] = Trans[2][1] = Trans[2][2] = 0;

	/* setup scaling transformation */
	Trans[0][0] = (mesh->maxGlobalCrd[0]-mesh->minGlobalCrd[0]) / (maxInputCrd[0] - minInputCrd[0]);
	Trans[1][1] = (mesh->maxGlobalCrd[2]-mesh->minGlobalCrd[2]) / (maxInputCrd[1] - minInputCrd[1]);

	/* setup translation transformation */
	Trans[0][2] = mesh->minGlobalCrd[0]- (Trans[0][0]*minInputCrd[0]);
	Trans[1][2] = mesh->minGlobalCrd[2]- (Trans[1][1]*minInputCrd[1]);

	/* re-open the input file */
	if( (iFile=fopen(iFileName, "r" )) == NULL ) {
			printf("ERROR in %s: couldn't open file %s\n", __func__, iFileName );
			exit(0);
	}
	
	/* open file if proc 0 */
	if( myRank == 0 ) 
		oFile = fopen( oFileName, "w" );
	else 
		oFile = fopen( oFileName, "a" );

	assert( mesh );

	for( line_I = 0 ; line_I < coordCount ; line_I++ ) {
		fscanf( iFile, "%lf %lf %lf %lf", &(wantCoord[0]), &(wantCoord[2]), &originalHeight, &bc);
		/* do top surface interpolation at the (x-z) - where 
		 * x = wantCoord[0]
		 * z = wantCoord[2]
		 * and the interpolated y will be evaluate in wantCoord[1] */

		/* do an transformation here, into the Underworld coords */
		uw_Coord[0] = Trans[0][0]*wantCoord[0] + Trans[0][1]*wantCoord[2] + Trans[0][2]*1;
		uw_Coord[2] = Trans[1][0]*wantCoord[0] + Trans[1][1]*wantCoord[2] + Trans[1][2]*1;

		if (Underworld_OutputForSPModel_InterpolateHeightInXZ( mesh, uw_Coord ) == 0 ) {
			/* report error */
			printf("\n\nError with the %d-th coord = (%g, %g, %g)\n\n", line_I, wantCoord[0], wantCoord[1], wantCoord[2] );
		}
		/* output is in format (x,y,z). Where y is the interpolated coord */
		fprintf( oFile, "%lf %lf %lf %lf\n", wantCoord[0], uw_Coord[1], wantCoord[2], bc );
	}

	printf("*******************************\n");
	printf("Time taken for OutputForSPModel\n");
	printf("            %g\n", MPI_Wtime()-startTime);
	printf("*******************************\n");

	fclose( iFile );
	fclose( oFile );
	Memory_Free( oFileName );
}
		
Index Underworld_OutputForSPModel_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager,
			Underworld_OutputForSPModel_Type,
			"0",
			_Underworld_OutputForSPModel_DefaultNew );
}

void* _Underworld_OutputForSPModel_DefaultNew( Name name ) {
	return Codelet_New( 
			Underworld_OutputForSPModel_Type,
			_Underworld_OutputForSPModel_DefaultNew,
			_Underworld_OutputForSPModel_AssignFromXML,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
}

void _Underworld_OutputForSPModel_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
	UnderworldContext* context = Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data );

	outputTimestep = Dictionary_GetUnsignedInt_WithDefault( context->dictionary, "OutputForSPModel_OutputTimestep", 10 );

	ContextEP_Append( context, AbstractContext_EP_Dump, Underworld_OutputForSPModelDo );
}



