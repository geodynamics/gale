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

Index Underworld_OutputForSPModel_Register( PluginsManager* pluginsManager );
void _Underworld_OutputForSPModel_Construct( void* component, Stg_ComponentFactory* cf, void* data );
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

void Underworld_OutputForSPModel_InterpolateHeightInXZ( FeMesh* mesh, double* wantCoord ) {
	Grid* elGrid;
	IArray* inc;
	ElementType* elType;
	double *coord, tmpMin[3], tmpMax[3], lProjCoord[2], Nx[4];
	unsigned elID, nodes[4], node_I, minIJK[3], maxIJK[3], nodesPerFace;
	int ii, jj, kk, ijk[3];

	MPI_Comm comm;
	int      myRank;
	int      nProcs;

	comm = Comm_GetMPIComm( Mesh_GetCommTopology( mesh, MT_VERTEX ) );

	MPI_Comm_size( comm, (int*)&nProcs );
	MPI_Comm_rank( comm, (int*)&myRank );

	inc = IArray_New();

	assert( mesh );
	elGrid = *(Grid**)ExtensionManager_Get( mesh->info, mesh,
  	ExtensionManager_GetHandle( mesh->info, "elementGrid" ) );

	nodesPerFace = 4;
	/* set up IJK parametrisation of nodes */
	minIJK[0] = minIJK[1] = minIJK[2] = 0;

	maxIJK[0] = elGrid->sizes[0]-1;
	maxIJK[1] = elGrid->sizes[1]-1;
	maxIJK[2] = elGrid->sizes[2]-1;

	jj = maxIJK[1];
	for( ii = 0 ; ii < maxIJK[0] ; ii++ ) {

		for( kk = 0 ; kk < maxIJK[1] ; kk++ ) {
			ijk[0] = ii; ijk[1] = jj; ijk[2]= kk; 

			/* get element */
			elID = RegularMeshUtils_Element_3DTo1D( mesh, ijk );
			elType = FeMesh_GetElementType( mesh, elID );
			/* get element's top surface nodes */
			ElementType_GetFaceNodes( elType, (Mesh*)mesh, elID, 1, nodesPerFace, nodes );
			/* test if wantCoord is in that */

			/* initialise tmp min and max coord, in x-z only */
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

			/* ASSUMPTION: if the wantCoord isn't in this x-range, then assume it's not 
			 * in any of these ii elements */
			if( wantCoord[0] < tmpMax[0] && wantCoord[0] > tmpMin[0] )  {
				/* lets get interested */
				if( wantCoord[2] < tmpMax[2] && wantCoord[2] > tmpMin[2] ) {
					/* wantCoord is in this face, now take the shapeFunctions of this element
					 * to aid interpolation */
					myConvertGlobalCoordToElLocal( tmpMax, tmpMin, wantCoord, lProjCoord );
					/* evaluate the shape functions at that local spot */
					_BilinearElementType_SF_allNodes( elType, lProjCoord, Nx ); 

					Underworld_DoHeightInterpolation( mesh, elType, nodes, Nx, nodesPerFace, &wantCoord[1] );
				} else /* this break implements the above assumption */
						break;
			}
		}
	}
}

void Underworld_OutputForSPModelDo( UnderworldContext* context ) {
	FeMesh* mesh = context->velocityField->feMesh;
	double wantCoord[3] = { 0.33333, 0.0, 0.3333 };
	
	MPI_Comm comm;
	int      myRank;
	int      nProcs;

	FILE* oFile = NULL;
	char *oFileName = Memory_Alloc_Array_Unnamed( char, 
			/*		outputPath/SPModelMeshOutput.time.dat */ 
			strlen(context->outputPath)+1+strlen("spmMeshOutput")+1+5+1+3 );

	comm = Comm_GetMPIComm( Mesh_GetCommTopology( mesh, MT_VERTEX ) );

	MPI_Comm_size( comm, (int*)&nProcs );
	MPI_Comm_rank( comm, (int*)&myRank );

	sprintf( oFileName, "%s/spmMeshOutput.%.5u.dat", context->outputPath, context->timeStep );
	/* open file if proc 0 */
	if( myRank == 0 ) 
		oFile = fopen( oFileName, "w" );
	else 
		oFile = fopen( oFileName, "a" );

	assert( mesh );

	/* do top surface interpolation at the (x-z) - where 
	 * x = wantCoord[0]
	 * z = wantCoord[2]
	 * and the interpolated y will be evaluate in wantCoord[1] */
	Underworld_OutputForSPModel_InterpolateHeightInXZ( mesh, wantCoord );

	printf(" the interpolated height at (x,z) - (%.5g, %.5g) = %.5g\n", wantCoord[0], wantCoord[2], wantCoord[1] );
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
			_Underworld_OutputForSPModel_Construct,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
}

void _Underworld_OutputForSPModel_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	UnderworldContext* context = Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data );
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput, Underworld_OutputForSPModelDo );
}

