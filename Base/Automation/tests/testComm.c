/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 
110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
** $Id: IMap.c 3952 2007-01-09 06:24:06Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include "Base/Foundation/Foundation.h"
#include "Base/IO/IO.h"
#include "Base/Container/Container.h"
#include "Base/Automation/Automation.h"

#include "Base/Foundation/TestBegin.h"


void testSetup( int* argc, char** argv[] ) {
   BaseFoundation_Init( argc, argv );
   BaseIO_Init( argc, argv );
   BaseContainer_Init( argc, argv );
   BaseAutomation_Init( argc, argv );
}

void testTeardown() {
   BaseAutomation_Finalise();
   BaseContainer_Finalise();
   BaseIO_Finalise();
   BaseFoundation_Finalise();
}

TestBegin( Construct ) {
   Comm* comm;

   TestNoAssert( comm = Comm_New() );
   TestTrue( comm );
   TestTrue( comm->mpiComm == MPI_COMM_WORLD );
   TestTrue( comm->recvs == NULL );
   TestTrue( comm->sends == NULL );
   TestTrue( comm->stats == NULL );

  done:
   NewClass_Delete( comm );
}
TestEnd

TestBegin( SetMPIComm ) {
   Comm* comm;

   comm = Comm_New();
   TestNoAssert( Comm_SetMPIComm( comm, MPI_COMM_WORLD ) );
   TestTrue( comm->mpiComm == MPI_COMM_WORLD );

  done:
   NewClass_Delete( comm );
}
TestEnd

TestBegin( SetNbrs ) {
   Comm* comm;
   int nRanks, rank;
   int nNbrs, nbrs[2];

   MPI_Comm_size( MPI_COMM_WORLD, &nRanks );
   MPI_Comm_rank( MPI_COMM_WORLD, &rank );
   if( nRanks > 1 ) {
      if( rank == 0 ) {
	 nNbrs = 1;
	 nbrs[0] = rank + 1;
      }
      else if( rank == nRanks - 1 ) {
	 nNbrs = 1;
	 nbrs[0] = rank - 1;
      }
      else {
	 nNbrs = 2;
	 nbrs[0] = rank - 1;
	 nbrs[1] = rank + 1;
      }
   }
   else
      nNbrs = 0;

   comm = Comm_New();
   Comm_SetMPIComm( comm, MPI_COMM_WORLD );
   TestNoAssert( Comm_SetNeighbours( comm, nNbrs, nbrs ) );

   if( nRanks > 1 ) {
      if( rank == 0 ) {
	 TestTrue( comm->nbrs.size == 1 );
	 TestTrue( comm->nbrs.ptr[0] == rank + 1 );
	 TestTrue( 1 );
      }
      else if( rank == nRanks - 1 ) {
	 TestTrue( comm->nbrs.size == 1 );
	 TestTrue( comm->nbrs.ptr[0] == rank - 1 );
	 TestTrue( 1 );
      }
      else {
	 TestTrue( comm->nbrs.size == 2 );
	 TestTrue( comm->nbrs.ptr[0] == rank - 1 );
	 TestTrue( comm->nbrs.ptr[1] == rank + 1 );
      }
      TestTrue( comm->recvs != NULL );
      TestTrue( comm->sends != NULL );
      TestTrue( comm->stats != NULL );
   }
   else {
      TestTrue( comm->nbrs.size == 0 );
      TestTrue( comm->nbrs.ptr == NULL );
      TestTrue( comm->recvs == NULL );
      TestTrue( comm->sends == NULL );
      TestTrue( comm->stats == NULL );
   }

  done:
   NewClass_Delete( comm );
}
TestEnd

TestBegin( AddNbrs ) {
   Comm* comm;
   int nRanks, rank;
   int nNbrs, nbrs[2];

   MPI_Comm_size( MPI_COMM_WORLD, &nRanks );
   MPI_Comm_rank( MPI_COMM_WORLD, &rank );
   if( nRanks > 1 ) {
      if( rank == 0 ) {
	 nNbrs = 1;
	 nbrs[0] = rank + 1;
      }
      else if( rank == nRanks - 1 ) {
	 nNbrs = 1;
	 nbrs[0] = rank - 1;
      }
      else {
	 nNbrs = 1;
	 nbrs[0] = rank - 1;
      }
   }
   else
      nNbrs = 0;

   comm = Comm_New();
   Comm_SetMPIComm( comm, MPI_COMM_WORLD );
   Comm_SetNeighbours( comm, nNbrs, nbrs );
   if( rank > 0 && rank < nRanks - 1 ) {
      nNbrs = 1;
      nbrs[0] = rank + 1;
   }
   else
      nNbrs = 0;
   TestNoAssert( Comm_AddNeighbours( comm, nNbrs, nbrs ) );

   if( nRanks > 1 ) {
      if( rank == 0 ) {
	 TestTrue( comm->nbrs.size == 1 );
	 TestTrue( comm->nbrs.ptr[0] == rank + 1 );
	 TestTrue( 1 );
      }
      else if( rank == nRanks - 1 ) {
	 TestTrue( comm->nbrs.size == 1 );
	 TestTrue( comm->nbrs.ptr[0] == rank - 1 );
	 TestTrue( 1 );
      }
      else {
	 TestTrue( comm->nbrs.size == 2 );
	 TestTrue( comm->nbrs.ptr[0] == rank - 1 );
	 TestTrue( comm->nbrs.ptr[1] == rank + 1 );
      }
      TestTrue( comm->recvs != NULL );
      TestTrue( comm->sends != NULL );
      TestTrue( comm->stats != NULL );
   }
   else {
      TestTrue( comm->nbrs.size == 0 );
      TestTrue( comm->nbrs.ptr == NULL );
      TestTrue( comm->recvs == NULL );
      TestTrue( comm->sends == NULL );
      TestTrue( comm->stats == NULL );
   }

  done:
   NewClass_Delete( comm );
}
TestEnd

TestBegin( RemNbrs ) {
   Comm* comm;
   int nRanks, rank;
   int nNbrs, nbrs[2];
   IMap mapObj, *map = &mapObj;

   MPI_Comm_size( MPI_COMM_WORLD, &nRanks );
   MPI_Comm_rank( MPI_COMM_WORLD, &rank );
   if( nRanks > 1 ) {
      if( rank == 0 ) {
	 nNbrs = 1;
	 nbrs[0] = rank + 1;
      }
      else if( rank == nRanks - 1 ) {
	 nNbrs = 1;
	 nbrs[0] = rank - 1;
      }
      else {
	 nNbrs = 2;
	 nbrs[0] = rank - 1;
	 nbrs[1] = rank + 1;
      }
   }
   else
      nNbrs = 0;

   comm = Comm_New();
   Comm_SetMPIComm( comm, MPI_COMM_WORLD );
   Comm_SetNeighbours( comm, nNbrs, nbrs );
   if( rank > 0 && rank < nRanks - 1 ) {
      nNbrs = 1;
      nbrs[0] = 0; /* local index */
   }
   else
      nNbrs = 0;
   IMap_Init( map );
   TestNoAssert( Comm_RemoveNeighbours( comm, nNbrs, nbrs, map ) );

   if( nRanks > 1 ) {
      if( rank == 0 ) {
	 TestTrue( comm->nbrs.size == 1 );
	 TestTrue( comm->nbrs.ptr[0] == rank + 1 );
      }
      else if( rank == nRanks - 1 ) {
	 TestTrue( comm->nbrs.size == 1 );
	 TestTrue( comm->nbrs.ptr[0] == rank - 1 );
      }
      else {
	 TestTrue( comm->nbrs.size == 1 );
	 TestTrue( comm->nbrs.ptr[0] == rank + 1 );
      }
      TestTrue( comm->recvs != NULL );
      TestTrue( comm->sends != NULL );
      TestTrue( comm->stats != NULL );
   }
   else {
      TestTrue( comm->nbrs.size == 0 );
      TestTrue( comm->nbrs.ptr == NULL );
      TestTrue( comm->recvs == NULL );
      TestTrue( comm->sends == NULL );
      TestTrue( comm->stats == NULL );
   }

  done:
   NewClass_Delete( comm );
   IMap_Destruct( map );
}
TestEnd

TestBegin( Allgather ) {
   Comm* comm;
   int nRanks, rank;
   int nNbrs, nbrs[2];

   MPI_Comm_size( MPI_COMM_WORLD, &nRanks );
   MPI_Comm_rank( MPI_COMM_WORLD, &rank );
   if( nRanks > 1 ) {
      if( rank == 0 ) {
	 nNbrs = 1;
	 nbrs[0] = rank + 1;
      }
      else if( rank == nRanks - 1 ) {
	 nNbrs = 1;
	 nbrs[0] = rank - 1;
      }
      else {
	 nNbrs = 2;
	 nbrs[0] = rank - 1;
	 nbrs[1] = rank + 1;
      }
   }
   else
      nNbrs = 0;

   comm = Comm_New();
   Comm_SetMPIComm( comm, MPI_COMM_WORLD );
   Comm_SetNeighbours( comm, nNbrs, nbrs );
   /* TODO */

  done:
   NewClass_Delete( comm );
}
TestEnd


#define nTests 5
TestSuite_Test tests[nTests] = {{"construct", testConstruct}, 
				{"set MPI communicator", testSetMPIComm}, 
				{"set neighbours", testSetNbrs}, 
				{"add neighbours", testAddNbrs}, 
				{"remove neighbours", testRemNbrs}};


#include "Base/Foundation/TestEnd.h"
