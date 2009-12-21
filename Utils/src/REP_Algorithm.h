/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2007, Monash Cluster Computing 
** All rights reserved.
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
** 		* Redistributions of source code must retain the above copyright notice, 
** 			this list of conditions and the following disclaimer.
** 		* Redistributions in binary form must reproduce the above copyright 
**			notice, this list of conditions and the following disclaimer in the 
**			documentation and/or other materials provided with the distribution.
** 		* Neither the name of the Monash University nor the names of its contributors 
**			may be used to endorse or promote products derived from this software 
**			without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
** THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
** PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS 
** BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
** CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
** SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
** HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
** LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
** OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
**
** Contact:
*%		Louis Moresi - Louis.Moresi@sci.monash.edu.au
*%
**
** Author:
**		Julian Giordani - Julian.Giordani@sci.monash.edu.au
**
** Contributors:
*+		Robert Turnbull
*+		Vincent Lemiale
*+		Louis Moresi
*+		David May
*+		David Stegman
*+		Mirko Velic
*+		Patrick Sunter
*+
** $Id:  $
*/


#ifndef __Underworld_Utils_REP_Algorithm_h__
#define __Underworld_Utils_REP_Algorithm_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type REP_Algorithm_Type;

	/** REP_Algorithm class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */

	#define __REP_Algorithm \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__Stg_Component \
		/* Virtual functions go here */ \
		/* REP_Algorithm Data Structures */ \
		UnderworldContext*	   context; \
		FeVariable*		   velocityField; \
		FeVariable*		   pressureField; \
		RecoveredFeVariable**      repFieldList; \
	  int                        repFieldCount; \
		IntegrationPointsSwarm*    IPswarm; \
		BoundaryNodesInfo*         boundaryNodesInfo; \
		ConstitutiveMatrix*        constitutiveMatrix; \
		/* Entry Point Setup */ \
		Name                       recoveryEPName; \
		EntryPoint*                recoveryEP; \
		EntryPoint_Register*       entryPoint_Register; \
		/* General crap */ \
		Stream*                    myStream; \
		IArray*		                 incArray;
		
	struct REP_Algorithm { __REP_Algorithm };

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define REP_ALGORITHM_DEFARGS \
                STG_COMPONENT_DEFARGS

	#define REP_ALGORITHM_PASSARGS \
                STG_COMPONENT_PASSARGS

	REP_Algorithm* _REP_Algorithm_New(  REP_ALGORITHM_DEFARGS  );
	
	void _REP_Algorithm_Delete( void* rep );
	void _REP_Algorithm_Print( void* rep, Stream* stream );

	void* _REP_Algorithm_DefaultNew( Name name );
	void _REP_Algorithm_AssignFromXML( void* rep, Stg_ComponentFactory* cf, void* data ) ;
	void _REP_Algorithm_Build( void* rep, void* data ) ;
	void _REP_Algorithm_Initialise( void* rep, void* data ) ;
	void _REP_Algorithm_Execute( void* _context, void* data ) ;
	void _REP_Algorithm_Destroy( void* rep, void* data ) ;
	/* Copy */
	#define REP_Algorithm_Copy( self ) \
		(REP_Algorithm*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define REP_Algorithm_DeepCopy( self ) \
		(REP_Algorithm*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
	
	void* _REP_Algorithm_Copy( void* rep, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

	/* Private Functions */
	void _REP_Algorithm_Init( void* rep,
				  UnderworldContext* context );

	void _REP_Algorithm_FlagBoundaryNodes( REP_Algorithm* rep ) ;
	void _REP_Algorithm_CalculateRecoveredField( void* rep );

	void _REP_Algorithm_NormaliseCoord2D( double* nodePosition, double* minXYZ, double* maxXYZ, double* pVec );
	void _REP_Algorithm_NormaliseCoord3D( double* nodePosition, double* minXYZ, double* maxXYZ, double* pVec );

	void _REP_Algorithm_PopulateElementNeighbourList( FeMesh* mesh, 
						Node_LocalIndex patchNodeID, 
						Index neighbourElCount, 
						Index* el_List );

	void REP_Algorithm_AssembleAllLocalElement( REP_Algorithm* self );
	void REP_Add2LmStruct( LmStruct* patchList, Index nodeID );

	void _REP_Algorithm_Solver(double **array, double* bVec, int n);
	void _REP_Algorithm_SwapOutAssemblyTerms( void* rep);
	void _REP_Algorithm_ZeroHiFi( REP_Algorithm* self, double*** elHi_Mat, double*** elFi_Mat );
	void _REP_Algorithm_AssembleElement( REP_Algorithm* self, int lElement_I, double**** elHi_Mat, double*** elFi_Mat, double** p_pVec, double** GNx);
	void _REP_Algorithm_AssembleAllLocalElement( REP_Algorithm* self );
	void _REP_Algorithm_PutElementObjectsIntoProcObject( REP_Algorithm* self, int lElement_I, double*** Hi_Mat, double*** Fi_Mat );
	void _REP_Algorithm_Communicate( REP_Algorithm* self );
	void _REP_Algorithm_Boundaries( REP_Algorithm* self );
	void _REP_Algorithm_SolvePatches( REP_Algorithm* self );
	int _REP_Algorithm_locateInPatchListStruct( LmStruct* lmStruct, int nodeID );
	void REP_Algorithm_MakeLM( FeMesh* mesh, int lPatchID, int *elList, int* nEl, LmStruct* lmStruct );
	void _REP_Algorithm_CommunicateBoundaries( REP_Algorithm* self );
	void _REP_Algorithm_DoBoundaries( RecoveredFeVariable* self, int lNodeID, BoundaryNodesInfo* bNodeInfo );


	void dgesv_( int* N, int* NRHS, double* AT, int* LDA, int* IPIV, double* bVec, int* LDB, int* info);

#endif

