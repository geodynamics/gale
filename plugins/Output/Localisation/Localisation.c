/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2005, Monash Cluster Computing 
** All rights reserved.
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
**              * Redistributions of source code must retain the above copyright notice, 
**                      this list of conditions and the following disclaimer.
**              * Redistributions in binary form must reproduce the above copyright 
**                      notice, this list of conditions and the following disclaimer in the 
**                      documentation and/or other materials provided with the distribution.
**              * Neither the name of the Monash University nor the names of its contributors 
**                      may be used to endorse or promote products derived from this software 
**                      without specific prior written permission.
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
*%              Louis Moresi - Louis.Moresi@sci.monash.edu.au
*%
** Contributors:
*+              Robert Turnbull
*+              Vincent Lemiale
*+              Louis Moresi
*+              David May
*+              David Stegman
*+              Mirko Velic
*+              Patrick Sunter
*+              Julian Giordani
*+
** $Id:  $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <StgFEM/FrequentOutput/FrequentOutput.h>

#include "Localisation.h"

unsigned     globPartInc;
MPI_Datatype Localisation_MPI_Datatype;

const Type Underworld_Localisation_Type = "Underworld_Localisation";

void _Underworld_Localisation_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
        UnderworldContext* context;

        context = Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data ); 

        /* Add functions to entry points */
        Underworld_Localisation_Setup( context );
        ContextEP_Append( context, AbstractContext_EP_FrequentOutput, Underworld_Localisation_Output );
}

void _Underworld_Localisation_Build( void* component, void* data ) {
        Underworld_Localisation*        self = (Underworld_Localisation*)component;

        Stg_Component_Build( self->reducedStrainRateFieldInvariantRoot, data, False );

}

void* _Underworld_Localisation_DefaultNew( Name name ) {
        return _Codelet_New(
                sizeof(Underworld_Localisation),
                Underworld_Localisation_Type,
                _Codelet_Delete,
                _Codelet_Print,
                _Codelet_Copy,
                _Underworld_Localisation_DefaultNew,
                _Underworld_Localisation_Construct,
                _Underworld_Localisation_Build,
                _Codelet_Initialise,
                _Codelet_Execute,
                _Codelet_Destroy,
                name );
}

Index Underworld_Localisation_Register( PluginsManager* pluginsManager ) {
        return PluginsManager_Submit( pluginsManager, Underworld_Localisation_Type, "0", _Underworld_Localisation_DefaultNew );
}

void Underworld_Localisation_Setup( UnderworldContext* context ) {
        FieldVariable_Register*              fV_Register               = context->fieldVariable_Register;
        FieldVariable*                       strainRateField;
        Func_Ptr                             _carryOut;
        Dof_Index                            resultDofs;
        Dof_Index                            operandDofs;
        Index                                numberOfOperands;
        Operator*                            ownOperator;
        Dimension_Index                      dim;
        
        Underworld_Localisation* self;

        /* create datatype for MPI communications */
        Localisation_Create_MPI_Datatype();
        
        self = (Underworld_Localisation*)LiveComponentRegister_Get(
                                        context->CF->LCRegister,
                                        Underworld_Localisation_Type );
        
        StgFEM_FrequentOutput_PrintString( context, "Localisation" );

	/* get localisation parameter */
	self->deformationFactor = Dictionary_GetDouble_WithDefault( context->dictionary, "LocalisationDeformationFactor", 0.8 );

        Journal_Firewall( 
                        context->gaussSwarm != NULL, 
                        Underworld_Error,
                        "Cannot find gauss swarm. Cannot use %s.\n", CURR_MODULE_NAME );

        /* create Some FeVariables to determine Localisation */
        strainRateField  = FieldVariable_Register_GetByName( fV_Register, "StrainRateField" );

        dim = strainRateField->dim;
        /* setup required parameters to create new operate */
        resultDofs       = 1 ;
        numberOfOperands = 1 ;
        operandDofs      = ( dim == 2 ? 3 : 6 );
        _carryOut        = ( dim == 2 ? Underworld_Localisation_SymmetricTensor_LowerDimension_InvariantRoot_2d : Underworld_Localisation_SymmetricTensor_LowerDimension_InvariantRoot_3d );
        
        ownOperator = Operator_New( "Localisation_SymmetricTensor_LowerDimension_InvariantRoot", _carryOut, numberOfOperands, operandDofs, resultDofs, dim );

        self->reducedStrainRateFieldInvariantRoot = OperatorFeVariable_NewUnary_OwnOperator(
                        "ReducedStrainRateFieldInvariantRoot",
                        strainRateField, 
                        ownOperator );

        /* add the variables to register so we can checkpoint & examine if necessary */
        FieldVariable_Register_Add( fV_Register, self->reducedStrainRateFieldInvariantRoot );
        
}

void Underworld_Localisation_Output( UnderworldContext* context ) {
        Underworld_Localisation* self;
        double   totalIntegral;
        double   integralSoFar;
        double   weightSoFar;
        double   weightSoFar2;
        double   *min, *max;
        double   translate;
        int meshSizeJ = Dictionary_GetInt_WithDefault( context->dictionary, "elementResJ", 1 );
        int ii;
        int aa;
        int numPoints;
        int numPointsGlob;
        int numPointsProc;
	int myRank;
	int nProcs;
	MPI_Status         status;        
	IntegrandWeightStruct* integrandWeightGlob;
        const int FINISHED_WRITING_TAG = 100;

        MPI_Comm_size( context->communicator, (int*)&nProcs );
        MPI_Comm_rank( context->communicator, (int*)&myRank );

        globPartInc = 0;
        
        self = (Underworld_Localisation*)LiveComponentRegister_Get(
                                context->CF->LCRegister,
                                Underworld_Localisation_Type );

        /* find the size of the domain */
        min = Memory_Alloc_Array_Unnamed( double, Mesh_GetDimSize( self->reducedStrainRateFieldInvariantRoot->feMesh ) );
        max = Memory_Alloc_Array_Unnamed( double, Mesh_GetDimSize( self->reducedStrainRateFieldInvariantRoot->feMesh ) );

        Mesh_GetGlobalCoordRange( self->reducedStrainRateFieldInvariantRoot->feMesh, min, max );
        
        translate = ( max[ J_AXIS ] - min[ J_AXIS ] ) / (2*meshSizeJ);

        /* create arrays to store the integral sum integrands and weights */
        /* first determine how large array needs to be */
        /* factor of 2 or 4 allows for number of particles used to integrate surface */
        if(self->reducedStrainRateFieldInvariantRoot->dim == 2)
           numPoints = 2*Localisation_FindTotElements( self->reducedStrainRateFieldInvariantRoot , J_AXIS, max[ J_AXIS ] - translate );
        else
           numPoints = 2*Localisation_FindTotElements( self->reducedStrainRateFieldInvariantRoot , J_AXIS, max[ J_AXIS ] - translate );
        
        self->integrandWeight  = Memory_Alloc_Array_Unnamed( IntegrandWeightStruct, numPoints );
                
        /* use modified integrate plane routine to do the work!
           these routines store the integrand values, and corresponding weights, for each particle used to 
           construct the integral.. these are then qsort(ed) by integrand value, and summed starting at the top, until 
           a set percentage of the total integral is reached.. the ratio of the weight of all the thus far counted particles, 
           to all the particles, is then the localisation value. */
        totalIntegral = Localisation_IntegratePlane( self->reducedStrainRateFieldInvariantRoot, J_AXIS, max[ J_AXIS ] - translate, self );
	
	/* find the total global number of integration points */
	MPI_Allreduce( &globPartInc, &numPointsGlob, 1, MPI_INT, MPI_SUM, context->communicator );
        
        /* collect data */
        if ( myRank == 0 ) {
                int count = 0;
                /* allocate global list array */
                integrandWeightGlob = Memory_Alloc_Array_Unnamed( IntegrandWeightStruct, numPointsGlob );
                
                /* first copy the rank=0 data into the global array */
                for(ii = 0; ii < globPartInc; ii++){
                        integrandWeightGlob[ii].integrand = self->integrandWeight[ii].integrand;
                        integrandWeightGlob[ii].weight    = self->integrandWeight[ii].weight;
                        count += 1;
                }
                
                /* collect data from other processors. */
                for(aa = 1; aa < nProcs; aa++){
                        /* first get the size of the processors array */
                        MPI_Recv( &numPointsProc, 1, MPI_INT, aa, FINISHED_WRITING_TAG, context->communicator, &status );
                        /* get the data */
                        MPI_Recv( integrandWeightGlob+count, numPointsProc, Localisation_MPI_Datatype, aa, FINISHED_WRITING_TAG, context->communicator, &status );
                        /*find position for next dataset to be placed*/
                        count += numPointsProc;
                }
        } else {
                        /* send the size of this processors array */
                        MPI_Ssend( &globPartInc, 1, MPI_INT, 0, FINISHED_WRITING_TAG, context->communicator );
                        /* send the data */
                        MPI_Ssend( self->integrandWeight, globPartInc, Localisation_MPI_Datatype, 0, FINISHED_WRITING_TAG, context->communicator );
        }
        /* now determine localisation */
        if ( myRank == 0 ) {
                /* sort data array in ascending order (smallest integral contributor placed first in array) */ 
                qsort( integrandWeightGlob, numPointsGlob, sizeof( IntegrandWeightStruct ), _QsortIntegrand );
                /* sum integral contributors, starting from the largest values, and stopping once integral contributions exceed deformationFactor*totalIntegral */
                integralSoFar = 0.;
                weightSoFar   = 0.;
                for(ii = numPointsGlob - 1; ii >= 0 ; ii--){
                        integralSoFar += integrandWeightGlob[ii].integrand;
                        weightSoFar   += integrandWeightGlob[ii].weight;
                        if (integralSoFar > self->deformationFactor*totalIntegral) break;
                }
                /* do sum again this time to find total area */
                weightSoFar2   = 0.;
                for(ii = numPointsGlob - 1; ii >= 0 ; ii--){
                        weightSoFar2   += integrandWeightGlob[ii].weight;
                }
                /* set Localisation to be ratio of area for which 80% of the deformation occurs, to the total area */
                self->Localisation = weightSoFar/weightSoFar2;   
                StgFEM_FrequentOutput_PrintValue( context, self->Localisation );
                
                FreeArray( integrandWeightGlob );
        }
        
        FreeArray( min );
        FreeArray( max );
        FreeArray( self->integrandWeight );
}

void Underworld_Localisation_SymmetricTensor_LowerDimension_InvariantRoot_2d( void* operator, double* operand0, double* result ) {
        Operator* self = (Operator*) operator;
        double    temp;
        
        Operator_FirewallUnary( self );
        Operator_FirewallResultDofs( self, self->operandDofs );
        
        temp      = operand0[0];
        result[0] = fabs( temp );
}

void Underworld_Localisation_SymmetricTensor_LowerDimension_InvariantRoot_3d( void* operator, double* operand0, double* result ) {
        Operator* self = (Operator*) operator;
        SymmetricTensor temp;
        
        Operator_FirewallUnary( self );
        Operator_FirewallResultDofs( self, self->operandDofs );
        
        temp[0] = operand0[0];  /* strainrate_{xx} */
        temp[1] = operand0[2];  /* strainrate_{zz} */
        temp[2] = operand0[4];  /* strainrate_{xz} */
        
        *result = SymmetricTensor_2ndInvariant( temp, 2 );
}


double Localisation_IntegratePlane( void* feVariable, Axis planeAxis, double planeHeight, void* _Localisation  ){
        FeVariable*                fevariable = (FeVariable*) feVariable;
        Underworld_Localisation*      Localisation  = (Underworld_Localisation*) _Localisation;
        IJK                        planeIJK;
        Element_LocalIndex         lElement_I;
        Element_GlobalIndex        gElement_I;
        Element_LocalIndex         elementLocalCount  = FeMesh_GetElementLocalSize( fevariable->feMesh );
        Axis                       aAxis              = ( planeAxis == I_AXIS ? J_AXIS : I_AXIS );
        Axis                       bAxis              = ( planeAxis == K_AXIS ? J_AXIS : K_AXIS );
        double                     integral;
        /* Swarm Stuff */
        Swarm*                     tmpSwarm;
        Bool                       dimExists[]        = { False, False, False };
        ExtensionManager_Register* extensionMgr_Register;
        SingleCellLayout*          singleCellLayout;
        GaussParticleLayout*       gaussParticleLayout;
        Particle_Index             lParticle_I;
        IntegrationPoint*          particle;
        /* Plane location stuff */
        double                     storedXi_J_AXIS;
        Coord                      planeCoord;
        double                     planeXi           = -1;
        double                     planeXiGlobal;
        Index                      planeLayer        = 0;
        Index                      planeLayerGlobal;
        Particle_InCellIndex       particlesPerDim[] = {2,2,2};

        /* Find Elements which plane cuts through */
        memcpy( planeCoord, Mesh_GetVertex( fevariable->feMesh, 0 ), sizeof( Coord ) );
        planeCoord[ planeAxis ] = planeHeight;

        if( Mesh_Algorithms_SearchElements( fevariable->feMesh->algorithms, planeCoord, &lElement_I ) && 
            lElement_I < elementLocalCount )
        {
                Coord           planeXiCoord;

                gElement_I = FeMesh_ElementDomainToGlobal( fevariable->feMesh, lElement_I );
                RegularMeshUtils_Element_1DTo3D( fevariable->feMesh, gElement_I, planeIJK );
                planeLayer = planeIJK[ planeAxis ];
                
                /* Find Local Coordinate of plane */
                FeMesh_CoordGlobalToLocal( fevariable->feMesh, lElement_I, planeCoord, planeXiCoord );
                planeXi = planeXiCoord[ planeAxis ];
        }
        
        /* Should be broadcast */
        MPI_Allreduce( &planeXi,    &planeXiGlobal, 1, MPI_DOUBLE, MPI_MAX, fevariable->communicator );
        MPI_Allreduce( &planeLayer, &planeLayerGlobal, 1, MPI_UNSIGNED, MPI_MAX, fevariable->communicator );

        /* Create Swarm in plane */
        extensionMgr_Register = ExtensionManager_Register_New();
        dimExists[ aAxis ] = True;
        if (fevariable->dim == 3)
                dimExists[ bAxis ] = True;
        
        singleCellLayout = SingleCellLayout_New( "cellLayout", dimExists, NULL, NULL );
        particlesPerDim[ planeAxis ] = 1;
        gaussParticleLayout = GaussParticleLayout_New( "particleLayout", fevariable->dim - 1, particlesPerDim );
        tmpSwarm = Swarm_New( 
                        "tmpgaussSwarm",
                        singleCellLayout, 
                        gaussParticleLayout,
                        fevariable->dim,
                        sizeof(IntegrationPoint), 
                        extensionMgr_Register, 
                        NULL,
                        fevariable->communicator,
                        NULL    );
        Stg_Component_Build( tmpSwarm, NULL, False );

        /* Change Positions of the particles */
        Stg_Component_Initialise( tmpSwarm, NULL, False );
        for ( lParticle_I = 0 ; lParticle_I < tmpSwarm->particleLocalCount ; lParticle_I++ ) {
                particle = (IntegrationPoint*) Swarm_ParticleAt( tmpSwarm, lParticle_I );

                storedXi_J_AXIS = particle->xi[ J_AXIS ];
                particle->xi[ aAxis ]     = particle->xi[ I_AXIS ];
                particle->xi[ bAxis ]     = storedXi_J_AXIS;
                particle->xi[ planeAxis ] = planeXiGlobal;
        }
        
        integral = Localisation_IntegrateLayer_AxisIndependent( fevariable, tmpSwarm, planeAxis, planeLayerGlobal, 
                        fevariable->dim - 1, aAxis, bAxis, planeAxis, Localisation ); 

        /* Delete */
        Stg_Class_Delete( tmpSwarm );
        Stg_Class_Delete( gaussParticleLayout );
        Stg_Class_Delete( singleCellLayout );
        Stg_Class_Delete( extensionMgr_Register );
        
        return integral;
}

double Localisation_IntegrateLayer_AxisIndependent( 
                void* feVariable, void* _swarm,
                Axis layerAxis, Index layerIndex, Dimension_Index dim, 
                Axis axis0, Axis axis1, Axis axis2, void* _Localisation )
{ 
        FeVariable*                fevariable = (FeVariable*)         feVariable;
        Swarm*                     swarm      = (Swarm*)              _swarm;
        Underworld_Localisation*      Localisation  = (Underworld_Localisation*) _Localisation;
        Element_LocalIndex         lElement_I;
        Element_GlobalIndex        gElement_I;
        IJK                        elementIJK;
        double                     elementIntegral;
        double                     integral;
        double                     integralGlobal;
        
        Journal_DPrintf( fevariable->debug, "In %s() for FeVariable \"%s\":\n", __func__, fevariable->name );

        /* Initialise Sumation of Integral */
        integral = 0.0;

        Stream_Indent( fevariable->debug );
        
        for ( gElement_I = 0 ; gElement_I < FeMesh_GetElementGlobalSize( fevariable->feMesh ); gElement_I++ ) {
                RegularMeshUtils_Element_1DTo3D( fevariable->feMesh, gElement_I, elementIJK );

                /* Check if element is in layer plane */
                if ( elementIJK[ layerAxis ] != layerIndex )
                        continue;

                /* Check if element is local */
                if( !FeMesh_ElementGlobalToDomain( fevariable->feMesh, gElement_I, &lElement_I ) || 
                    lElement_I >= FeMesh_GetElementLocalSize( fevariable->feMesh ) )
                {
                        continue;
                }

                elementIntegral = Localisation_IntegrateElement_AxisIndependent( fevariable, swarm, lElement_I, dim, axis0, axis1, axis2,  Localisation  );
                Journal_DPrintfL( fevariable->debug, 2, "Integral of element %d was %f\n", lElement_I, elementIntegral );
                integral += elementIntegral;
        }
        Stream_UnIndent( fevariable->debug );


        /* Gather and sum integrals from other processors */
        MPI_Allreduce( &integral, &integralGlobal, 1, MPI_DOUBLE, MPI_SUM, fevariable->communicator );

        Journal_DPrintf( fevariable->debug, "Calculated global integral of layer %d in Axis %d was %f\n", layerIndex, layerAxis, integralGlobal );
        return integralGlobal;
}

double Localisation_IntegrateElement_AxisIndependent( 
                void* feVariable, void* _swarm, 
                Element_DomainIndex dElement_I, Dimension_Index dim, 
                Axis axis0, Axis axis1, Axis axis2, void* _Localisation ) 
{
        FeVariable*           fevariable = (FeVariable*)           feVariable;
        Swarm*                swarm      = (Swarm*)                _swarm;
        Underworld_Localisation* Localisation  = (Underworld_Localisation*) _Localisation;
        FeMesh*                 feMesh             = fevariable->feMesh;
        FeMesh*                 mesh;
        ElementType*         elementType;
        Cell_LocalIndex      cell_I;
        Particle_InCellIndex cParticle_I;
        Particle_InCellIndex cellParticleCount;
        IntegrationPoint*    particle;
        double               detJac;
        double               integral;
        double               value;
        
        /* Initialise Summation of Integral */
        integral = 0.0;

        /* Use feVariable's mesh as geometry mesh if one isn't passed in */
        if( Stg_Class_IsInstance( feMesh->algorithms, Mesh_CentroidAlgorithms_Type ) )
                mesh = (FeMesh*)((Mesh_CentroidAlgorithms*)feMesh->algorithms)->elMesh;
        else
                mesh = feMesh;
        elementType = FeMesh_GetElementType( mesh, dElement_I );

        /* Determine number of particles in element */
        cell_I = CellLayout_MapElementIdToCellId( swarm->cellLayout, dElement_I );
        cellParticleCount = swarm->cellParticleCountTbl[ cell_I ];

        /* Loop over all particles in element to get integral*/
        for( cParticle_I = 0 ; cParticle_I < cellParticleCount ; cParticle_I++ ) {
                Particle_Index dParticle_I;
                /* Get Pointer to particle */
                particle = (IntegrationPoint*) Swarm_ParticleInCellAt( swarm, cell_I, cParticle_I );

                dParticle_I = swarm->cellParticleTbl[(cell_I)][(cParticle_I)];
                /* Interpolate Value of Field at Particle */
                FeVariable_InterpolateWithinElement( feVariable, dElement_I, particle->xi, &value );

                Journal_DPrintfL( fevariable->debug, 3, "%s: Integrating element %d - particle %d - Value = %g\n", fevariable->name, dElement_I, cParticle_I, value );

                /* Calculate Determinant of Jacobian */
                detJac = ElementType_JacobianDeterminant_AxisIndependent( 
                                elementType, mesh, dElement_I, particle->xi, dim, axis0, axis1, axis2 );
                
                /* store values for later sorting and usage */
                Localisation->integrandWeight[globPartInc].integrand = detJac * particle->weight * value;
                Localisation->integrandWeight[globPartInc].weight    = detJac * particle->weight;
                
                globPartInc += 1 ;
                
                /* Sum Integral */
                integral += detJac * particle->weight * value;
        }
        
        return integral;
}

unsigned Localisation_FindTotElements( void* feVariable, Axis layerAxis, double planeHeight )
{ 
        FeVariable*                fevariable = (FeVariable*)         feVariable;
        Element_LocalIndex         lElement_I;
        Element_GlobalIndex        gElement_I;
        unsigned                   totElements;
        IJK                        elementIJK;
        Coord                      planeCoord;
        Element_LocalIndex         elementLocalCount  = FeMesh_GetElementLocalSize( fevariable->feMesh );
        IJK                        planeIJK;
        double                     planeXi           = -1;
        double                     planeXiGlobal;
        Index                      planeLayer        = 0;
        Index                      planeLayerGlobal;
        
        /* Find Elements which plane cuts through */
        memcpy( planeCoord, Mesh_GetVertex( fevariable->feMesh, 0 ), sizeof( Coord ) );
        planeCoord[ layerAxis ] = planeHeight;

        if( Mesh_Algorithms_SearchElements( fevariable->feMesh->algorithms, planeCoord, &lElement_I ) && 
            lElement_I < elementLocalCount )
        {
                Coord           planeXiCoord;

                gElement_I = FeMesh_ElementDomainToGlobal( fevariable->feMesh, lElement_I );
                RegularMeshUtils_Element_1DTo3D( fevariable->feMesh, gElement_I, planeIJK );
                planeLayer = planeIJK[ layerAxis ];
                
                /* Find Local Coordinate of plane */
                FeMesh_CoordGlobalToLocal( fevariable->feMesh, lElement_I, planeCoord, planeXiCoord );
                planeXi = planeXiCoord[ layerAxis ];
        }
        
        /* Should be broadcast */
        MPI_Allreduce( &planeXi,    &planeXiGlobal, 1, MPI_DOUBLE, MPI_MAX, fevariable->communicator );
        MPI_Allreduce( &planeLayer, &planeLayerGlobal, 1, MPI_UNSIGNED, MPI_MAX, fevariable->communicator );
        
        totElements = 0.;
        
        for ( gElement_I = 0 ; gElement_I < FeMesh_GetElementGlobalSize( fevariable->feMesh ); gElement_I++ ) {
                RegularMeshUtils_Element_1DTo3D( fevariable->feMesh, gElement_I, elementIJK );

                /* Check if element is in layer plane */
                if ( elementIJK[ layerAxis ] != planeLayerGlobal )
                        continue;

                /* Check if element is local */
                if( !FeMesh_ElementGlobalToDomain( fevariable->feMesh, gElement_I, &lElement_I ) || 
                    lElement_I >= FeMesh_GetElementLocalSize( fevariable->feMesh ) )
                {
                        continue;
                }
                
                totElements += 1;

        }
        return totElements;
}


int _QsortIntegrand( const void* _a, const void* _b ) {
	IntegrandWeightStruct* a = (IntegrandWeightStruct*) _a;
	IntegrandWeightStruct* b = (IntegrandWeightStruct*) _b;
	if ( a->integrand > b->integrand )
		return 1;
	else
		return -1;
}
void Localisation_Create_MPI_Datatype() {
	MPI_Datatype          typeList[2] = { MPI_DOUBLE, MPI_DOUBLE };
	int                   blocklen[2] = { 1, 1 };
	MPI_Aint              displacement[2];
	IntegrandWeightStruct integrandWeightStruct;

	displacement[0] = GetOffsetOfMember( integrandWeightStruct, integrand );
	displacement[1] = GetOffsetOfMember( integrandWeightStruct, weight );
	
	MPI_Type_struct( 2, blocklen, displacement, typeList, &Localisation_MPI_Datatype );
	MPI_Type_commit( & Localisation_MPI_Datatype );
}
