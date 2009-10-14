/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2005, Monash Cluster Computing 
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
** Contributors:
*+		Robert Turnbull
*+		Vincent Lemiale
*+		Louis Moresi
*+		David May
*+		David Stegman
*+		Mirko Velic
*+		Patrick Sunter
*+		Julian Giordani
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

const Type Underworld_AverageTemperature_Type = "Underworld_AverageTemperature";
void Underworld_AverageTemperature_PrintHeaderToFile( void* context );
void Underworld_AverageTemperature_Output( void* _context );

void _Underworld_AverageTemperature_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	UnderworldContext*  context;

	context = Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data );

	Underworld_AverageTemperature_PrintHeaderToFile( context );
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput, Underworld_AverageTemperature_Output );
}

void* _Underworld_AverageTemperature_DefaultNew( Name name ) {
	return Codelet_New(
		Underworld_AverageTemperature_Type,
		_Underworld_AverageTemperature_DefaultNew,
		_Underworld_AverageTemperature_Construct,
		_Codelet_Build,
		_Codelet_Initialise,
		_Codelet_Execute,
		_Codelet_Destroy,
		name );
}

Index Underworld_AverageTemperature_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Underworld_AverageTemperature_Type, "0", _Underworld_AverageTemperature_DefaultNew );
}

void Underworld_AverageTemperature_Output( void* _context ) {
	UnderworldContext* context       = (UnderworldContext*) _context;
	FeVariable*        temperatureFe = (FeVariable*) LiveComponentRegister_Get( context->CF->LCRegister, "temperatureField" );
	FeMesh*		   mesh         = temperatureFe->feMesh;
	IntegrationPointsSwarm* swarm    = (IntegrationPointsSwarm*)LiveComponentRegister_Get( context->CF->LCRegister, "gaussSwarm" );
	IntegrationPoint*  particle;
	ElementType*       elementType;
	Element_LocalIndex lElement_I, lCell_I;
	Cell_LocalIndex    cParticle_I;
	Index              cellParticleCount, dim;
	double             particleTemperature;
	double             lTempIntegral = 0;
	double             lAvg = 0;
	double             simulationBuoyancyVolumeIntegral = 0;
	double             lBuoyancyVolumeIntegral = 0;
	double             simulationAvg = 0;
	double             detJac;

	dim = context->dim;

	for( lElement_I = 0 ; lElement_I < FeMesh_GetElementLocalSize( mesh ); lElement_I++ ) {
		elementType       = FeMesh_GetElementType( mesh, lElement_I );
		lCell_I           = CellLayout_MapElementIdToCellId( swarm->cellLayout, lElement_I );
		cellParticleCount = swarm->cellParticleCountTbl[ lCell_I ];
		/* get particles in each element that makes up the patch */
		for( cParticle_I = 0 ; cParticle_I < cellParticleCount ; cParticle_I++ ) {
			particle = (IntegrationPoint*) Swarm_ParticleInCellAt( swarm, lCell_I, cParticle_I );
			FeVariable_InterpolateWithinElement( temperatureFe, lElement_I, particle->xi, &particleTemperature );

			if( particleTemperature > 1e-6 ) {
				/* Calculate Determinant of Jacobian */
				detJac = ElementType_JacobianDeterminant( elementType, mesh, lElement_I, particle->xi, dim );
				/* Sum Integral */
				lBuoyancyVolumeIntegral += detJac * particle->weight;
				lTempIntegral += detJac * particle->weight * particleTemperature;
			}
		}
	}

	lAvg = lTempIntegral / lBuoyancyVolumeIntegral; 
	
	MPI_Reduce( &lAvg, &simulationAvg, 1, MPI_DOUBLE, MPI_SUM, 0, context->communicator );
	MPI_Reduce( &lBuoyancyVolumeIntegral, &simulationBuoyancyVolumeIntegral, 1, MPI_DOUBLE, MPI_SUM, 0, context->communicator );
	simulationAvg = simulationAvg / context->nproc;
	simulationBuoyancyVolumeIntegral = simulationBuoyancyVolumeIntegral / context->nproc;
	StgFEM_FrequentOutput_PrintValue( context, simulationAvg );
	StgFEM_FrequentOutput_PrintValue( context, simulationBuoyancyVolumeIntegral );
}

void Underworld_AverageTemperature_PrintHeaderToFile( void* context ) {
	StgFEM_FrequentOutput_PrintString( context, "AvgBuoyancyTemperature" );
	StgFEM_FrequentOutput_PrintString( context, "BuoyancyVolume" );
}

