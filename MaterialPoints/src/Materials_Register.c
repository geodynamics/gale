/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
** Copyright (c) 2005-2006, Monash Cluster Computing, Building 28, Monash University Clayton Campus,
**	Victoria, 3800, Australia
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**
** Contributors:
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Patrick D. Sunter, Software Engineer, VPAC. (patrick@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	David Stegman, Postdoctoral Fellow, Monash University. (david.stegman@sci.monash.edu.au)
**	Wendy Sharples, PhD Student, Monash University (wendy.sharples@sci.monash.edu.au)
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
** $Id: Materials_Register.c 532 2008-02-12 01:17:48Z DavidMay $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include <PICellerator/PopulationControl/PopulationControl.h>
#include <PICellerator/Weights/Weights.h>

#include "MaterialPoints.h"

#include <string.h>
#include <assert.h>
#include <math.h>

const Type Materials_Register_Type = "Materials_Register";

Materials_Register* Materials_Register_New( void ) {
	Materials_Register* self;

	self = (Materials_Register*) _NamedObject_Register_New(
		sizeof(Materials_Register),
		Materials_Register_Type,
		_Materials_Register_Delete,
		_Materials_Register_Print,
		_Materials_Register_Copy );

	return self;
}

void _Materials_Register_Delete( void* _materialsRegister ) {
	Materials_Register* self = (Materials_Register*) _materialsRegister;

	_NamedObject_Register_Delete( self );
}
void _Materials_Register_Print( void* _materialsRegister, Stream* stream ) {
	Materials_Register* self = (Materials_Register*) _materialsRegister;

	_NamedObject_Register_Print( self, stream );
}

void* _Materials_Register_Copy( void* _materialsRegister, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	Materials_Register* self                   = (Materials_Register*) _materialsRegister;
	Materials_Register* newMaterials_Register;

	newMaterials_Register = (Materials_Register*)_NamedObject_Register_Copy( self, dest, deep, nameExt, ptrMap );

	return newMaterials_Register;
}


Material* Materials_Register_GetByIndex( Materials_Register* self, Index materialIndex ) {
	if ( UNDEFINED_MATERIAL == materialIndex ) {
		return NULL;
	}

	return (Material*) NamedObject_Register_GetByIndex( self, materialIndex );
}


void Materials_Register_SetupSwarm( void* materialRegister, MaterialPointsSwarm* swarm )
{
	Materials_Register* self = (Materials_Register*)materialRegister;

	assert( swarm != NULL );

	_Materials_Register_LayoutGeometry( self, swarm );
/* 	Materials_Register_SetupParticleToMaterialMappings( self, swarm ); */
	Materials_Register_AssignParticleProperties( self, swarm, swarm->swarmVariable_Register->variable_Register );
}


ExtensionInfo_Index Materials_Register_AddMaterialExtension( void* materialsRegister, Type type, SizeT extensionSize ) {
	Materials_Register*     self         = (Materials_Register*) materialsRegister;
	Material*               material;
	ExtensionInfo_Index     firstResult  = 0;
	ExtensionInfo_Index     result       = 0;
	Material_Index          material_I   = 0;

	for ( material_I = 0 ; material_I < Materials_Register_GetCount( self ) ; material_I++) {
		material = Materials_Register_GetByIndex( self, material_I );

		result = ExtensionManager_Add( material->extensionMgr, type, extensionSize );

		if ( material_I == 0 )
			firstResult = result;
		else {
			Journal_Firewall( 
				firstResult == result, 
				Journal_Register( Error_Type, self->type ),
				"Material '%s' has a different number of extensions than eariler materials.\n", 
				material->name );
		}
	}

	return result;
}
			
void _Materials_Register_LayoutGeometry( void* materialsRegister, void* swarm ) {
	Materials_Register*     self         = (Materials_Register*) materialsRegister;
	Material_Index          material_I;
	Material*               material;
	
	for ( material_I = 0 ; material_I < Materials_Register_GetCount( self ) ; material_I++) {
		material = Materials_Register_GetByIndex( self, material_I );
		Material_Layout( material, swarm );
	}
}

void Materials_Register_AssignParticleProperties( 
		void*                   materialRegister,
		MaterialPointsSwarm*    swarm,
		Variable_Register*      variableRegister )
{
	Materials_Register* self               = (Materials_Register*)materialRegister;
	Material*           material;
	Particle_Index      lParticle_I;
	Particle_Index      particleLocalCount = swarm->particleLocalCount;
	Particle_Index      particleGlobalCount = 0;
	Stream*             stream = Journal_Register( Info_Type, self->type );
	double              setupStartTime = 0;
	double              setupTime = 0, tmin, tmax;
	Processor_Index     formerStreamPrintingRank;
	unsigned int        numberOfCompletionPrintIncrements=10;
	double              completionRatioIncrement= 1 / (double)numberOfCompletionPrintIncrements;
	double              nextCompletionRatioToPrint=0;
	Particle_Index      nextCompletedParticleCountToPrint=0;
	Particle_Index      nextPlusOneCompletedParticleCountToPrint=0;
	char                *title;

        formerStreamPrintingRank = Stream_GetPrintingRank( stream );
        Stream_SetPrintingRank( stream, 0 );

	Journal_Printf( stream, "In func %s(): for swarm \"%s\"\n", __func__, swarm->name );
	Stream_Indent( stream );
	setupStartTime = MPI_Wtime();
	MPI_Reduce( &particleLocalCount, &particleGlobalCount, 1, MPI_UNSIGNED, MPI_SUM, 0, swarm->comm );
	Journal_Printf( stream, "Assigning initial particle properties to the %u global particles\n",
		particleGlobalCount );
	Stream_Indent( stream );

	nextCompletionRatioToPrint = completionRatioIncrement;
	nextCompletedParticleCountToPrint = ceil(particleLocalCount * nextCompletionRatioToPrint - 0.001 );

	for ( lParticle_I = 0 ; lParticle_I < particleLocalCount ; lParticle_I++ ) {
		material = MaterialPointsSwarm_GetMaterialAt( swarm, lParticle_I );

		Journal_Firewall( 
				material != NULL, 
				Journal_Register( Error_Type, self->type ),
				"In func %s: Cannot find material for particle '%u'\n", 
				__func__, 
				lParticle_I );
	
		/* Loop through material's dictionary assigning values to the variables of this particle */
		Variable_Register_SetAllVariablesFromDictionary( variableRegister, lParticle_I, material->dictionary );

		
		if ( /*(swarm->myRank == 0) && */ ((lParticle_I+1) >= nextCompletedParticleCountToPrint ) ) {
/* 			 TODO: parallelise : non-master CPUs send a non-blocking update to the master to report */
/* 			 status. Master does blocking receive on all updates before printing */

			/* Special case for really small swarms, or really small increments - may cross more than one
				at once */
			nextPlusOneCompletedParticleCountToPrint = ceil(( particleLocalCount
				* (nextCompletionRatioToPrint + completionRatioIncrement )) - 0.001 );

			while ( (lParticle_I+1) >= nextPlusOneCompletedParticleCountToPrint )
			{
				nextCompletionRatioToPrint += completionRatioIncrement;
				nextPlusOneCompletedParticleCountToPrint = ceil(( particleLocalCount
					* (nextCompletionRatioToPrint + completionRatioIncrement )) - 0.001 );
				if ( nextCompletionRatioToPrint >= 1.0 ) {
					nextCompletionRatioToPrint = 1.0;
					break;
				}
			}
			Journal_Printf( stream, "done %.0f%% (%u particles)...\n", 
				(nextCompletionRatioToPrint * 100),
				lParticle_I+1 );
			nextCompletionRatioToPrint += completionRatioIncrement;
			nextCompletedParticleCountToPrint = ceil(particleLocalCount * nextCompletionRatioToPrint - 0.001);
		}
	}

	Stream_UnIndent( stream );
	/* Need this barrier so the time is accurate */
	MPI_Barrier( swarm->comm );
	setupTime = MPI_Wtime() - setupStartTime;

	Stream_UnIndent( stream );

	MPI_Reduce( &setupTime, &tmin, 1, MPI_DOUBLE, MPI_MIN, 0, swarm->comm );
	MPI_Reduce( &setupTime, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, swarm->comm );
	Journal_Printf( stream, "%s(): finished setup of material properties for swarm \"%s\"\n"
		"\ttook %g [min] / %g [max] secs\n", __func__, swarm->name, tmin, tmax );
	Stream_SetPrintingRank( stream, formerStreamPrintingRank );
}


void Variable_SetValueFromDictionary( void* _variable, Index index, Dictionary* dictionary ) {
	Variable*          variable = (Variable*)_variable;

	Variable_Update( variable );
	
	/* Assign variable from dictionary according to data type */
	switch (variable->dataTypes[0]) {
		case Variable_DataType_Char:
			Variable_SetValueChar(  variable, index, Dictionary_GetUnsignedInt( dictionary, variable->name ));
			break;
		case Variable_DataType_Short:
			Variable_SetValueShort( variable, index, Dictionary_GetUnsignedInt( dictionary, variable->name ));
			break;
		case Variable_DataType_Int:
			Variable_SetValueShort( variable, index, Dictionary_GetInt( dictionary, variable->name ));
			break;
		case Variable_DataType_Float:
			Variable_SetValueFloat(  variable, index, Dictionary_GetDouble( dictionary, variable->name ));
			break;
		case Variable_DataType_Double:
			Variable_SetValueDouble( variable, index, Dictionary_GetDouble( dictionary, variable->name ));
			break;
		default: {
			Journal_Printf( 
				Journal_MyStream( Error_Type, variable ), 
				"In func %s: Unable to set value of %s from dictionary.", 
				__func__, 
				variable->name );
		}
	}
}

void Variable_Register_SetAllVariablesFromDictionary( void* _variable_Register, Index index, Dictionary* dictionary ) {
	Variable_Register* variable_Register = (Variable_Register*) _variable_Register;
	Variable*          variable;
	Dictionary_Index   dictionary_I;
	Dictionary_Entry*  entry;

	for ( dictionary_I = 0 ; dictionary_I < dictionary->count ; dictionary_I++ ) {
		entry = dictionary->entryPtr[ dictionary_I ];

		/* Find variable on particle associated with this name */
		variable = Variable_Register_GetByName( variable_Register, entry->key );
		
		/* If there is no variable, continue through loop */
		if (variable == NULL) 
			continue;

		/* Assign Value */
		Variable_SetValueFromDictionary( variable, index, dictionary );
	}
}


