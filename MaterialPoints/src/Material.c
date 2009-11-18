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
** $Id: Material.c 532 2008-02-12 01:17:48Z DavidMay $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include <PICellerator/PopulationControl/PopulationControl.h>
#include <PICellerator/Weights/Weights.h>

#include "MaterialPoints.h"

#include <assert.h>
#include <string.h>

const Type Material_Type = "Material";

const Index UNDEFINED_MATERIAL = (unsigned)-1;

/* Public Constructor */
Material* Material_New( 
	Name						name,
	PICelleratorContext*	context,
	Stg_Shape*				shape,
	Dictionary*				materialDictionary,
	Materials_Register*	materialRegister )
{
	Material* self = _Material_DefaultNew( name );

	self->isConstructed = True;
	_Material_Init( self, context, shape, materialDictionary, materialRegister );

	return self;
}


void* _Material_DefaultNew( Name name ) {
	return (void*) _Material_New(
			sizeof(Material),
			Material_Type,
			_Material_Delete, 
			_Material_Print, 
			_Material_Copy, 
			_Material_DefaultNew, 
			_Material_AssignFromXML, 
			_Material_Build, 
			_Material_Initialise, 
			_Material_Execute, 
			_Material_Destroy,
			name,
			NULL,
			NULL,
			NULL );
}


/* Private Constructor */
Material* _Material_New(
	SizeT                                           _sizeOfSelf,
	Type                                            type,
	Stg_Class_DeleteFunction*                       _delete,
	Stg_Class_PrintFunction*                        _print,
	Stg_Class_CopyFunction*                         _copy,
	Stg_Component_DefaultConstructorFunction*       _defaultConstructor,
	Stg_Component_ConstructFunction*                _construct,
	Stg_Component_BuildFunction*                    _build,
	Stg_Component_InitialiseFunction*               _initialise,
	Stg_Component_ExecuteFunction*                  _execute,
	Stg_Component_DestroyFunction*                  _destroy,
	Name                                            name,
	Stg_Shape*                                      shape,
	Dictionary*                                     materialDictionary,
	Materials_Register*                             materialRegister )
{
	Material* self;
	
	assert( _sizeOfSelf >= sizeof(Material) );
	self = (Material*) _Stg_Component_New( 
			_sizeOfSelf,
			type,
			_delete,
			_print,
			_copy,
			_defaultConstructor,
			_construct,
			_build,
			_initialise,
			_execute,
			_destroy,		
			name,
			NON_GLOBAL );

	return self;
}


void _Material_Init( 
	void*						material,
	PICelleratorContext*	context,
	Stg_Shape*				shape,
	Dictionary*				materialDictionary,
	Materials_Register*	materialRegister )
{
	Material* self = (Material*)material;
	
	/* Set Values */
	self->context = context;
	self->dictionary = materialDictionary;
	self->shape = shape;

	/* Register material */
	if (materialRegister != NULL)	
		self->index = Materials_Register_Add( materialRegister, self );	
	else 
		self->index = 0;

	self->extensionMgr = ExtensionManager_New_OfExistingObject( self->name, self );
}

void _Material_Delete( void* material ) {
	Material* self = (Material*) material;

	_Stg_Component_Delete( material );
}


void _Material_Print( void* material, Stream* stream ) {
	Material* self = (Material*) material;

	_Stg_Component_Print( self, stream );

	Journal_PrintPointer( stream, self->dictionary );
	Journal_PrintUnsignedInt( stream, self->index );
}


void* _Material_Copy( void* material, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	Material* self = (Material*) material;
	Material* newMaterial;

	newMaterial = _Stg_Component_Copy( self, dest, deep, nameExt, ptrMap );

	newMaterial->dictionary = self->dictionary;
	newMaterial->index      = self->index;

	return newMaterial;
}


void _Material_AssignFromXML( void* material, Stg_ComponentFactory* cf, void* data ) {
	Material*				self = (Material*) material;
	Dictionary*				materialDictionary;
	Stg_Shape*				shape;
	Materials_Register*	materials_Register;
	PICelleratorContext*	context;

	context = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Context", PICelleratorContext, False, data );
	if( !context ) 
		context = Stg_ComponentFactory_ConstructByName( cf, "context", PICelleratorContext, True, data );

	materialDictionary = Dictionary_GetDictionary( cf->componentDict, self->name );
	shape =  Stg_ComponentFactory_ConstructByKey(  cf,  self->name,  "Shape", Stg_Shape,  True, data  ) ;
	materials_Register = context->materials_Register;

	_Material_Init( self, context, shape, materialDictionary, materials_Register );
}


void _Material_Build( void* material, void* data ) {}
void _Material_Initialise( void* material, void* data ) {}
void _Material_Execute( void* material, void* data ) {}

void _Material_Destroy( void* material, void* data ) {
	Material* self = (Material*) material;

	Stg_Class_Delete( self->extensionMgr );
}


void Material_Layout( void* material, MaterialPointsSwarm* swarm ) {	
	Material*             self               = (Material*) material;
	Stg_Shape*            shape              = self->shape;
	MaterialPoint*        particle;
	Particle_Index        lParticle_I;
	Particle_Index        particleLocalCount = swarm->particleLocalCount;
	Stream*               stream             = Journal_MyStream( Info_Type, self );

	Journal_RPrintf( stream, "Laying out material '%s' within %s '%s':\n", self->name, shape->type, shape->name );
	
	for ( lParticle_I = 0 ; lParticle_I < particleLocalCount ; lParticle_I++ ) {
		particle = (MaterialPoint*)Swarm_ParticleAt( swarm, lParticle_I );
		
		if ( Stg_Shape_IsCoordInside( shape, particle->coord ) ) {
			particle->materialIndex = self->index;
		}
	}
}

double Material_Volume( void* material, IntegrationPointsSwarm* swarm, Coord centroid ) {
	Material*            self               = (Material*)material;
	FeMesh*  	     feMesh             = swarm->mesh;
	ElementType*         elementType;
	IntegrationPoint*    particle;
	Coord                globalCoord;
	Cell_Index           lCell_I;
	Cell_Index           cellLocalCount     = swarm->cellLocalCount;
	Particle_InCellIndex cParticle_I;
	Particle_Index       lParticle_I;
	Material_Index       material_I         = self->index;
	Dimension_Index      dim                = swarm->dim;
	Coord                localCentroid;
	double               detJac;
	double               volume; 
	double               volumeGlobal;

	/* Initialise Sumation of Integral */
	volume = 0.0;	
	memset( localCentroid, 0, sizeof(Coord) );

	/* Loop over all cells in domain */
	for ( lCell_I = 0 ; lCell_I < cellLocalCount ; lCell_I++ ) {
		elementType = FeMesh_GetElementType( feMesh, lCell_I );
		for( cParticle_I = 0 ; cParticle_I < swarm->cellParticleCountTbl[lCell_I] ; cParticle_I++ ) {
			lParticle_I = swarm->cellParticleTbl[lCell_I][cParticle_I];
			particle = (IntegrationPoint*)Swarm_ParticleAt( swarm, lParticle_I );

			/* If particle isn't of the same material type as the material passed in - 
			 * then there is no contribution to the volume */
			if ( material_I != IntegrationPointsSwarm_GetMaterialIndexOn( swarm, particle ) )
				continue;

			/* Calculate Determinant of Jacobian */
			detJac = ElementType_JacobianDeterminant( elementType, feMesh, lCell_I, particle->xi, dim );

			/* Sum Volume */
			volume += detJac * particle->weight;

			FeMesh_CoordLocalToGlobal( feMesh, lCell_I, particle->xi, globalCoord );
				
			/* Sum centroid */
			localCentroid[ I_AXIS ] += detJac * particle->weight * globalCoord[ I_AXIS ];
			localCentroid[ J_AXIS ] += detJac * particle->weight * globalCoord[ J_AXIS ];
			if ( dim == 3 )
				localCentroid[ K_AXIS ] += detJac * particle->weight * globalCoord[ K_AXIS ];
		}
	}

	/* Gather and sum volumes from other processors */
	MPI_Allreduce( &volume, &volumeGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &localCentroid, centroid, dim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

	/* Divide by volume to obtain true centroid */
	centroid[ I_AXIS ] /= volumeGlobal;
	centroid[ J_AXIS ] /= volumeGlobal;
	if ( dim == 3 )
		centroid[ K_AXIS ] /= volumeGlobal;
	else
		centroid[ K_AXIS ] = 0.0;

	return volumeGlobal;
}


void Material_IntegrateField( 
		void*                   material, 
		IntegrationPointsSwarm* swarm, 
		FeVariable*             field, 
		double*                 volumeGlobal, 
		double*                 result ) 
{
	Material*            self               = (Material*)material;
	FeMesh*  	     feMesh             = swarm->mesh;
	ElementType*         elementType;
	IntegrationPoint*    particle;
	Cell_Index           lCell_I;
	Cell_Index           cellLocalCount     = swarm->cellLocalCount;
	Particle_InCellIndex cParticle_I;
	Particle_Index       lParticle_I;
	Material_Index       material_I         = self->index;
	double               detJac;
	double               volume; 
	double*              fieldValue;
	double*              localResult;
	Dof_Index            fieldComponentCount = field->fieldComponentCount;
	Dof_Index            component_I;
	
	/* Allocate memory */
	fieldValue  = Memory_Alloc_Array( double, fieldComponentCount, "fieldValue" );
	localResult = Memory_Alloc_Array( double, fieldComponentCount, "localResult" );

	/* Initialise Sumation of Integral */
	volume = 0.0;	
	memset( localResult, 0, sizeof(double) * fieldComponentCount );

	/* Loop over all cells in domain */
	for ( lCell_I = 0 ; lCell_I < cellLocalCount ; lCell_I++ ) {
		elementType = FeMesh_GetElementType( feMesh, lCell_I );
		for( cParticle_I = 0 ; cParticle_I < swarm->cellParticleCountTbl[lCell_I] ; cParticle_I++ ) {
			lParticle_I = swarm->cellParticleTbl[lCell_I][cParticle_I];
			particle = (IntegrationPoint*) Swarm_ParticleAt( swarm, lParticle_I );

			/* If particle isn't of the same material type as the material passed in - 
			 * then there is no contribution to the volume */
			if ( material_I != IntegrationPointsSwarm_GetMaterialIndexOn( swarm, particle ) )
				continue;

			/* Calculate Determinant of Jacobian */
			detJac = ElementType_JacobianDeterminant( elementType, feMesh, lCell_I, particle->xi, swarm->dim );

			/* Sum Volume */
			volume += detJac * particle->weight;

			/* Integrate field */ 
			FeVariable_InterpolateWithinElement( field, lCell_I, particle->xi, fieldValue );
			for ( component_I = 0 ; component_I < fieldComponentCount ; component_I++ ) {
				localResult[ component_I ] += detJac * particle->weight * fieldValue[ component_I ];
			}
		}
	}

	/* Gather and sum volumes from other processors */
	MPI_Allreduce( localResult, result, fieldComponentCount, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	if ( volumeGlobal )
		MPI_Allreduce( &volume, volumeGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

	Memory_Free( fieldValue );
	Memory_Free( localResult );
}
