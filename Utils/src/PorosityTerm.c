

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
/*#include <PICellerator/Voronoi/Voronoi.h>
#include <PICellerator/PopulationControl/PopulationControl.h>
#include <PICellerator/Weights/Weights.h>
#include <PICellerator/MaterialPoints/MaterialPoints.h>
need to include the buoyancy force term (from which the porosity term inherits) also
*/

#include <PICellerator/PICellerator.h>

#include "PorosityTerm.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type PorosityTerm_Type = "PorosityTerm";

PorosityTerm* PorosityTerm_New( 
		Name                                                name,
		ForceVector*                                        forceVector,
		Swarm*                                              integrationSwarm,
		AbstractContext*                                    context,
		Materials_Register*                                 materials_Register )
		/*ExtensionManager_Register*                          particles_Register )*/
{
	PorosityTerm* self = (PorosityTerm*) _PorosityTerm_DefaultNew( name );

	PorosityTerm_InitAll( 
			self,
			forceVector,
			integrationSwarm,
			context,
			materials_Register );
			/*particles_Register );*/

	return self;
}

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
PorosityTerm* _PorosityTerm_New( 
		SizeT                                               sizeOfSelf,  
		Type                                                type,
		Stg_Class_DeleteFunction*                           _delete,
		Stg_Class_PrintFunction*                            _print,
		Stg_Class_CopyFunction*                             _copy, 
		Stg_Component_DefaultConstructorFunction*           _defaultConstructor,
		Stg_Component_ConstructFunction*                    _construct,
		Stg_Component_BuildFunction*                        _build,
		Stg_Component_InitialiseFunction*                   _initialise,
		Stg_Component_ExecuteFunction*                      _execute,
		Stg_Component_DestroyFunction*                      _destroy,
		ForceTerm_AssembleElementFunction*                   _assembleElement,		
		Name                                                name )
{
	PorosityTerm* self;
	
	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(PorosityTerm) );
	self = (PorosityTerm*) _ForceTerm_New( 
		sizeOfSelf, 
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
		_assembleElement,
		name );
	
	/* Function pointers for this class that are not on the parent class should be set here */
	
	return self;
}

void _PorosityTerm_Init( 
		PorosityTerm*                              self, 
		AbstractContext*                                    context,
		Materials_Register*                                 materials_Register )
		/*ExtensionManager_Register*                          particles_Register )*/
{
	self->context             = context;
	/*self->materials_Register  = materials_Register;*/
	/*self->particles_Register  = particles_Register;*/
}

void PorosityTerm_InitAll( 
		void*                                               forceTerm,
		ForceVector*                                        forceVector,
		Swarm*                                              integrationSwarm,
		AbstractContext*                                    context,
		Materials_Register*                                 materials_Register )
		/*ExtensionManager_Register*                          particles_Register )*/
{
	PorosityTerm* self = (PorosityTerm*) forceTerm;

	ForceTerm_InitAll( self, forceVector, integrationSwarm, NULL );
	// TODO: Buoyancy
	_PorosityTerm_Init( self, context, materials_Register );
	/*_PorosityTerm_Init( self, context, particles_Register );*/
}

void _PorosityTerm_Delete( void* forceTerm ) {
	PorosityTerm* self = (PorosityTerm*)forceTerm;

	_ForceTerm_Delete( self );
}

void _PorosityTerm_Print( void* forceTerm, Stream* stream ) {
	PorosityTerm* self = (PorosityTerm*)forceTerm;
	
	_ForceTerm_Print( self, stream );

	/* General info */
	/*Journal_PrintPointer( stream, self->materials_Register );*/
	/*Journal_PrintPointer( stream, self->particles_Register );*/
	Journal_PrintPointer( stream, self->context );
}

void* _PorosityTerm_DefaultNew( Name name ) {
	return (void*)_PorosityTerm_New( 
		sizeof(PorosityTerm), 
		PorosityTerm_Type,
		_PorosityTerm_Delete,
		_PorosityTerm_Print,
		NULL,
		_PorosityTerm_DefaultNew,
		_PorosityTerm_Construct,
		_PorosityTerm_Build,
		_PorosityTerm_Initialise,
		_PorosityTerm_Execute,
		_PorosityTerm_Destroy,
		_PorosityTerm_AssembleElement,
		name );
}

void _PorosityTerm_Construct( void* forceTerm, Stg_ComponentFactory* cf, void* data ) {
	PorosityTerm*      	    self             	= (PorosityTerm*)forceTerm;
	AbstractContext*            context;
	Materials_Register*         materials_Register;
	/*ExtensionManager_Register*  particles_Register;*/

	/* Construct Parent */
	_ForceTerm_Construct( self, cf, data );

	context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data ) ;

	materials_Register = Stg_ObjectList_Get( cf->registerRegister, "Materials_Register" );
	/*particles_Register = Stg_ObjectList_Get( cf->registerRegister, "Particles_Register" );*/
	assert( materials_Register );
	/*assert( particles_Register );*/

	_PorosityTerm_Init( self, context, materials_Register );
	/*_PorosityTerm_Init( self, context, particles_Register );*/

	/* can do this better using the specific component dictionary */ 
	self->referencePorosity = Stg_ComponentFactory_GetDouble( cf, self->name, "phi0", 0.0 );
	self->gravity           = Stg_ComponentFactory_GetDouble( cf, self->name, "porosityGrav", 0.0 );
	self->materialSwarm 	= Stg_ComponentFactory_ConstructByKey( cf, self->name, "MaterialSwarm", MaterialPointsSwarm, False, NULL );
	
	/* shohld really be done in build phase.... dave, 10.08.07 */
	MaterialPointsSwarm*      materialPointSwarm	= self->materialSwarm;
	PorosityExt*     	  particleExt;
	SwarmVariable*	          swarmVariable;
	self->particleExtHandle = ExtensionManager_Add( materialPointSwarm->particleExtensionMgr, "Porosity", sizeof(PorosityExt) );
	particleExt = ExtensionManager_Get( materialPointSwarm->particleExtensionMgr, NULL, self->particleExtHandle );
	swarmVariable = Swarm_NewScalarVariable(
		       	materialPointSwarm,
			"PorositySwarmVariable", 
			(ArithPointer)particleExt, 
			Variable_DataType_Double );
}

void _PorosityTerm_Build( void* forceTerm, void* data ) {
	PorosityTerm*             self                  = (PorosityTerm*)forceTerm;
	MaterialPointsSwarm*      materialPointSwarm	= self->materialSwarm;
	StandardParticle          particle;
	PorosityExt*     	  particleExt;
	SwarmVariable*	          swarmVariable; 

	_ForceTerm_Build( self, data );

	Stg_Component_Build( materialPointSwarm, data, False ); /* ...new... */

	/* Register porosity extensionm */
	/*self->particleExtHandle = ExtensionManager_Add( materialPointSwarm->particleExtensionMgr, "Porosity", sizeof(PorosityExt) );*/

	/* Get a dummy extension offset to allow creating the SwarmVariable */
	/*particleExt = ExtensionManager_Get( materialPointSwarm->particleExtensionMgr, &particle, self->particleExtHandle );*/

	/*swarmVariable = Swarm_NewScalarVariable(
		       	materialPointSwarm,
			"PorositySwarmVariable", 
			(ArithPointer) &particleExt->phi - (ArithPointer) &particle,
			Variable_DataType_Double );*/
}


void _PorosityTerm_Initialise( void* forceTerm, void* data ) {
	PorosityTerm*             self          	= (PorosityTerm*)forceTerm;
	Index                     particle_I;
	Particle*		  particle;
	MaterialPointsSwarm*      swarm			= self->materialSwarm;
	/*Materials_Register*       materials_Register 	= self->materials_Register;*/
	PorosityExt*		  porosityExt;

	_ForceTerm_Initialise( self, data );

	Stg_Component_Initialise( swarm, data, False ); /* ...new... */

	/* have to do all of this in the initialise phase, after the particle layout has been initialised */	
	/*
	for ( particle_I = 0 ; particle_I < swarm->particleLocalCount ; particle_I++) {
		particle = (Particle*)Swarm_ParticleAt( swarm, particle_I );
		porosityExt = (PorosityExt*)ExtensionManager_Get( swarm->particleExtensionMgr, particle, self->particleExtHandle );
		porosityExt->phi = 0.0;
	}*/
}

void _PorosityTerm_Execute( void* forceTerm, void* data ) {
	_ForceTerm_Execute( forceTerm, data );
}

void _PorosityTerm_Destroy( void* forceTerm, void* data ) {
	_ForceTerm_Destroy( forceTerm, data );
}


void _PorosityTerm_AssembleElement( void* forceTerm, ForceVector* forceVector, Element_LocalIndex lElement_I, double* elForceVec ) {
	PorosityTerm*               	     self               = (PorosityTerm*) forceTerm;
	IntegrationPoint*                    intParticle;
	MaterialPoint*                       matParticle;
	/*PorosityTerm_MaterialExt*	     materialExt;*/
	Particle_InCellIndex                 cParticle_I;
	Particle_InCellIndex                 cellParticleCount;
	Element_NodeIndex                    elementNodeCount;
	Dof_Index                        nodeDofCount;
	Dimension_Index                      dim                = forceVector->dim;
	IntegrationPointsSwarm*              intSwarm           = self->integrationSwarm;
	MaterialPointsSwarm*                 matSwarm           = self->materialSwarm;
	FeMesh*       		             mesh               = forceVector->feVariable->feMesh;
	Node_ElementLocalIndex               eNode_I;
	Cell_Index                           cell_I;
	ElementType*                         elementType;
	double                               detJac             = 0.0;
	double                               factor;
	double                               Ni[27];
	double                               porosityBodyForce;
	double			             referenceForce;
	double*                              xi;
	/*Material*                            material;*/
	/*PorosityElement*                     porosityElement;*/
	/*Index                                elementCount;
	Index                                element_I;*/
	PorosityExt*			     porosityExt;
	double				     phi;

	elementType       = FeMesh_GetElementType( mesh, lElement_I );
	elementNodeCount  = elementType->nodeCount;
	nodeDofCount      = dim;
	/*cell_I            = CellLayout_MapElementIdToCellId( swarm->cellLayout, lElement_I );*/
	/*cellParticleCount = swarm->cellParticleCountTbl[cell_I];*/
	cell_I            = CellLayout_MapElementIdToCellId( intSwarm->cellLayout, lElement_I );
	cellParticleCount = intSwarm->cellParticleCountTbl[cell_I];

	referenceForce = self->gravity * self->referencePorosity;

	for( cParticle_I = 0 ; cParticle_I < cellParticleCount ; cParticle_I++ ) {
		/*particle = (IntegrationPoint*) Swarm_ParticleInCellAt( swarm, cell_I, cParticle_I );*/
		intParticle = (IntegrationPoint*) Swarm_ParticleInCellAt( intSwarm, cell_I, cParticle_I );
		matParticle = OneToOneMapper_GetMaterialPoint( intSwarm->mapper, intParticle, &matSwarm );
		/*matParticle = (IntegrationPoint*) Swarm_ParticleInCellAt( matSwarm, cell_I, cParticle_I );*/

		/*porosityExt = (PorosityExt*)ExtensionManager_Get( swarm->particleExtensionMgr, particle, self->particleExtHandle );*/
		porosityExt = (PorosityExt*)ExtensionManager_Get( matSwarm->particleExtensionMgr, matParticle, self->particleExtHandle );
		phi = porosityExt->phi;

		porosityBodyForce = phi * referenceForce;
		/* Get parameters */
		/*material = IntegrationPointsSwarm_GetMaterialOn( (IntegrationPointsSwarm*)swarm, particle );
		materialExt = ExtensionManager_Get( material->extensionMgr, material, self->materialExtHandle );*/
		
		/* Check if this material has porosity term */
		/*elementCount = materialExt->porosityElementCount;
		if ( elementCount == 0 )
			continue;*/

		/* Calculate porosity body force */
		/*porosityBodyForce = 0.0;
        	for ( element_I = 0 ; element_I < elementCount ; element_I++ ) {
			porosityElement = &materialExt->porosityElementList[ element_I ];*/
			/*porosityBodyForce = porosityElement->phi;*/ /* needs to incorporate gravity and reference porosity also */
		/*}*/

		/* Get Values to det integration */
		/*xi  = particle->xi;*/
		xi  = intParticle->xi;
		detJac = ElementType_JacobianDeterminant( elementType, mesh, lElement_I, xi, dim );
		ElementType_EvaluateShapeFunctionsAt( elementType, xi, Ni );

		/* Apply porosity term */
		/*factor = detJac * particle->weight * porosityBodyForce;*/
		
		
		factor = detJac * intParticle->weight * porosityBodyForce;
		for( eNode_I = 0 ; eNode_I < elementNodeCount; eNode_I++ ) {
			elForceVec[ eNode_I * nodeDofCount + J_AXIS ] += factor * Ni[ eNode_I ] ;
		}
	}
}

