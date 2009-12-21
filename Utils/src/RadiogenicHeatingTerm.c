

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PopulationControl/PopulationControl.h>
#include <PICellerator/Weights/Weights.h>
#include <PICellerator/MaterialPoints/MaterialPoints.h>
#include <PICellerator/Utils/Utils.h>

#include "types.h"
#include "RadiogenicHeatingTerm.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type RadiogenicHeatingTerm_Type = "RadiogenicHeatingTerm";

RadiogenicHeatingTerm* RadiogenicHeatingTerm_New( 
		Name                                                name,
		ForceVector*                                        forceVector,
		Swarm*                                              integrationSwarm,
		AbstractContext*                                    context,
		Materials_Register*                                 materials_Register )
{
	RadiogenicHeatingTerm* self = (RadiogenicHeatingTerm*) _RadiogenicHeatingTerm_DefaultNew( name );

	RadiogenicHeatingTerm_InitAll( 
			self,
			forceVector,
			integrationSwarm,
			context,
			materials_Register );

	return self;
}

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
RadiogenicHeatingTerm* _RadiogenicHeatingTerm_New( 
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
	RadiogenicHeatingTerm* self;
	
	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(RadiogenicHeatingTerm) );
	self = (RadiogenicHeatingTerm*) _ForceTerm_New( 
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

void _RadiogenicHeatingTerm_Init( 
		RadiogenicHeatingTerm*                              self, 
		AbstractContext*                                    context,
		Materials_Register*                                 materials_Register )
{
	self->context             = context;
	self->materials_Register  = materials_Register;
}

void RadiogenicHeatingTerm_InitAll( 
		void*                                               forceTerm,
		ForceVector*                                        forceVector,
		Swarm*                                              integrationSwarm,
		AbstractContext*                                    context,
		Materials_Register*                                 materials_Register )
{
	RadiogenicHeatingTerm* self = (RadiogenicHeatingTerm*) forceTerm;

	ForceTerm_InitAll( self, forceVector, integrationSwarm, NULL );
	_RadiogenicHeatingTerm_Init( self, context, materials_Register );
}

void _RadiogenicHeatingTerm_Delete( void* forceTerm ) {
	RadiogenicHeatingTerm* self = (RadiogenicHeatingTerm*)forceTerm;

	_ForceTerm_Delete( self );
}

void _RadiogenicHeatingTerm_Print( void* forceTerm, Stream* stream ) {
	RadiogenicHeatingTerm* self = (RadiogenicHeatingTerm*)forceTerm;
	
	_ForceTerm_Print( self, stream );

	/* General info */
	Journal_PrintPointer( stream, self->materials_Register );
	Journal_PrintPointer( stream, self->context );
}

void* _RadiogenicHeatingTerm_DefaultNew( Name name ) {
	return (void*)_RadiogenicHeatingTerm_New( 
		sizeof(RadiogenicHeatingTerm), 
		RadiogenicHeatingTerm_Type,
		_RadiogenicHeatingTerm_Delete,
		_RadiogenicHeatingTerm_Print,
		NULL,
		_RadiogenicHeatingTerm_DefaultNew,
		_RadiogenicHeatingTerm_Construct,
		_RadiogenicHeatingTerm_Build,
		_RadiogenicHeatingTerm_Initialise,
		_RadiogenicHeatingTerm_Execute,
		_RadiogenicHeatingTerm_Destroy,
		_RadiogenicHeatingTerm_AssembleElement,
		name );
}

void _RadiogenicHeatingTerm_Construct( void* forceTerm, Stg_ComponentFactory* cf, void* data ) {
	RadiogenicHeatingTerm*      self             = (RadiogenicHeatingTerm*)forceTerm;
	AbstractContext*            context;
	Materials_Register*         materials_Register;

	/* Construct Parent */
	_ForceTerm_Construct( self, cf, data );

	context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data ) ;

	materials_Register = Stg_ObjectList_Get( cf->registerRegister, "Materials_Register" );
	assert( materials_Register );

	_RadiogenicHeatingTerm_Init( self, context, materials_Register );
}

void _RadiogenicHeatingTerm_Build( void* forceTerm, void* data ) {
	RadiogenicHeatingTerm*               self               = (RadiogenicHeatingTerm*)forceTerm;
	RadiogenicHeatingTerm_MaterialExt*   materialExt;
	Material_Index                       material_I;
	Material*                            material;
	HeatingElement*                      heatingElement;
	Index                                heatingElement_I;
	Index                                heatingElementCount;
	Dictionary_Entry_Value*              list;
	Dictionary_Entry_Value*              entry;
	Materials_Register*                  materials_Register = self->materials_Register;

	_ForceTerm_Build( self, data );

	/* Sort out material extension stuff */
	self->materialExtHandle = Materials_Register_AddMaterialExtension( 
			self->materials_Register, self->type, sizeof(RadiogenicHeatingTerm_MaterialExt) );
	for ( material_I = 0 ; material_I < Materials_Register_GetCount( materials_Register ) ; material_I++) {
		material = Materials_Register_GetByIndex( materials_Register, material_I );
		materialExt = ExtensionManager_GetFunc( material->extensionMgr, material, self->materialExtHandle );

		/* Get List of Heating Elements from material's dictionary */
        	list = Dictionary_Get( material->dictionary, "heatingElements" );
    		heatingElementCount = Dictionary_Entry_Value_GetCount( list );
        	materialExt->heatingElementList = Memory_Alloc_Array( HeatingElement, heatingElementCount, "Heating Element" );
    		memset( materialExt->heatingElementList, 0, heatingElementCount * sizeof(HeatingElement) );
    		materialExt->heatingElementCount = heatingElementCount;
    	
    		for ( heatingElement_I = 0 ; heatingElement_I < heatingElementCount ; heatingElement_I++) { 
			heatingElement = &materialExt->heatingElementList[ heatingElement_I ];
    			entry = Dictionary_Entry_Value_GetElement( list, heatingElement_I );
	    	
	    		heatingElement->Q = Dictionary_Entry_Value_AsDouble( Dictionary_Entry_Value_GetMember( entry, "Q"));
	    		heatingElement->lambda = Dictionary_Entry_Value_AsDouble( Dictionary_Entry_Value_GetMember( entry, "lambda"));
	   }
	}
}

void _RadiogenicHeatingTerm_Initialise( void* forceTerm, void* data ) {
	RadiogenicHeatingTerm*             self             = (RadiogenicHeatingTerm*)forceTerm;

	_ForceTerm_Initialise( self, data );
}

void _RadiogenicHeatingTerm_Execute( void* forceTerm, void* data ) {
	_ForceTerm_Execute( forceTerm, data );
}

void _RadiogenicHeatingTerm_Destroy( void* forceTerm, void* data ) {
	_ForceTerm_Destroy( forceTerm, data );
}


void _RadiogenicHeatingTerm_AssembleElement( void* forceTerm, ForceVector* forceVector, Element_LocalIndex lElement_I, double* elForceVec ) {
	RadiogenicHeatingTerm*               self               = (RadiogenicHeatingTerm*) forceTerm;
	IntegrationPoint*                    particle;
	RadiogenicHeatingTerm_MaterialExt*   materialExt;
	Particle_InCellIndex                 cParticle_I;
	Particle_InCellIndex                 cellParticleCount;
	Element_NodeIndex                    elementNodeCount;
	Dimension_Index                      dim                = forceVector->dim;
	Swarm*                               swarm              = self->integrationSwarm;
	FeMesh*       		             mesh               = forceVector->feVariable->feMesh;
	Node_ElementLocalIndex               eNode_I;
	Cell_Index                           cell_I;
	ElementType*                         elementType;
	double                               detJac             = 0.0;
	double                               factor;
	double                               Ni[27];
	double                               radiogenicHeating;
	double*                              xi;
	Material*                            material;
	double                               time               = self->context->currentTime;
	HeatingElement*                      heatingElement;
	Index                                heatingElementCount;
	Index                                heatingElement_I;

	elementType       = FeMesh_GetElementType( mesh, lElement_I );
	elementNodeCount  = elementType->nodeCount;
	cell_I            = CellLayout_MapElementIdToCellId( swarm->cellLayout, lElement_I );
	cellParticleCount = swarm->cellParticleCountTbl[cell_I];

	for( cParticle_I = 0 ; cParticle_I < cellParticleCount ; cParticle_I++ ) {
		particle = (IntegrationPoint*) Swarm_ParticleInCellAt( swarm, cell_I, cParticle_I );

		/* Get parameters */
		material = IntegrationPointsSwarm_GetMaterialOn( (IntegrationPointsSwarm*)swarm, particle );
		materialExt = ExtensionManager_Get( material->extensionMgr, material, self->materialExtHandle );
		
		/* Check if this material has heating term */
		heatingElementCount = materialExt->heatingElementCount;
		if ( heatingElementCount == 0 )
			continue;

		/* Calculate Radiogenic Heating */
		radiogenicHeating = 0.0;
        	for ( heatingElement_I = 0 ; heatingElement_I < heatingElementCount ; heatingElement_I++ ) {
			heatingElement = &materialExt->heatingElementList[ heatingElement_I ];
	        	radiogenicHeating += heatingElement->Q * exp(-heatingElement->lambda * (time));
		}

		/* Get Values to det integration */
		xi  = particle->xi;
		detJac = ElementType_JacobianDeterminant( elementType, mesh, lElement_I, xi, dim );
		ElementType_EvaluateShapeFunctionsAt( elementType, xi, Ni );

		/* Apply heating term */
		factor = detJac * particle->weight * radiogenicHeating;
		for( eNode_I = 0 ; eNode_I < elementNodeCount; eNode_I++ ) { 		
			elForceVec[ eNode_I ] += factor * Ni[ eNode_I ] ;
		}
	}
}

