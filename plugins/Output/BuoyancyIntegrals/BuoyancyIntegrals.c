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
** $Id: BuoyancyIntegrals.c 487 2007-06-07 05:48:32Z LukeHodkinson $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <StgFEM/FrequentOutput/FrequentOutput.h>


/* prototypes */
Index Underworld_BuoyancyIntegrals_Register( PluginsManager *pluginsManager );
void* _Underworld_BuoyancyIntegrals_DefaultNew( Name name );
void _Underworld_BuoyancyIntegrals_CTX_Delete( void *component );
void _Underworld_BuoyancyIntegrals_AssignFromXML( void *component, Stg_ComponentFactory *cf, void *data );
void Underworld_BuoyancyIntegrals_Output( UnderworldContext *context );
void Underworld_BuoyancyIntegrals_Setup( void* _context );

typedef struct {
	__Codelet
			double int_w_bar_dt;
	double beta;
	double gravity;
	int dim;
	double y_b_initial;
	FieldVariable *temperatureField;
	double x_b, z_b;
	MaterialPointsSwarm *cob_swarm; /* center of buouyancy swarm */
} Underworld_BuoyancyIntegrals_CTX;



/*---------------------------------- PLUGIN SOURCE CODE --------------------------------------*/

const Type Underworld_BuoyancyIntegrals_Type = "Underworld_BuoyancyIntegrals";

Index Underworld_BuoyancyIntegrals_Register( PluginsManager *pluginsManager ) 
{
	return PluginsManager_Submit( pluginsManager, Underworld_BuoyancyIntegrals_Type, "0", _Underworld_BuoyancyIntegrals_DefaultNew );
}

void* _Underworld_BuoyancyIntegrals_DefaultNew( Name name ) 
{
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(Underworld_BuoyancyIntegrals_CTX);
	Type                                                      type = Underworld_BuoyancyIntegrals_Type;
	Stg_Class_DeleteFunction*                              _delete = _Underworld_BuoyancyIntegrals_CTX_Delete;
	Stg_Class_PrintFunction*                                _print = _Codelet_Print;
	Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Underworld_BuoyancyIntegrals_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _Underworld_BuoyancyIntegrals_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _Codelet_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _Codelet_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _Codelet_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _Codelet_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;

	
	return _Codelet_New(  CODELET_PASSARGS  );
}

void _Underworld_BuoyancyIntegrals_CTX_Delete( void *component ) 
{
	Underworld_BuoyancyIntegrals_CTX  *ctx = (Underworld_BuoyancyIntegrals_CTX*)component;
	
	_Codelet_Delete( ctx );
}


void _Underworld_BuoyancyIntegrals_AssignFromXML( void *component, Stg_ComponentFactory *cf, void *data ) 
{
	UnderworldContext *context;
	Underworld_BuoyancyIntegrals_CTX *ctx;
	MaterialPointsSwarm *cob_swarm; /* center of buouyancy swarm */
	
	
	
	context = Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data ); 
	
	/* Add functions to entry points */
	ContextEP_Append( context, AbstractContext_EP_AssignFromXMLExtensions, Underworld_BuoyancyIntegrals_Setup );
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput, Underworld_BuoyancyIntegrals_Output );
	
	
	ctx = (Underworld_BuoyancyIntegrals_CTX*)LiveComponentRegister_Get(
			context->CF->LCRegister,
			Underworld_BuoyancyIntegrals_Type );
	
	/* Look for a swarm to which we will assign the calculated center of buoyancy to */
	ctx->cob_swarm = NULL;
	cob_swarm = Stg_ComponentFactory_ConstructByName( cf, "center_buoyancy_swarm", MaterialPointsSwarm, False, data );
	if( cob_swarm != NULL ) {
		ctx->cob_swarm = cob_swarm;
	}
	
}

void Underworld_BuoyancyIntegrals_Setup( void *_context )
{
	UnderworldContext *context = (UnderworldContext*) _context;
	FieldVariable_Register *fV_Register = context->fieldVariable_Register;
	Underworld_BuoyancyIntegrals_CTX *ctx;
	Stg_Shape *shape;
	
	
	/* allocate memory */
	ctx = (Underworld_BuoyancyIntegrals_CTX*)LiveComponentRegister_Get(
			context->CF->LCRegister,
			Underworld_BuoyancyIntegrals_Type );
	
	/* init values */
	ctx->dim = Stg_ComponentFactory_GetRootDictInt( context->CF, "dim", -1 );
	if( (int)ctx->dim == -1 ) {
		printf("******************** ERROR dim IS UNINITIALISED ******************************** \n");
	}
	
	
	ctx->beta = Stg_ComponentFactory_GetRootDictDouble( context->CF, "alpha", -1 );
	if( (int)ctx->beta == -1 ) {
		printf("******************** ERROR ALPHA IS UNINITIALISED ******************************** \n");
	}
	
	ctx->gravity = Stg_ComponentFactory_GetRootDictDouble( context->CF, "gravity", -1 );
	if( (int)ctx->gravity == -1 ) {
		printf("******************** ERROR GRAVITY IS UNINITIALISED ******************************** \n");
	}
	
	ctx->y_b_initial = Stg_ComponentFactory_GetRootDictDouble( context->CF, "y_b_initial", -1 );
	if( (int)ctx->y_b_initial == -1 ) {
		printf("********************* ERROR Y_B_INITIAL IS NOT SET *********************** \n");
	}

	ctx->int_w_bar_dt = ctx->y_b_initial;
	
	if (ctx->dim ==3){
		shape = (Stg_Shape*)Stg_ComponentFactory_ConstructByName( context->CF, "cylinder", Stg_Shape, True, 0 /* dummy */ );
		ctx->x_b = shape->centre[0];
		ctx->z_b = shape->centre[2];
	} else if (ctx->dim==2){
		shape = (Stg_Shape*)Stg_ComponentFactory_ConstructByName( context->CF, "disk", Stg_Shape, True, 0 /* dummy */ );
		ctx->x_b = shape->centre[0];
	}
	
	ctx->temperatureField = FieldVariable_Register_GetByName( fV_Register, "TemperatureField" );
	
	
	
#if 0
	/* Create Some FeVariables to calculate nusselt number */
	temperatureField          = FieldVariable_Register_GetByName( fV_Register, "TemperatureField" );
	velocityField             = FieldVariable_Register_GetByName( fV_Register, "VelocityField" );
	
	
	/* To compute B */
	/* perturbed_temperatre := T - T_0 */
	/*
	self->perturbed_temperatre_field = OperatorFeVariable_NewUnary(
	"AdvectiveHeatFluxField",
	temperatureField,
	"VectorScale" );
	*/
	/* 
	Does not appear to be anything to shift fields by a scalar. 
	Will just assume that T0 = 0
	*/
	self->perturbed_temperatre_field = temperatureField;
	
	
	/* To compute w_bar */
	/* grab the second component of the velocity field ez \dot u */
	self->vertical_velocity_field = (FeVariable*) OperatorFeVariable_NewUnary(
			"VerticalVelocityField",
			velocityField,
			"TakeSecondComponent" );
	self->vertical_velocity_field->feMesh = ((FeVariable*)velocityField)->feMesh;
	
	self->vT_field = OperatorFeVariable_NewBinary(
			"VerticalAdvectiveTemperatureField",
			temperatureField,
			self->vertical_velocity_field,
			"VectorScale" );
	self->vT_field->feMesh = ((FeVariable*)velocityField)->feMesh;
#endif
	
	
	
	/* Add names to the integrals in the frequent output */
	StgFEM_FrequentOutput_PrintString( context, "B" );
	StgFEM_FrequentOutput_PrintString( context, "w_bar" );
	StgFEM_FrequentOutput_PrintString( context, "x_b" );
	StgFEM_FrequentOutput_PrintString( context, "y_b" );
	if (ctx->dim==3){
	StgFEM_FrequentOutput_PrintString( context, "z_b" );
	}
	StgFEM_FrequentOutput_PrintString( context, "int_w_bar_dt" );
	StgFEM_FrequentOutput_PrintString( context, "temp_b" );
	StgFEM_FrequentOutput_PrintString( context, "temp_max" );
	
	
	/* Set the initial coordinates for the cob swarm */
	if( ctx->cob_swarm != NULL ) {
		Swarm *swarm;
		int point_count;
		GlobalParticle *particle;
		
		/* Force swarm to allocate points so that I can set initial values */
		Stg_Component_Build( ctx->cob_swarm, NULL, False );
		Stg_Component_Initialise( ctx->cob_swarm, NULL, False );
		
		swarm = (Swarm*)ctx->cob_swarm;
		/* Set initial coords based on existing data from input xml */
		particle = (GlobalParticle*)Swarm_ParticleAt( swarm, 0 );
		particle->coord[0] = ctx->x_b;
		particle->coord[1] = ctx->y_b_initial;
		particle->coord[2] = ctx->z_b;
	}
	
	
	
	
}


void perform_integrals( UnderworldContext *context, double *B, double *w_bar, double *y_b, double *int_w_bar_dt )
{
	Underworld_BuoyancyIntegrals_CTX *ctx;
	IntegrationPoint *ip;
	double *xi;
	double weight;
	Particle_InCellIndex p, ngp;
	Cell_Index cell_I;
	double velocity[3], global_coord[3];
	FeMesh* mesh;
	Element_LocalIndex e;
	
	double i_T, i_v, i_y, i_vT; /* interpolated quantity */
	double sum_T, sum_vT, sum_yT; /* integral sum */
	double _sum_T, _sum_vT, _sum_yT; /* element sums */
	double g_sum_T, g_sum_vT, g_sum_yT; /* element sums */
	int n_elements;
	ElementType *elementType;
	Node_Index elementNodeCount;
	Dimension_Index dim;
	double det_jac, dt;
	double **GNx;
	double _sum_vol, sum_vol;
	
	FeVariable *velocityField, *temperatureField;
	Swarm* gaussSwarm;
	
	ctx = (Underworld_BuoyancyIntegrals_CTX*)LiveComponentRegister_Get(
			context->CF->LCRegister,
			Underworld_BuoyancyIntegrals_Type );
	
	velocityField = (FeVariable*)LiveComponentRegister_Get( context->CF->LCRegister, "VelocityField" );
	temperatureField = (FeVariable*)LiveComponentRegister_Get( context->CF->LCRegister, "temperatureField" );
	gaussSwarm = (Swarm*)LiveComponentRegister_Get( context->CF->LCRegister, "gaussSwarm" );

	/* initialise values to compute */
	*B = *w_bar = *y_b = -1.0;
	*int_w_bar_dt = ctx->int_w_bar_dt;
	
	
	
#if 0
	/* test 1 - get beta */
	*B = ctx->beta;	
	/* test 2 - get gravity */
	*w_bar = ctx->gravity;
	/* test 3 get dt */
	*y_b = context->dt;
#endif	
	
	/* evaluate some integrals */
	
#if 0
	ngp = context->gaussSwarm->particleLocalCount;
	cell_I = 0;
	for( p = 0; p<ngp; p++ ) {
		ip = (IntegrationPoint*)Swarm_ParticleInCellAt( context->gaussSwarm, cell_I, p );
		xi = ip->xi;
		weight = ip->weight;
		///*
		printf("[%d] w = %f : xi = %f %f \n", p, weight, xi[0], xi[1] );
		//*/
	}
#endif
	
	/* assuming all elements are the same */
	dim = ctx->dim;
	elementType = FeMesh_GetElementType( velocityField->feMesh, 0 );
	elementNodeCount = elementType->nodeCount;
	GNx = Memory_Alloc_2DArray( double, dim, elementNodeCount, "Global Shape Function Derivatives for mayhem" );

	mesh = temperatureField->feMesh;
	
	
	_sum_T = _sum_vT = _sum_yT = 0.0;
	_sum_vol = 0.0;
	
	ngp = gaussSwarm->particleLocalCount;
	n_elements = FeMesh_GetElementLocalSize( mesh );
	//	printf("n_elements = %d \n", n_elements );
	
	for( e=0; e<n_elements; e++ ) {
		cell_I = CellLayout_MapElementIdToCellId( gaussSwarm->cellLayout, e );
		elementType = FeMesh_GetElementType( mesh, e );
		
		sum_T  = 0.0;
		sum_vT = 0.0;
		sum_yT = 0.0;
		sum_vol = 0.0;
		
		i_T = i_v = i_vT = i_y = 0.0;
		
		for( p=0; p<ngp; p++ ) {
			ip = (IntegrationPoint*)Swarm_ParticleInCellAt( gaussSwarm, cell_I, p );
			xi = ip->xi;
			weight = ip->weight;
			
			ElementType_ShapeFunctionsGlobalDerivs(
					elementType,
					mesh, e,
					xi, dim, &det_jac, GNx );
			
			FeVariable_InterpolateFromMeshLocalCoord( temperatureField, mesh, e, xi, &i_T );
			FeVariable_InterpolateFromMeshLocalCoord( velocityField,    mesh, e, xi, velocity );
			
			i_v = velocity[1];
			i_vT = i_v * i_T;
			
			FeMesh_CoordLocalToGlobal( mesh, e, xi, global_coord );
			i_y = global_coord[1];
			
			//printf("%f %f %f %f J=%f \n", i_T, i_v, i_vT, i_y, det_jac );
			
			sum_T  = sum_T  + weight * i_T * det_jac;
			sum_vT = sum_vT + weight * i_vT * det_jac;
			sum_yT = sum_yT + weight * i_y * i_T * det_jac;
			sum_vol = sum_vol + weight * det_jac;
		}
		
		_sum_T = _sum_T + sum_T;
		_sum_vT = _sum_vT + sum_vT;
		_sum_yT = _sum_yT + sum_yT;
		_sum_vol = _sum_vol + sum_vol;		
		
	}
	//printf("%f %f %f \n", _sum_T, _sum_vT, _sum_yT );
	
	
	
	/* all reduce */
	MPI_Allreduce ( &_sum_T, &g_sum_T, 1, MPI_DOUBLE, MPI_SUM, context->communicator );
	MPI_Allreduce ( &_sum_vT, &g_sum_vT, 1, MPI_DOUBLE, MPI_SUM, context->communicator );
	MPI_Allreduce ( &_sum_yT, &g_sum_yT, 1, MPI_DOUBLE, MPI_SUM, context->communicator );
	/*
	MPI_Allreduce ( &_sum_vol, &g_sum_vol, 1, MPI_DOUBLE, MPI_SUM, context->communicator );
	printf("l_sum_vol = %f \n", _sum_vol );
	printf("g_sum_vol = %f \n", g_sum_vol );
	*/
	
	*B     = ctx->gravity * ctx->beta * g_sum_T;
	*w_bar = (ctx->gravity * ctx->beta/ *B) * g_sum_vT;
	*y_b   = (ctx->gravity * ctx->beta/ *B) * g_sum_yT;
	
	dt = context->dt;
	*int_w_bar_dt = *int_w_bar_dt + (*w_bar) * dt;
	ctx->int_w_bar_dt = *int_w_bar_dt;
	
	Memory_Free( GNx );
	
}



void eval_temperature( UnderworldContext *context, double y_b, double *temp_b )
{
	Underworld_BuoyancyIntegrals_CTX *ctx;
	double global_coord[3];
	InterpolationResult result;
	double T;
	FeVariable* temperatureField;

	T = -66.99;
	
	ctx = (Underworld_BuoyancyIntegrals_CTX*)LiveComponentRegister_Get(
			context->CF->LCRegister,
			Underworld_BuoyancyIntegrals_Type );
	
	temperatureField = (FeVariable*)LiveComponentRegister_Get( context->CF->LCRegister, "temperatureField" );
	/* Get x_b, and z_b from xml */
	/* "cylinder" z_b = CentreZ (0.5), x_b = CentreX (1.0) */
	if (ctx->dim==3){
		global_coord[0] = ctx->x_b;
		global_coord[1] = y_b;
		global_coord[2] = ctx->z_b;
	}
    	if (ctx->dim==2){
		global_coord[0] = ctx->x_b;
		global_coord[1] = y_b;
	}

	result = FieldVariable_InterpolateValueAt( temperatureField, global_coord, &T );
	MPI_Allreduce ( &T, temp_b, 1, MPI_DOUBLE, MPI_MAX, context->communicator );

}


void assign_coords_to_swarm( double x_b, double y_b, double z_b, MaterialPointsSwarm *cob_swarm )
{
	Swarm *swarm;
	int point_count;
	GlobalParticle *particle;
	Particle_Index lParticle_I;
	Stream *errorStream = Journal_Register( Error_Type, "Underworld_BuoyancyIntegrals: assign_coords_to_swarm" );
	int rank;
	
	
	/* Cast to get parent */
	swarm = (Swarm*)cob_swarm;
	
	
	/* check cob swarm only has one point */
	point_count = swarm->particleLocalCount;
	/*
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	printf("rank[%d]: np (cob) = %d \n", rank, point_count );
	*/
	Journal_Firewall(
			point_count == 1,
			errorStream,
			"Error in %s:\n"
			"Swarm with name %s is be used to plot center of buoyancy.\n"
			"This swarm contains more then 1 point. This is unexpected!!", __func__, swarm->name );
	
	
	
	
	/* get point to the first particle */
	lParticle_I = 0;
	particle = (GlobalParticle*)Swarm_ParticleAt( swarm, lParticle_I );
	
	/* Set the coorindates */
	/*printf("*********** ASSIGNING COORDS TO COB ****************** \n");*/
	particle->coord[0] = x_b;
	particle->coord[1] = y_b;
	particle->coord[2] = z_b;
	
}


void Underworld_BuoyancyIntegrals_Output( UnderworldContext *context ) 
{
	Underworld_BuoyancyIntegrals_CTX *ctx;
	double B, w_bar, y_b, int_w_bar_dt;
	double temp_b; /* the temperature at (x_b,y_b,z_b) */
	double temp_max;
	FeVariable* temperatureField;
	
	perform_integrals( context, &B, &w_bar, &y_b, &int_w_bar_dt );
	eval_temperature( context, y_b, &temp_b );
	
	ctx = (Underworld_BuoyancyIntegrals_CTX*)LiveComponentRegister_Get(
			context->CF->LCRegister,
			Underworld_BuoyancyIntegrals_Type );
	
	temperatureField = (FeVariable*)LiveComponentRegister_Get( context->CF->LCRegister, "temperatureField" );
	
	StgFEM_FrequentOutput_PrintValue( context, B );
	StgFEM_FrequentOutput_PrintValue( context, w_bar );
	StgFEM_FrequentOutput_PrintValue( context, ctx->x_b );
	StgFEM_FrequentOutput_PrintValue( context, y_b );
	if (ctx->dim==3){
		StgFEM_FrequentOutput_PrintValue( context, ctx->z_b );
	}
	StgFEM_FrequentOutput_PrintValue( context, int_w_bar_dt );
	StgFEM_FrequentOutput_PrintValue( context, temp_b );
	
	temp_max = _FeVariable_GetMaxGlobalFieldMagnitude( temperatureField );
	StgFEM_FrequentOutput_PrintValue( context, temp_max );
	
	
	/* assign coords of center of buoyancy to swarm (if it exists) */
	if( ctx->cob_swarm != NULL ) {
		assign_coords_to_swarm( ctx->x_b, y_b, ctx->z_b, ctx->cob_swarm );
	}
	
	
	
}


