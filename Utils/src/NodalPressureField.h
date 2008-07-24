/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Monash Cluster Computing, Australia
** (C) 2003-2004 All Rights Reserved
**
** Primary Authors:
** Luke Hodkinson
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Underworld_Utils_NodalPressureField_h__
#define __Underworld_Utils_NodalPressureField_h__

/*
** Class definition.
*/

extern const Type NodalPressureField_Type;

#define __NodalPressureField			\
   __ParticleFeVariable				\
						\
   /* Virtual functions. */			\
						\
   /* Members. */				\
   Variable_Register* variable_Register;	\
   FeVariable* pressureField;

struct NodalPressureField { __NodalPressureField };

#define NODALPRESSUREFIELD_ARGS			\
   PARTICLEFEVARIABLE_ARGS

#define NODALPRESSUREFIELD_PASSARGS		\
   PARTICLEFEVARIABLE_PASSARGS

/*
** Constructors/Destructors.
*/

NodalPressureField* _NodalPressureField_New( NODALPRESSUREFIELD_ARGS );
void _NodalPressureField_Init( NodalPressureField* self,
			       Variable_Register* variable_Register );
void* _NodalPressureField_DefaultNew( Name name );
void _NodalPressureField_Delete( void* _self );

/*
** Methods.
*/

void _NodalPressureField_Print( void* _self, Stream* stream );
void* _NodalPressureField_Copy( void* _self, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

void _NodalPressureField_Construct( void* _self, Stg_ComponentFactory* cf, void* data ) ;
void _NodalPressureField_Build( void* _self, void* data ) ;
void _NodalPressureField_Initialise( void* _self, void* data ) ;
void _NodalPressureField_Execute( void* _self, void* data ) ;
void _NodalPressureField_Destroy( void* _self, void* data ) ;
void _NodalPressureField_ValueAtParticle( void* _self, IntegrationPointsSwarm* swarm,
					  Element_LocalIndex lElement_I, void* particle,
					  double* pressure ) ;

#endif
