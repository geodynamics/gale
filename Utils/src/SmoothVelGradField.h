/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Monash Cluster Computing, Australia
** (C) 2003-2004 All Rights Reserved
**
** Primary Authors:
** Luke Hodkinson
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Underworld_Utils_SmoothVelGradField_h__
#define __Underworld_Utils_SmoothVelGradField_h__

/*
** Class definition.
*/

extern const Type SmoothVelGradField_Type;

#define __SmoothVelGradField			\
   __ParticleFeVariable				\
						\
   /* Virtual functions. */			\
						\
   /* Members. */				\
   Variable_Register* variable_Register;	\
   Variable* dataVariableList[9];		\
   FeVariable* velField;

struct SmoothVelGradField { __SmoothVelGradField };

#define SMOOTHVELGRADFIELD_ARGS			\
   PARTICLEFEVARIABLE_ARGS

#define SMOOTHVELGRADFIELD_PASSARGS		\
   PARTICLEFEVARIABLE_PASSARGS

/*
** Constructors/Destructors.
*/

SmoothVelGradField* _SmoothVelGradField_New( SMOOTHVELGRADFIELD_ARGS );
void _SmoothVelGradField_Init( SmoothVelGradField* self,
			       Variable_Register* variable_Register,
			       FeVariable* velField,
			       SystemLinearEquations* sle );
void* _SmoothVelGradField_DefaultNew( Name name );
void _SmoothVelGradField_Delete( void* _self );
void SmoothVelGradField_NonLinearUpdate( void* _sle, void* _ctx );
/*
** Methods.
*/

void _SmoothVelGradField_Print( void* _self, Stream* stream );
void* _SmoothVelGradField_Copy( void* _self, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

void _SmoothVelGradField_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data ) ;
void _SmoothVelGradField_Build( void* _self, void* data ) ;
void _SmoothVelGradField_Initialise( void* _self, void* data ) ;
void _SmoothVelGradField_Execute( void* _self, void* data ) ;
void _SmoothVelGradField_Destroy( void* _self, void* data ) ;
void _SmoothVelGradField_ValueAtParticle( void* _self, IntegrationPointsSwarm* swarm,
					  Element_LocalIndex lElement_I, void* particle,
					  double* pressure ) ;

#endif
