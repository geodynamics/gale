/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student Monash University, VPAC. (davidm@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
*/
/** \file
**  Role:
**	Finite Element Constitutive Matrix object.
**
** Assumptions:
**
** Comments:
**
** $Id: TimeIntegrand.h 3851 2006-10-12 08:57:22Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgDomain_Utils_TimeIntegrand_h__
#define __StgDomain_Utils_TimeIntegrand_h__
	
	typedef Bool (TimeIntegrand_CalculateTimeDerivFunction) ( void* timeIntegrator, Index array_I, double* timeDeriv );
	typedef void (TimeIntegrand_IntermediateFunction) ( void* timeIntegrator, Index array_I );

	extern const Type TimeIntegrand_Type;
	
	/* TimeIntegrand information */
	#define __TimeIntegrand  \
		/* General info */ \
		__Stg_Component \
		\
		DomainContext*				   context;		 \
		/* Virtual info */ \
		TimeIntegrand_CalculateTimeDerivFunction* _calculateTimeDeriv;  \
		TimeIntegrand_IntermediateFunction*       _intermediate;  \
		/* Other info */ \
		TimeIntegrator*                            timeIntegrator;       \
		Variable*                                  variable;             \
		Index                                      dataCount;            \
		Stg_Component**                            data;                 \
		Bool                                       allowFallbackToFirstOrder; \
		Stream*                                    debug;                \
		
	struct TimeIntegrand { __TimeIntegrand };
	
	/* Creation implementation / Virtual constructor */
	TimeIntegrand* TimeIntegrand_New( 
		Name                                   name,
		DomainContext*                         context,
		TimeIntegrator*                        timeIntegrator, 
		Variable*                              variable, 
		Index                                  dataCount, 
		Stg_Component**                        data,
		Bool                                   allowFallbackToFirstOrder );

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define TIMEINTEGRAND_DEFARGS \
                STG_COMPONENT_DEFARGS, \
                TimeIntegrand_CalculateTimeDerivFunction*  _calculateTimeDeriv, \
                TimeIntegrand_IntermediateFunction*              _intermediate

	#define TIMEINTEGRAND_PASSARGS \
                STG_COMPONENT_PASSARGS, \
	        _calculateTimeDeriv, \
	        _intermediate      

	TimeIntegrand* _TimeIntegrand_New(  TIMEINTEGRAND_DEFARGS  );

void _TimeIntegrand_Init( 
		void*                                      timeIntegrand,
		DomainContext*                             context,
		TimeIntegrator*                            timeIntegrator, 
		Variable*                                  variable, 
		Index                                      dataCount, 
		Stg_Component**                            data,
		Bool                                       allowFallbackToFirstOrder );
		
	/* 'Class' Virtual Functions */
	void _TimeIntegrand_Delete( void* timeIntegrator );
	void _TimeIntegrand_Print( void* timeIntegrator, Stream* stream );
	#define TimeIntegrand_Copy( self ) \
		(TimeIntegrand*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define TimeIntegrand_DeepCopy( self ) \
		(TimeIntegrand*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _TimeIntegrand_Copy( void* timeIntegrator, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	/* 'Stg_Component' Virtual Functions */
	void* _TimeIntegrand_DefaultNew( Name name ) ;
	void _TimeIntegrand_AssignFromXML( void* timeIntegrand, Stg_ComponentFactory* cf, void* data ) ;
	void _TimeIntegrand_Build( void* timeIntegrator, void* data );
	void _TimeIntegrand_Initialise( void* timeIntegrator, void* data );
	void _TimeIntegrand_Execute( void* timeIntegrator, void* data );
	void _TimeIntegrand_Destroy( void* timeIntegrand, void* data );

	/* +++ Virtual Functions +++ */
	#define TimeIntegrand_CalculateTimeDeriv( timeIntegrand, array_I, timeDeriv ) \
		( ((TimeIntegrand*) timeIntegrand )->_calculateTimeDeriv( timeIntegrand, array_I, timeDeriv ) )
	#define TimeIntegrand_Intermediate( timeIntegrand, array_I ) \
		( ((TimeIntegrand*) timeIntegrand )->_intermediate( timeIntegrand, array_I ) )

	/* +++ Private Functions +++ */
	Bool _TimeIntegrand_AdvectionTimeDeriv( void* timeIntegrand, Index array_I, double* timeDeriv ) ;
	void _TimeIntegrand_Intermediate( void* timeIntegrand, Index array_I );
	void _TimeIntegrand_RewindToStartAndApplyFirstOrderUpdate( 
		TimeIntegrand* self,
		double*         arrayDataPtr,
		double*         startData,
		double          startTime,
		double          dt,
		double*         timeDeriv,
		Index           array_I );

	/* +++ Public Functions +++ */
	void TimeIntegrand_FirstOrder( void* timeIntegrator, Variable* startValue, double dt );
	void TimeIntegrand_SecondOrder( void* timeIntegrator, Variable* startValue, double dt );
	void TimeIntegrand_FourthOrder( void* timeIntegrator, Variable* startValue, double dt );

	void TimeIntegrand_StoreTimeDeriv( void* timeIntegrand, Variable* timeDeriv ) ;
	void TimeIntegrand_Add2TimesTimeDeriv( void* timeIntegrand, Variable* timeDerivVariable ) ;
	void TimeIntegrand_FourthOrderFinalStep( void* timeIntegrand, Variable* startData, Variable* timeDerivVariable, double dt ) ;

	#define TimeIntegrand_GetTime( timeIntegrand ) \
		TimeIntegrator_GetTime( ((TimeIntegrand*) timeIntegrand)->timeIntegrator ) 

#endif 

