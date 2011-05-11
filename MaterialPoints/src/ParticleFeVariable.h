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
*/

#ifndef __PICellerator_MaterialPoints_ParticleFeVariable_h__
#define __PICellerator_MaterialPoints_ParticleFeVariable_h__

	/** @see ParticleFeVariable_ValueAtParticle */
	typedef void ParticleFeVariable_ValueAtParticleFunction( 
		void*							particleFeVariable, 
		IntegrationPointsSwarm*	swarm, 
		Element_LocalIndex		lElement_I, 
		void*							particle, 
		double*						particleValue );
	
	/** Textual name of this class */
	extern const Type ParticleFeVariable_Type;

	/** ParticleFeVariable class contents */
	#define __ParticleFeVariable \
		/* General info */ \
		__FeVariable \
		\
		/* Virtual info */ \
		ParticleFeVariable_ValueAtParticleFunction*	_valueAtParticle; \
		\
		/* ParticleFeVariable info */ \
		double*									data; \
		Variable*								dataVariable; \
		char*									assemblyVectorName; \
		ForceVector*								assemblyVector; \
		ForceTerm*								assemblyTerm; \
		char*									massMatrixName; \
		ForceVector*								massMatrix; \
		ForceTerm*								massMatrixForceTerm; \
		int									currentParticleIndex; \
		Bool									useDeriv; 
		
	struct ParticleFeVariable { __ParticleFeVariable };

	/* --- Contstructors / Destructors --- */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define PARTICLEFEVARIABLE_DEFARGS \
                FEVARIABLE_DEFARGS, \
                ParticleFeVariable_ValueAtParticleFunction*  _valueAtParticle

	#define PARTICLEFEVARIABLE_PASSARGS \
                FEVARIABLE_PASSARGS, \
	        _valueAtParticle

	ParticleFeVariable* _ParticleFeVariable_New(  PARTICLEFEVARIABLE_DEFARGS  );
	
	void _ParticleFeVariable_Init( ParticleFeVariable* self, IntegrationPointsSwarm* swarm );
	
	void _ParticleFeVariable_Delete( void* variable );

	void _ParticleFeVariable_Print( void* variable, Stream* stream );

	void* _ParticleFeVariable_Copy( const void* feVariable, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

	void _ParticleFeVariable_AssignFromXML( void* variable, Stg_ComponentFactory* cf, void* data );

	void _ParticleFeVariable_Build( void* variable, void* data );

	void _ParticleFeVariable_Initialise( void* variable, void* data );

	void _ParticleFeVariable_Execute( void* variable, void* data );

	void _ParticleFeVariable_Destroy( void* variable, void* data );

	/** Returns in particleValue the value represented on this particle */
	#define ParticleFeVariable_ValueAtParticle( particleFeVariable, swarm, lElement_I, particle, particleValue ) \
		( (ParticleFeVariable*)(particleFeVariable) )->_valueAtParticle(                                     \
			(particleFeVariable), \
			(swarm), \
			(lElement_I), \
			(particle), \
			(particleValue) )

	void ParticleFeVariable_Update( void* materialFeVariable );

	void ParticleFeVariable_AssembleElement( 
		void*						forceTerm, 
		ForceVector*			forceVector, 
		Element_LocalIndex	lElement_I, 
		double*					elForceVector );

	void ParticleFeVariable_AssembleElement_Deriv(
		void*						_forceTerm,
		ForceVector*			forceVector,
		Element_LocalIndex	lElement_I,
		double*					elForceVector ) ;

	void ParticleFeVariable_AssembleElementShapeFunc( 
		void*						forceTerm, 
		ForceVector*			forceVector, 
		Element_LocalIndex	lElement_I, 
		double*					elForceVector );

#endif 

