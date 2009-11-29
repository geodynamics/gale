/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**	Computational Infrastructure for Geodynamics - http://www.geodynamics.org
**
** Contributors:
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
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
**	Fe specific entry points.
**
** Assumptions:
**
** Comments:
**
** $Id: EntryPoint.h 733 2007-02-07 00:55:26Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __FiniteElement_FeEntryPoint_h__
#define __FiniteElement_FeEntryPoint_h__
	
	/* Templates for feEntry point type */
	typedef void		(FeEntryPoint_AssembleStiffnessMatrix_Function)	( 
					void* stiffnessMatrix, 
					Bool bcRemoveQuery,
					void* _sle,
					void* _context
					);
	typedef void		(FeEntryPoint_AssembleStiffnessMatrix_CallFunction)	( 
					void* feEntryPoint, 
					void* stiffnessMatrix, 
					Bool bcRemoveQuery,
					void* _sle,
					void* _context
					);

	typedef void		(FeEntryPoint_AssembleForceVector_Function)	( 
					void* forceVector );
	typedef void		(FeEntryPoint_AssembleForceVector_CallFunction)	( 
					void* feEntryPoint, 
					void* forceVector ); 
		
	#define 			FeEntryPoint_AssembleStiffnessMatrix_CastType	(EntryPoint_CastType_MAX+1)
	#define 			FeEntryPoint_AssembleForceVector_CastType	(FeEntryPoint_AssembleStiffnessMatrix_CastType+1)
	#define 			FeEntryPoint_CastType_MAX			(FeEntryPoint_AssembleForceVector_CastType+1)
	
	/** Textual name of this class */
	extern const Type FeEntryPoint_Type;
	
	/** FeEntryPoint info */
	#define __FeEntryPoint \
		/* General info */ \
		__EntryPoint \
		\
		/* Virtual info */ \
		\
		/* FeEntryPoint info */
	struct FeEntryPoint { __FeEntryPoint };
	
	/* Create a new FeEntryPoint */
	FeEntryPoint* FeEntryPoint_New( Name name, unsigned int castType );
	
	/* Initialise an FeEntryPoint */
	void FeEntryPoint_Init( void* feEntryPoint, Name name, unsigned int castType );
	
	/* Creation implementation */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define FEENTRYPOINT_DEFARGS \
                ENTRYPOINT_DEFARGS

	#define FEENTRYPOINT_PASSARGS \
                ENTRYPOINT_PASSARGS

	FeEntryPoint* _FeEntryPoint_New(  FEENTRYPOINT_DEFARGS  );
	
	/* Initialisation implementation */
	void _FeEntryPoint_Init( FeEntryPoint* self );
	
	/* Copy */
	#define FeEntryPoint_Copy( self ) \
		(FeEntryPoint*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define FeEntryPoint_DeepCopy( self ) \
		(FeEntryPoint*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
	
	/* modified GetRun implementation */
	Func_Ptr _FeEntryPoint_GetRun( void* feEntryPoint );
	
	void _FeEntryPoint_Run_AssembleStiffnessMatrix( 
		void* feEntryPoint, 
		void* stiffnessMatrix, 
		Bool bcRemoveQuery,
		void* _sle,
		void* _context
		);

	void _FeEntryPoint_Run_AssembleForceVector( 
		void* feEntryPoint, 
		void* forceVector );
	
#endif /* __Fe_FeEntryPoint_h__ */

