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
/** \file
**  Role:
**	Class allowing the user to construct and solve a hybrid FEM/PIC problem
**
** Assumptions:
**
** Comments:
**
** $Id: Context.h 374 2006-10-12 08:59:41Z SteveQuenette $
*
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __PICellerator_MaterialPoints_Context_h__
#define __PICellerator_MaterialPoints_Context_h__
	
	/* Textual name of this class */
	extern const Type PICelleratorContext_Type;
	
	#define __PICelleratorContext \
		/* General info */ \
		__FiniteElementContext \
		\
		/* Virtual info */ \
		\
		/* PICelleratorContext info */ \
		Materials_Register*	materials_Register;

	struct PICelleratorContext { __PICelleratorContext };


	
	/* Constructors ----------------------------------------------------------------------------------------------------*/
	
	/** Constructor */
	void* _PICelleratorContext_DefaultNew( Name name );
	
	PICelleratorContext* PICelleratorContext_New( 
		Name			name,
		double		start,
		double		stop,
		MPI_Comm		communicator,
		Dictionary*	dictionary );
	
	/** Creation implementation / Virtual constructor */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define PICELLERATORCONTEXT_DEFARGS \
                FINITEELEMENTCONTEXT_DEFARGS

	#define PICELLERATORCONTEXT_PASSARGS \
                FINITEELEMENTCONTEXT_PASSARGS

	PICelleratorContext* _PICelleratorContext_New(  PICELLERATORCONTEXT_DEFARGS  );
	
	/** Initialisation implementation */
	void _PICelleratorContext_Init( void* context );

	/* Virtual Functions -----------------------------------------------------------------------------------------------*/

	void _PICelleratorContext_AssignFromXML( void* context, Stg_ComponentFactory* cf, void* data );
	
	/* Stg_Class_Delete implementation */
	void _PICelleratorContext_Delete( void* context );

	/* Destroy implementation */	
	void _PICelleratorContext_Destroy( void* context );

	/* Print implementation */
	void _PICelleratorContext_Print( void* context, Stream* stream );
	
	/* Set the dt */
	void _PICelleratorContext_SetDt( void* context, double dt );

	void _PICelleratorContext_AssignFromXML( void* context, Stg_ComponentFactory* cf, void* data );
	
	/* Public functions -----------------------------------------------------------------------------------------------*/
	void PICelleratorContext_CreateDefaultMaterial( void* context ) ;
	
	/* Private functions -----------------------------------------------------------------------------------------------*/
#endif /* __PICelleratorContext_h__ */

