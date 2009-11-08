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
** Role:
**
** Assumptions:
**
** Comments:
**	The original sizes need to be manually set by the user.... this whole system needs rethinking... it can be done better
**
** $Id: ElementType_Register.h 656 2006-10-18 06:45:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_Discretisation_ElementType_Register_h__
#define __StgFEM_Discretisation_ElementType_Register_h__
	
	
	/* Textual name of this class */
	extern const Type ElementType_Register_Type;

	/* global default instantiation of this register (created in Init.c) */
	extern ElementType_Register* elementType_Register;
	
	/* ElementType_Register info */
	#define __ElementType_Register \
		/* General info */ \
		__Stg_Component \
		\
		DomainContext*			context; \
		/* Virtual info */ \
		\
		/* ElementType_Register info */ \
		Stream*					debug; \
		ElementType_Index		count; \
		SizeT						_size; \
		SizeT						_delta; \
		ElementType**			elementType;

	struct ElementType_Register { __ElementType_Register };

	#define ELEMENTTYPEREGISTER_DEFARGS \
		STG_COMPONENT_DEFARGS \

	#define ELEMENTTYPEREGISTER_PASSARGS \
  		STG_COMPONENT_PASSARGS \
	
	/* Create a new ElementType_Register */
	void* ElementType_Register_DefaultNew( Name name );
	
	ElementType_Register* ElementType_Register_New( Name name );
	
	/* Creation implementation / Virtual constructor */
	ElementType_Register* _ElementType_Register_New( ELEMENTTYPEREGISTER_DEFARGS );
	
	/* Initialisation implementation */
	void _ElementType_Register_Init( void* elementType_Register );
	
	/* Stg_Class_Delete implementation */
	void _ElementType_Register_Delete( void* elementType_Register );
	
	/* Print implementation */
	void _ElementType_Register_Print( void* elementType_Register, Stream* stream );
	
	void _ElementType_Register_AssignFromXML( void* elementType_Register, Stg_ComponentFactory *cf, void* data );
	
	void _ElementType_Register_Build( void* elementType_Register, void *data );
	
	void _ElementType_Register_Initialise( void* elementType_Register, void *data );
	
	void _ElementType_Register_Execute( void* elementType_Register, void *data );

	void _ElementType_Register_Destroy( void* elementType_Register, void *data );
	
	
	/* Add a new elementType */
	ElementType_Index ElementType_Register_Add( void* elementType_Register, void* elementType );
	
	/* Get the handle to an elementType */
	ElementType_Index ElementType_Register_GetIndex( void* elementType_Register, Type type );
	
	/* Get an element type from the register */
	#define ElementType_Register_At( elementType_Register, handle )		((elementType_Register)->elementType[(handle)] )
	ElementType* _ElementType_Register_At( void* elementType_Register, ElementType_Index handle );
	
#endif /* __StgFEM_Discretisation_ElementType_Register_h__ */
