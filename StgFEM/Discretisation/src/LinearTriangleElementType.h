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
**	Implementation of the FE Element basis functions for the standard Triangle coordinate system.
**
** Assumptions:
**
** Comments:
**
** $Id: LinearTriangleElementType.h 1177 2008-07-15 01:29:58Z DavidLee $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_Discretisation_LinearTriangleElementType_h__
#define __StgFEM_Discretisation_LinearTriangleElementType_h__
	
	/* Textual name of this class */
	extern const Type LinearTriangleElementType_Type;
	
	/* LinearTriangleElementType information */
	#define __LinearTriangleElementType \
		/* General info */ \
		__ElementType \
		\
		/* Virtual info */ \
		\
		/* LinearTriangleElementType info */

	struct LinearTriangleElementType { __LinearTriangleElementType };
	

	
	/* Create a new LinearTriangleElementType and initialise */
	void* _LinearTriangleElementType_DefaultNew( Name name );

	LinearTriangleElementType* LinearTriangleElementType_New( Name name );
	
	/* Creation implementation / Virtual constructor */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define LINEARTRIANGLEELEMENTTYPE_DEFARGS \
                ELEMENTTYPE_DEFARGS

	#define LINEARTRIANGLEELEMENTTYPE_PASSARGS \
                ELEMENTTYPE_PASSARGS

	LinearTriangleElementType* _LinearTriangleElementType_New(  LINEARTRIANGLEELEMENTTYPE_DEFARGS  );			
	
	/* Initialise implementation */
	void _LinearTriangleElementType_Init( LinearTriangleElementType* self );
	
	/* Stg_Class_Delete a LinearTriangleElementType construst */
	void _LinearTriangleElementType_Delete( void* elementType );
	
	/* Print the contents of an LinearTriangleElementType construct */
	void _LinearTriangleElementType_Print( void* elementType, Stream* stream );
	
	void _LinearTriangleElementType_AssignFromXML( void* elementType, Stg_ComponentFactory *cf, void* data );
	
	/* LinearTriangle element type build implementation */
	void _LinearTriangleElementType_Build( void* elementType, void *data );
	
	void _LinearTriangleElementType_Initialise( void* elementType, void *data );
	
	void _LinearTriangleElementType_Execute( void* elementType, void *data );
	
	void _LinearTriangleElementType_Destroy( void* elementType, void *data );
	
	/** ElementType_EvaluateShapeFunctionsAt implementation. */
	void _LinearTriangleElementType_SF_allNodes( void* elementType, const double localCoord[], double* const evaluatedValues );
	
	/** ElementType_EvaluateShapeFunctionsAt implementation. */
	void _LinearTriangleElementType_SF_allLocalDerivs_allNodes( void* elementType, const double localCoord[],
		double** const evaluatedDerivatives );
	
	int _LinearTriangularElementType_SurfaceNormal( void* elementType, unsigned element_I, unsigned dim, double* xi, double* normal );

#endif /* __StgFEM_Discretisation_LinearTriangleElementType_h__ */

