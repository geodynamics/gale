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
**
** Assumptions:
**
** Comments:
**
** $Id: TrilinearElementType.h 1178 2008-07-15 04:12:09Z DavidLee $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_Discretisation_TrilinearElementType_h__
#define __StgFEM_Discretisation_TrilinearElementType_h__
	
	/* Textual name of this class */
	extern const Type TrilinearElementType_Type;
	
	/* TrilinearElementType information */
	#define __TrilinearElementType \
		/* General info */ \
		__ElementType \
		\
		/* Virtual info */ \
		\
		/* TrilinearElementType info */ \
		Coord	minElLocalCoord; /** Bottom corner in elLocal mathematical space */ \
		Coord	maxElLocalCoord; /** Top corner in elLocal mathematical space */ \
		double	elLocalLength[3]; /** Length of element in elLocal space */ \
		\
		unsigned**	tetInds;
		
	struct TrilinearElementType { __TrilinearElementType };

	#define TRILINEARELEMENTTYPE_DEFARGS \
    	ELEMENTTYPE_DEFARGS

	#define TRILINEARELEMENTTYPE_PASSARGS \
    	ELEMENTTYPE_PASSARGS
	
	/* Create a new TrilinearElementType and initialise */
	void* _TrilinearElementTypea_DefaultNew( Name name );

	TrilinearElementType* TrilinearElementType_New( Name name );
	
	/* Creation implementation / Virtual constructor */
	TrilinearElementType* _TrilinearElementType_New( TRILINEARELEMENTTYPE_DEFARGS );
	
	/* Initialise implementation */
	void _TrilinearElementType_Init( TrilinearElementType* self );
	
	/* Stg_Class_Delete a TrilinearElementType construst */
	void _TrilinearElementType_Delete( void* elementType );
	
	/* Print the contents of an TrilinearElementType construct */
	void _TrilinearElementType_Print( void* elementType, Stream* stream );
	
	/* Trilinear element type build implementation */
	void _TrilinearElementType_AssignFromXML( void* elementType, Stg_ComponentFactory *cf, void* data );
	
	void _TrilinearElementType_Build( void* elementType, void *data );
	
	void _TrilinearElementType_Initialise( void* elementType, void *data );
	
	void _TrilinearElementType_Execute( void* elementType, void *data );
	
	void _TrilinearElementType_Destroy( void* elementType, void *data );
	
	/** ElementType_EvaluateShapeFunctionsAt implementation. */
	void _TrilinearElementType_SF_allNodes( void* elementType, const double localCoord[], double* const evaluatedValues );
	
	/** ElementType_EvaluateShapeFunctionsDerivsAt implementation. */
	void _TrilinearElementType_SF_allLocalDerivs_allNodes( void* elementType, const double localCoord[],
		double** const evaluatedDerivatives );
	
	/** ElementType_ConvertGlobalCoordToElLocal implementation.
	Uses a shortcut approach if using "box" elements - otherwise uses the general function. */
	void _TrilinearElementType_ConvertGlobalCoordToElLocal(
		void*		elementType,
		void*		mesh, 
		unsigned	element, 
		const double*	globalCoord,
		double*		elLocalCoord );
	
	double _TrilinearElementType_JacobianDeterminantSurface(
		void*		elementType,
		void*		mesh,
		unsigned	element_I,
		const double	localCoord[],
		unsigned	face_I,
		unsigned	norm );

#endif /* __StgFEM_Discretisation_TrilinearElementType_h__ */
