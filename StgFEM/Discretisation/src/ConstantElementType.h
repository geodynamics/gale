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
**	Defines the basis functions for a "constant" (discontinuous)
**	FE Element.
**
** Assumptions:
**
** Comments:
**
** $Id: ConstantElementType.h 1177 2008-07-15 01:29:58Z DavidLee $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_Discretisation_ConstantElementType_h__
#define __StgFEM_Discretisation_ConstantElementType_h__
	
	/* Textual name of this class */
	extern const Type ConstantElementType_Type;
	
	/* ConstantElementType information */
	#define __ConstantElementType \
		/* General info */ \
		__ElementType \
		\
		/* Virtual info */ \
		\
		/* ConstantElementType info */
	struct ConstantElementType { __ConstantElementType };



	
	void* ConstantElementType_DefaultNew( Name name );

	/* Create a new ConstantElementType and initialise */
	ConstantElementType* ConstantElementType_New( Name name );
	
	/* Creation implementation / Virtual constructor */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define CONSTANTELEMENTTYPE_DEFARGS \
                ELEMENTTYPE_DEFARGS

	#define CONSTANTELEMENTTYPE_PASSARGS \
                ELEMENTTYPE_PASSARGS

	ConstantElementType* _ConstantElementType_New(  CONSTANTELEMENTTYPE_DEFARGS  );
	
	/* Initialise a ConstantElementType construct */
	void _ConstantElementType_Init( ConstantElementType* self );
	
	/* Stg_Class_Delete a ConstantElementType construst */
	void _ConstantElementType_Delete( void* elementType );
	
	/* Print the contents of an ConstantElementType construct */
	void _ConstantElementType_Print( void* elementType, Stream* stream );
	
	/** Constant element type build implementation */
	
	void _ConstantElementType_AssignFromXML( void* elementType, Stg_ComponentFactory *cf, void* data );
	
	void _ConstantElementType_Build( void* elementType, void *data );
	
	void _ConstantElementType_Initialise( void* elementType, void *data );
	
	void _ConstantElementType_Execute( void* elementType, void *data );
	
	void _ConstantElementType_Destroy( void* elementType, void *data );
	
	/** Calculate the shape function for all degrees of freedom for all nodes */
	void _ConstantElementType_SF_allNodes( void* elementType, const double localCoord[], double* const evaluatedValues );
	
	/** Calculate the shape function local derivatives for all degrees of freedom for all nodes */
	void _ConstantElementType_SF_allLocalDerivs_allNodes( void* elementType, const double localCoord[],
		double** const evaluatedDerivatives );
	
	/** Constant element type implementation of ElementType_ConvertGlobalCoordToElLocal().
	Nothing to do... the "constant element" doesn't really have local coordinates
	since the shape function is independent of Xi and Eta ... just return a 
	centroid {0,0,0} which is consistent with the centre of the co-ordinate system
	for the bilinear, trilinear etc element types.
	*/
	void _ConstantElementType_ConvertGlobalCoordToElLocal(
		void*		elementType,
		void*		mesh, 
		unsigned	element, 
		const double*	globalCoord,
		double*		elLocalCoord );

	double _ConstantElementType_JacobianDeterminantSurface(
		void*		elementType,
		void*		mesh,
		const double	localCoord[],
		unsigned*	nodes,
		unsigned	norm );

	int _ConstantElementType_SurfaceNormal(
		void*		elementType,
		unsigned	element_I,
		unsigned	dim,
		double*		xi,
		double*		normal );

#endif /* __StgFEM_Discretisation_ConstantElementType_h__ */

