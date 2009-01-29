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
** $Id: TrilinearInnerElType.h 822 2007-04-27 06:20:35Z DaveLee $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_Discretisation_TrilinearInnerElType_h__
#define __StgFEM_Discretisation_TrilinearInnerElType_h__
	
	/* Textual name of this class */
	extern const Type TrilinearInnerElType_Type;
	
	/* TrilinearInnerElType information */
	#define __TrilinearInnerElType \
		/* General info */ \
		__ElementType \
		\
		/* Virtual info */ \
		\
		/* TrilinearInnerElType info */ \
		Coord	minElLocalCoord; /** Bottom corner in elLocal mathematical space */ \
		Coord	maxElLocalCoord; /** Top corner in elLocal mathematical space */ \
		double	elLocalLength[3]; /** Length of element in elLocal space */ \
		\
		unsigned**	tetInds;
		
	struct TrilinearInnerElType { __TrilinearInnerElType };
	
	
	/* Create a new TrilinearInnerElType and initialise */
	void* TrilinearInnerElTypea_DefaultNew( Name name );

	TrilinearInnerElType* TrilinearInnerElType_New( Name name );
	
	/* Creation implementation / Virtual constructor */
	TrilinearInnerElType* _TrilinearInnerElType_New(
		SizeT								_sizeOfSelf,
		Type								type,
		Stg_Class_DeleteFunction*					_delete,
		Stg_Class_PrintFunction*					_print,
		Stg_Class_CopyFunction*						_copy, 
		Stg_Component_DefaultConstructorFunction*			_defaultConstructor,
		Stg_Component_ConstructFunction*				_construct,
		Stg_Component_BuildFunction*					_build,
		Stg_Component_InitialiseFunction*				_initialise,
		Stg_Component_ExecuteFunction*					_execute,
		Stg_Component_DestroyFunction*					_destroy,
		Name								name,
		Bool								initFlag,
		ElementType_EvaluateShapeFunctionsAtFunction*			_evaluateShapeFunctionsAt,
		ElementType_EvaluateShapeFunctionLocalDerivsAtFunction*		_evaluateShapeFunctionLocalDerivsAt,
		ElementType_ConvertGlobalCoordToElLocalFunction*		_convertGlobalCoordToElLocal,
		ElementType_JacobianDeterminantSurfaceFunction*			_jacobianDeterminantSurface,
		ElementType_SurfaceNormalFunction*				_surfaceNormal,
		Index								nodeCount );
	
	/* Initialise implementation */
	void _TrilinearInnerElType_Init( TrilinearInnerElType* self );
	
	/* Stg_Class_Delete a TrilinearInnerElType construst */
	void _TrilinearInnerElType_Delete( void* elementType );
	
	/* Print the contents of an TrilinearInnerElType construct */
	void _TrilinearInnerElType_Print( void* elementType, Stream* stream );
	
	/* Trilinear inner element type build implementation */
	void _TrilinearInnerElType_Construct( void* elementType, Stg_ComponentFactory *cf, void* data );
	
	void _TrilinearInnerElType_Build( void* elementType, void *data );
	
	void _TrilinearInnerElType_Initialise( void* elementType, void *data );
	
	void _TrilinearInnerElType_Execute( void* elementType, void *data );
	
	void _TrilinearInnerElType_Destroy( void* elementType, void *data );
	
	/** ElementType_EvaluateShapeFunctionsAt implementation. */
	void _TrilinearInnerElType_SF_allNodes( void* elementType, const double localCoord[], double* const evaluatedValues );
	
	/** ElementType_EvaluateShapeFunctionsDerivsAt implementation. */
	void _TrilinearInnerElType_SF_allLocalDerivs_allNodes( void* elementType, const double localCoord[],
		double** const evaluatedDerivatives );
	
	/** ElementType_ConvertGlobalCoordToElLocal implementation.
	Uses a shortcut approach if using "box" elements - otherwise uses the general function. */
	/*void _TrilinearInnerElType_ConvertGlobalCoordToElLocal(
		void*		elementType,
		void*		mesh, 
		unsigned	element, 
		const double*	globalCoord,
		double*		elLocalCoord );
	*/

	void _TrilinearInnerElType_SurfaceNormal( void* elementType, unsigned element_I, unsigned dim, double* xi, double* normal );

#endif /* __StgFEM_Discretisation_TrilinearInnerElType_h__ */
