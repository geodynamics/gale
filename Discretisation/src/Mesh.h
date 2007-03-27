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
**	Brings together the geometry and node interconnectivity of a mesh, plus the element type (i.e. shape functions) needed by the finite element method.
**
** Assumptions:
**	Elements used with this mesh inherit from FiniteElement_Element.
**	The contents of the mesh's elements must be initialised by the user.
**
** Comments:
**	This mesh inherits from the parent Discretisation-level Mesh class, which defines how
**	the mesh is laid out, manages local to global mappings etc.
**	The extra complexity this class adds from an FE perspective is the relationship between
**	the geometric element /node types to the F.E. ElementType determining shape functions.
**
**	Situations involving this class become interesting when using multiple FeVariable objects
**	in a simulation. See the Twiki notes at 
**	csd.vpac.org/twiki/bin/view/Stgermain/FiniteElement_Mesh for more on this.
**
** $Id $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_Discretisation_FiniteElement_Mesh_h__
#define __StgFEM_Discretisation_FiniteElement_Mesh_h__
	
	/** Textual name of this class */
	extern const Type FiniteElement_Mesh_Type;
	
	/** FiniteElement_Mesh information */
	#define __FiniteElement_Mesh \
		/* General info */ \
		__Mesh \
		\
		/* Virtual info */ \
		\
		/* FiniteElement_Mesh info */ \
		/** The basic discretisation Mesh */ \
		ElementType_Register*			elementType_Register; \
		
	/** Brings together and manages the life-cycle of a a mesh and all the 
	info relevant to it for the Finite Element Method. See Mesh.h for more. */
	struct FiniteElement_Mesh { __FiniteElement_Mesh };
	
	/* --- Constructors / Destructor --- */
	
	/** Create a new FiniteElement_Mesh and initialise */

	void* FiniteElement_Mesh_DefaultNew( Name name );
	
	FiniteElement_Mesh* FiniteElement_Mesh_New(
		Name					name,
		void*					layout,
		SizeT					nodeSize,
		SizeT					elementSize,
		void*					extension_Register,
		void*					elementType_Register,
		Dictionary*				dictionary );

	void FiniteElement_Mesh_LoadFromDict( void* mesh, Dictionary* subDict, Dictionary* dictionary, Stg_ObjectList* objList);
	
	/* Initialise a FiniteElement_Mesh construct */
	void FiniteElement_Mesh_Init(
		FiniteElement_Mesh*			self,
		Name					name,
		void*					layout,
		SizeT					nodeSize,
		SizeT					elementSize,
		void*					extension_Register,
		void*					elementType_Register,
		Dictionary*				dictionary );
	
	/* Creation implementation / Virtual constructor */
	FiniteElement_Mesh* _FiniteElement_Mesh_New(
		SizeT					_sizeOfSelf,
		Type					type,
		Stg_Class_DeleteFunction*			_delete,
		Stg_Class_PrintFunction*			_print,
		Stg_Class_CopyFunction*			_copy,
		Stg_Component_DefaultConstructorFunction*	_defaultConstructor,
		Stg_Component_ConstructFunction*	_construct,
		Stg_Component_BuildFunction*		_build,
		Stg_Component_InitialiseFunction*		_initialise,
		Stg_Component_ExecuteFunction*		_execute,
		Stg_Component_DestroyFunction*		_destroy,
		Name					name,
		Bool					initFlag,
		Mesh_Node_IsLocalFunction*		_nodeIsLocal,
		Mesh_Node_IsShadowFunction*		_nodeIsShadow,
		Mesh_Element_IsLocalFunction*		_elementIsLocal,
		Mesh_Element_IsShadowFunction*		_elementIsShadow,
		void*					layout,
		SizeT					nodeSize,
		SizeT					elementSize, 
		void*					extension_Register,
		void*					elementType_Register,
		Dictionary*				dictionary );

	/* Initialise implementation */
	void _FiniteElement_Mesh_Init( 
		FiniteElement_Mesh*			self,
		void*					elementType_Register );
	
	/* Stg_Class_Delete a FiniteElement_Mesh construst */
	void _FiniteElement_Mesh_Delete( void* feMesh );

	/* --- Virtual Function implementations --- */
	
	/* Print the contents of an FiniteElement_Mesh construct */
	void _FiniteElement_Mesh_Print( void* feMesh, Stream* stream );
	
	/* Copy */
	#define FiniteElement_Mesh_Copy( self ) \
		(FiniteElement_Mesh*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define FiniteElement_Mesh_DeepCopy( self ) \
		(FiniteElement_Mesh*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
	
	void* _FiniteElement_Mesh_Copy( void* feMesh, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	/* Build implementation */
	void _FiniteElement_Mesh_Build( void* feMesh, void* data );
	
	/* Construct implementation */
	void _FiniteElement_Mesh_Construct( void* feMesh, Stg_ComponentFactory* cf, void* data );
	
	/* Initialisation implementation */
	void _FiniteElement_Mesh_Initialise( void* feMesh, void* data );
	
	/* Execution implementation */
	void _FiniteElement_Mesh_Execute( void* feMesh, void* data );
	
	/* Destruct Implementation */
	void _FiniteElement_Mesh_Destroy( void* feMesh, void* data );

	/* --- Public Functions --- */
	
	/* For use when Node is not yet a complete type */
	#define FiniteElement_Mesh_CoordAt( self, node_dI ) \
		Mesh_CoordAt( self, node_dI )
	double* _FiniteElement_Mesh_CoordAt( void* feMesh, Node_DomainIndex node_dI );
	
	/* For use when Node is not yet a complete type */
	#define FiniteElement_Mesh_NodeAt( self, node_dI ) \
		Mesh_NodeAt( self, node_dI )
	Node* _FiniteElement_Mesh_NodeAt( void* feMesh, Node_DomainIndex node_I );
	
	/* For use when Element is not yet a complete type */
	#define FiniteElement_Mesh_ElementAt( self, element_I ) \
		((FiniteElement_Element*)Mesh_ElementAt( self, element_I ))
	FiniteElement_Element* _FiniteElement_Mesh_ElementAt( void* feMesh, Element_DomainIndex element_I );
	
	/** Get the element type of an element */
	#define FiniteElement_Mesh_ElementTypeAt( self, el_I ) \
		(ElementType_Register_At( (self)->elementType_Register, FiniteElement_Mesh_ElementAt( self, el_I )->elementType_I ))
	ElementType* _FiniteElement_Mesh_ElementTypeAt( void* feMesh, Element_DomainIndex element_I );
	
	/* --- Private Functions --- */
	
	void _FiniteElement_Mesh_SetupElementTypes_Default( FiniteElement_Mesh* self );

	/* --- Public Functions --- */
	Variable* FiniteElement_Mesh_RegisterNodeCoordsAsVariables( void* feMesh, void* _variable_Register, Variable** variableList ) ;
	Variable* FiniteElement_Mesh_RegisterLocalNodeCoordsAsVariables( void* feMesh, void* _variable_Register, Variable** variableList );
	void FiniteElement_Mesh_CalcGlobalCoordFromLocalCoord( void* feMesh, Dimension_Index dim, Element_LocalIndex lElement_I, Coord xi, Coord globalCoord ) ;

	/** This is a higher-level function to the actual one on the elementType that does the work.
		It was decided the user should pass in the current element, since they usually know this when doing this
		conversion */
	void FiniteElement_Mesh_CalcLocalCoordFromGlobalCoord(
		void* feMesh,
		Element_DomainIndex elementCoordIn,
		const Coord globalCoord,
		Coord elLocalCoord );

	
#endif /* __StgFEM_Discretisation_FiniteElement_Mesh_h__ */
