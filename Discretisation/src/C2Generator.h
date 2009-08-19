/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
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
**
** Assumptions:
**
** Invariants:
**
** Comments:
**
** $Id: C2Generator.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_Discretisaton_C2Generator_h__
#define __StgFEM_Discretisaton_C2Generator_h__

	/** Textual name of this class */
	extern const Type C2Generator_Type;

	/** Virtual function types */

	/** C2Generator class contents */
	#define __C2Generator			\
		/* General info */		\
		__CartesianGenerator		\
						\
		/* Virtual info */		\
						\
		/* C2Generator info */		\

	struct C2Generator { __C2Generator };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define C2GENERATOR_DEFARGS \
		CARTESIANGENERATOR_DEFARGS

	#define C2GENERATOR_PASSARGS \
		CARTESIANGENERATOR_PASSARGS

	C2Generator* C2Generator_New( Name name );
	C2Generator* _C2Generator_New( C2GENERATOR_DEFARGS );
	void _C2Generator_Init( C2Generator* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _C2Generator_Delete( void* meshGenerator );
	void _C2Generator_Print( void* meshGenerator, Stream* stream );
	void _C2Generator_Construct( void* meshGenerator, Stg_ComponentFactory* cf, void* data );
	void _C2Generator_Build( void* meshGenerator, void* data );
	void _C2Generator_Initialise( void* meshGenerator, void* data );
	void _C2Generator_Execute( void* meshGenerator, void* data );
	void _C2Generator_Destroy( void* meshGenerator, void* data );

	void C2Generator_SetTopologyParams( void* meshGenerator, unsigned* sizes, 
					    unsigned maxDecompDims, unsigned* minDecomp, unsigned* maxDecomp );
	void C2Generator_GenElementVertexInc( void* meshGenerator, IGraph* topo, Grid*** grids );
	void C2Generator_GenFaceVertexInc( void* meshGenerator, IGraph* topo, Grid*** grids );
	void C2Generator_GenEdgeVertexInc( void* meshGenerator, IGraph* topo, Grid*** grids );
	void C2Generator_GenElementTypes( void* meshGenerator, Mesh* mesh );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

#endif /* __StgFEM_Discretisaton_C2Generator_h__ */
