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
** Invariants:
**
** Comments:
**
** $Id: MGCoarsener.h 2225 1970-01-02 13:48:23Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_plugins_MultiGrid_MGCoarsener_h__
#define __StgFEM_plugins_MultiGrid_MGCoarsener_h__
	

	/* Textual name of this class */
	extern const Type MGCoarsener_Type;

	/* Virtual function types */
	typedef void	(MGCoarsener_SetMeshFunc)( void* coarsener, void* mesh );
	typedef void	(MGCoarsener_CoarsenFunc)( void* coarsener, unsigned level, MGMapping* dstMap );
	
	/* Class contents */
	#define __MGCoarsener \
		/* General info */ \
		__Stg_Class \
		\
		/* Virtual info */ \
		MGCoarsener_SetMeshFunc*		_setMesh; \
		MGCoarsener_CoarsenFunc*		_coarsen; \
		\
		/* MGCoarsener info ... */ \
		Mesh*					mesh;

	struct MGCoarsener { __MGCoarsener };
	
	
	/*-----------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/
	
	/* Creation implementation */
	MGCoarsener* _MGCoarsener_New(
		SizeT					_sizeOfSelf, 
		Type						type,
		Stg_Class_DeleteFunction*	_delete,
		Stg_Class_PrintFunction*		_print, 
		Stg_Class_CopyFunction*		_copy, 
		MGCoarsener_SetMeshFunc*		_setMesh, 
		MGCoarsener_CoarsenFunc*		_coarsen );
	
	/* Initialisation implementation functions */
	void _MGCoarsener_Init( MGCoarsener* self );
	
	
	/*-----------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/
	
	/* Stg_Class_Delete implementation */
	void _MGCoarsener_Delete( void* coarsener );
	
	/* Print implementation */
	void _MGCoarsener_Print( void* coarsener, Stream* stream );
	
	/* Copy implementation */
	void* _MGCoarsener_Copy( 
		void*		coarsener, 
		void*		dest, 
		Bool			deep, 
		Name			nameExt, 
		PtrMap*		ptrMap );
	
	
	/*-----------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/
	
	#define MGCoarsener_SetMesh( self, mesh ) \
		((self)->_setMesh( self, mesh ))
	
	#define MGCoarsener_Coarsen( self, level, dstMap ) \
		((self)->_coarsen( self, level, dstMap ))
	
	
	/*-----------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/
	
	void _MGCoarsener_SetMesh( void* coarsener, void* mesh );
	
	
#endif /* __StgFEM_plugins_MultiGrid_MGCoarsener_h__ */
