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
**	Utilities to get the sides of a regular mesh etc.
**
** Assumptions:
**
** Comments:
**
** $Id: RegularMeshUtils.h 4081 2007-04-27 06:20:07Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Discretisation_Utils_RegularMeshUtils_h__
#define __Discretisation_Utils_RegularMeshUtils_h__

	extern Index RegularMeshUtils_ascendingIJK_ToHughesNodeNumberMap[8];

	/*--------------------------------------------------------------------------------------------------------------------------
	** Mapping functions
	*/

	void RegularMeshUtils_Node_1DTo3D( void* mesh, unsigned global, unsigned* inds );
	unsigned RegularMeshUtils_Node_3DTo1D( void* mesh, unsigned* inds );

	void RegularMeshUtils_Element_1DTo3D( void* mesh, unsigned global, unsigned* inds );
	unsigned RegularMeshUtils_Element_3DTo1D( void* mesh, unsigned* inds );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Set functions
	*/
	
	/** Create a new set, based on node indices, of nodes on the top of the global regular mesh */
	IndexSet* RegularMeshUtils_CreateGlobalTopSet( void* _mesh );
	
	/** Create a new set, based on node indices, of nodes on the bottom of the global regular mesh */
	IndexSet* RegularMeshUtils_CreateGlobalBottomSet( void* _mesh );
	
	/** Create a new set, based on node indices, of nodes on the left of the global regular mesh */
	IndexSet* RegularMeshUtils_CreateGlobalLeftSet( void* _mesh );
	
	/** Create a new set, based on node indices, of nodes on the right of the global regular mesh */
	IndexSet* RegularMeshUtils_CreateGlobalRightSet( void* _mesh );
	
	/** Create a new set, based on node indices, of nodes on the front of the global regular mesh */
	IndexSet* RegularMeshUtils_CreateGlobalFrontSet( void* _mesh );
	
	/** Create a new set, based on node indices, of nodes on the back of the global regular mesh */
	IndexSet* RegularMeshUtils_CreateGlobalBackSet( void* _mesh );


	/** Create a new set, based on node indices, of nodes on the top without the corners of the global regular mesh */
	IndexSet* RegularMeshUtils_CreateGlobalInnerTopSet( void* _mesh );
	
	/** Create a new set, based on node indices, of nodes on the bottom without the corners of the global regular mesh */
	IndexSet* RegularMeshUtils_CreateGlobalInnerBottomSet( void* _mesh );
	
	/** Create a new set, based on node indices, of nodes on the left without the corners of the global regular mesh */
	IndexSet* RegularMeshUtils_CreateGlobalInnerLeftSet( void* _mesh );
	
	/** Create a new set, based on node indices, of nodes on the right without the corners of the global regular mesh */
	IndexSet* RegularMeshUtils_CreateGlobalInnerRightSet( void* _mesh );
	
	/** Create a new set, based on node indices, of nodes on the front without the corners of the global regular mesh */
	IndexSet* RegularMeshUtils_CreateGlobalInnerFrontSet( void* _mesh );
	
	/** Create a new set, based on node indices, of nodes on the back without the corners of the global regular mesh */
	IndexSet* RegularMeshUtils_CreateGlobalInnerBackSet( void* _mesh );
	

	/** Create a new set, based on node indices, of the node on the bottom front lefthand corner */
	IndexSet* RegularMeshUtils_CreateGlobalBottomLeftFrontSet( void* _mesh );

	/** Create a new set, based on node indices, of the node on the bottom front righthand corner */
	IndexSet* RegularMeshUtils_CreateGlobalBottomRightFrontSet( void* _mesh );

	/** Create a new set, based on node indices, of the node on the top front lefthand corner */
	IndexSet* RegularMeshUtils_CreateGlobalTopLeftFrontSet( void* _mesh );

	/** Create a new set, based on node indices, of the node on the top front righthand corner */
	IndexSet* RegularMeshUtils_CreateGlobalTopRightFrontSet( void* _mesh );

	/** Create a new set, based on node indices, of the node on the bottom back lefthand corner */
	IndexSet* RegularMeshUtils_CreateGlobalBottomLeftBackSet( void* _mesh );

	/** Create a new set, based on node indices, of the node on the bottom back righthand corner */
	IndexSet* RegularMeshUtils_CreateGlobalBottomRightBackSet( void* _mesh );

	/** Create a new set, based on node indices, of the node on the top back lefthand corner */
	IndexSet* RegularMeshUtils_CreateGlobalTopLeftBackSet( void* _mesh );

	/** Create a new set, based on node indices, of the node on the top back righthand corner */
	IndexSet* RegularMeshUtils_CreateGlobalTopRightBackSet( void* _mesh );


	/** Create a new set, based on node indices, of nodes on the top of the local regular mesh */
	IndexSet* RegularMeshUtils_CreateLocalInGlobalTopSet( void* _mesh );
	
	/** Create a new set, based on node indices, of nodes on the bottom of the local regular mesh */
	IndexSet* RegularMeshUtils_CreateLocalInGlobalBottomSet( void* _mesh );
	
	/** Create a new set, based on node indices, of nodes on the left of the local regular mesh */
	IndexSet* RegularMeshUtils_CreateLocalInGlobalLeftSet( void* _mesh );
	
	/** Create a new set, based on node indices, of nodes on the right of the local regular mesh */
	IndexSet* RegularMeshUtils_CreateLocalInGlobalRightSet( void* _mesh );
	
	/** Create a new set, based on node indices, of nodes on the front of the local regular mesh */
	IndexSet* RegularMeshUtils_CreateLocalInGlobalFrontSet( void* _mesh );
	
	/** Create a new set, based on node indices, of nodes on the back of the local regular mesh */
	IndexSet* RegularMeshUtils_CreateLocalInGlobalBackSet( void* _mesh );

	Node_DomainIndex RegularMeshUtils_GetDiagOppositeAcrossElementNodeIndex( void* _mesh, 
										 Element_DomainIndex refElement_dI, 
										 Node_DomainIndex refNode_dI );

#endif /* __Discretisation_Utils_RegularMeshUtils_h__ */
