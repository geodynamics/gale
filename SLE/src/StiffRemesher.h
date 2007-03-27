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
** $Id: StiffRemesher.h 2225 1970-01-02 13:48:23Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_SLE_StiffRemesher_h__
#define __StgFEM_SLE_StiffRemesher_h__

	/* Textual name of this class. */
	extern const Type StiffRemesher_Type;

	/* Virtual function types. */

	/* Class contents. */
	#define __StiffRemesher					\
		/* General info */				\
		__Remesher					\
								\
		/* Virtual info */				\
								\
		/* StiffRemesher info ... */			\
		unsigned		nDims;			\
		Swarm*			swarm;			\
		double***		weights;		\
								\
		Bool**			rests;			\
		unsigned		baseEqNum;		\
		unsigned**		eqNums;			\
								\
		Matrix*			stiffMat;		\
		Vector*			solVec;			\
		Vector*			rhsVec;			\
		MatrixSolver*		matSolver;

	struct StiffRemesher { __StiffRemesher };


	/*-----------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	/* Create a StiffRemesher */
	StiffRemesher* StiffRemesher_New( Name name );

	/* Creation implementation */
	StiffRemesher* _StiffRemesher_New( CLASS_ARGS, 
					   COMPONENT_ARGS, 
					   REMESHER_ARGS );

	/* Initialise a StiffRemesher */
	void StiffRemesher_Init( StiffRemesher* self );

	/* Initialisation implementation functions */
	void _StiffRemesher_Init( StiffRemesher* self );


	/*-----------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _StiffRemesher_Delete( void* stiffRemesher );
	void _StiffRemesher_Print( void* stiffRemesher, Stream* stream );
	StiffRemesher* _StiffRemesher_DefaultNew( Name name );
	void _StiffRemesher_Construct( void* stiffRemesher, Stg_ComponentFactory* cf, void* data );
	void _StiffRemesher_Build( void* stiffRemesher, void* data );
	void _StiffRemesher_Initialise( void* stiffRemesher, void* data );
	void _StiffRemesher_Execute( void* stiffRemesher, void* data );
	void _StiffRemesher_Destroy( void* stiffRemesher, void* data );

	void _StiffRemesher_SetMesh( void* stiffRemesher, Mesh* mesh );


	/*-----------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void StiffRemesher_SetRestriction( void* stiffRemesher, unsigned lNodeInd, unsigned dim, Bool state );
	void StiffRemesher_BuildSystem( void* stiffRemesher );


	/*-----------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	void _StiffRemesher_Free( StiffRemesher* self );
	void _StiffRemesher_FreeSystem( StiffRemesher* self );
	void _StiffRemesher_UpdateWeights( StiffRemesher* self );
	void _StiffRemesher_CalcWeights( StiffRemesher* self );

#endif
