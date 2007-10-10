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
** $Id: Decomposer.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Domain_Mesh_Decomposer_h__
#define __Domain_Mesh_Decomposer_h__

	/** Textual name of this class */
	extern const Type	Decomposer_Type;

	/** Virtual function types */
	typedef void (Decomposer_DecomposeFunc)( void* decomposer, unsigned nDomains, unsigned* domains, 
						 CommTopology** commTopology, Decomp** decomp, Decomp_Sync** sync );

	/** Class contents */
	#define __Decomposer					\
		/* General info */				\
		__Stg_Class					\
								\
		/* Virtual info */				\
		Decomposer_DecomposeFunc*	decomposeFunc;	\
								\
		/* Decomposer info */				\
		MPI_Comm		comm;

	struct Decomposer { __Decomposer };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define DECOMPOSER_DEFARGS				\
		STG_CLASS_DEFARGS,				\
		Decomposer_DecomposeFunc*	decomposeFunc

	#define DECOMPOSER_PASSARGS			\
		STG_CLASS_PASSARGS, decomposeFunc

	Decomposer* Decomposer_New();
	Decomposer* _Decomposer_New( DECOMPOSER_DEFARGS );
	void _Decomposer_Init( Decomposer* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _Decomposer_Delete( void* decomposer );
	void _Decomposer_Print( void* decomposer, Stream* stream );

	void _Decomposer_Decompose( void* decomposer, unsigned nDomains, unsigned* domains, 
				    CommTopology** commTopology, Decomp** decomp, Decomp_Sync** sync );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void Decomposer_SetComm( void* decomposer, MPI_Comm comm );

	#define Decomposer_Decompose( decomposer, nDomains, domains, commTopo, decomp, sync )				\
		(assert( (decomposer) && ((Decomposer*)decomposer)->decomposeFunc ),					\
		 ((Decomposer*)decomposer)->decomposeFunc( decomposer, nDomains, domains, commTopo, decomp, sync ))

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	void Decomposer_BuildCommTopology( Decomposer* self, unsigned nDomains, unsigned* domains, 
					   CommTopology** commTopo, RangeSet*** isects );
	void Decomposer_BuildLocalIntersections( Decomposer* self, unsigned nDomains, unsigned* domains, 
						 CommTopology* commTopo, RangeSet*** isects );
	void Decomposer_Claim( Decomposer* self, CommTopology* topo, RangeSet** isects, 
			       unsigned nDomains, unsigned* domains, 
			       Decomp** decomp, Decomp_Sync** sync );
	void Decomposer_BuildIndices( Decomposer* self, unsigned nDomains, unsigned* domains, RangeSet* claimed, 
				      CommTopology* commTopo, Decomp** decomp, Decomp_Sync** sync );

#endif /* __Domain_Mesh_Decomposer_h__ */
