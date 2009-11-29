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
** $Id: SROpGenerator.h 672 2006-12-14 00:58:36Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Experimental_Solvers_SROpGenerator_h__
#define __Experimental_Solvers_SROpGenerator_h__

	/** Textual name of this class */
	extern const Type SROpGenerator_Type;

	/** Virtual function types */

	/** SROpGenerator class contents */
	#define __SROpGenerator				\
		/* General info */			\
		__MGOpGenerator				\
							\
		/* Virtual info */			\
							\
		/* SROpGenerator info */		\
		FeVariable*		fineVar;	\
		FeEquationNumber*	fineEqNum;	\
		Mesh**			meshes;		\
		unsigned**		topMaps;	\
		unsigned***		eqNums;		\
		unsigned*		nLocalEqNums;	\
		unsigned*		eqNumBases;

	struct SROpGenerator { __SROpGenerator };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/



	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define SROPGENERATOR_DEFARGS \
                MGOPGENERATOR_DEFARGS

	#define SROPGENERATOR_PASSARGS \
                MGOPGENERATOR_PASSARGS

	SROpGenerator* SROpGenerator_New( Name name );
	SROpGenerator* _SROpGenerator_New(  SROPGENERATOR_DEFARGS  );
	void _SROpGenerator_Init( SROpGenerator* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _SROpGenerator_Delete( void* srOpGenerator );
	void _SROpGenerator_Print( void* srOpGenerator, Stream* stream );
	void _SROpGenerator_AssignFromXML( void* srOpGenerator, Stg_ComponentFactory* cf, void* data );
	void _SROpGenerator_Build( void* srOpGenerator, void* data );
	void _SROpGenerator_Initialise( void* srOpGenerator, void* data );
	void _SROpGenerator_Execute( void* srOpGenerator, void* data );
	void _SROpGenerator_Destroy( void* srOpGenerator, void* data );

	Bool SROpGenerator_HasExpired( void* srOpGenerator );
	//void SROpGenerator_Generate( void* srOpGenerator, Matrix*** pOps, Matrix*** rOps );
	void SROpGenerator_Generate( void* srOpGenerator, Mat** pOps, Mat** rOps );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void SROpGenerator_SetFineVariable( void* srOpGenerator, void* variable );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	void SROpGenerator_GenMeshes( SROpGenerator* self );
	void SROpGenerator_GenLevelMesh( SROpGenerator* self, unsigned level );
	void SROpGenerator_GenLevelTopMap( SROpGenerator* self, unsigned level );
	void SROpGenerator_GenLevelEqNums( SROpGenerator* self, unsigned level );

	//void SROpGenerator_GenOps( SROpGenerator* self, Matrix** pOps, Matrix** rOps );
	void SROpGenerator_GenOps( SROpGenerator* self, Mat* pOps, Mat* rOps );
	//void SROpGenerator_GenLevelOp( SROpGenerator* self, unsigned level, Matrix* P );
	void SROpGenerator_GenLevelOp( SROpGenerator* self, unsigned level, Mat P );
	void SROpGenerator_CalcOpNonZeros( SROpGenerator* self, unsigned level, 
					   unsigned** nDiagNonZeros, unsigned** nOffDiagNonZeros );
	void SROpGenerator_DestructMeshes( SROpGenerator* self );
	//void SROpGenerator_Simple( SROpGenerator *self, Matrix **pOps, Matrix **rOps );
	void SROpGenerator_Simple( SROpGenerator *self, Mat* pOps, Mat* rOps );
	//Matrix *SROpGenerator_SimpleFinestLevel( SROpGenerator *self );
	//Matrix *SROpGenerator_SimpleCoarserLevel( SROpGenerator *self, int level );
	Mat SROpGenerator_SimpleFinestLevel( SROpGenerator *self );
	Mat SROpGenerator_SimpleCoarserLevel( SROpGenerator *self, int level );

#endif /* __Experimental_Solvers_SROpGenerator_h__ */

