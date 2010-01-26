/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd,
** 110 Victoria Street, Melbourne, 3053, Australia.
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
**  Role:
**
** Assumptions:
**
** Invariants:
**
** Comments:
**
** $Id: MixedStabiliserTerm.h 2225 1970-01-02 13:48:23Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Experimental_Solvers_MixedStabiliserTerm_h__
#define __Experimental_Solvers_MixedStabiliserTerm_h__

extern const Type MixedStabiliserTerm_Type;

#define __MixedStabiliserTerm			\
   __StiffnessMatrixTerm                        \
   IntegrationPointsSwarm* picSwarm;          \
   StoreVisc* storeVisc;                 \
   double* Ni;                                  \
   double** elStiffMat;

struct MixedStabiliserTerm { __MixedStabiliserTerm };

	#ifndef ZERO
	#define ZERO 0
	#endif

	#define MIXEDSTABILISERTERM_DEFARGS \
                STIFFNESSMATRIXTERM_DEFARGS

	#define MIXEDSTABILISERTERM_PASSARGS \
                STIFFNESSMATRIXTERM_PASSARGS

        MixedStabiliserTerm* _MixedStabiliserTerm_New( MIXEDSTABILISERTERM_DEFARGS);

void _MixedStabiliserTerm_Init( MixedStabiliserTerm* self );

void _MixedStabiliserTerm_Delete( void* _self );

void _MixedStabiliserTerm_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data );
void _MixedStabiliserTerm_Build( void* _self, void* data );
void _MixedStabiliserTerm_Initialise( void* _self, void* data );
void _MixedStabiliserTerm_Execute( void* _self, void* data );
void _MixedStabiliserTerm_Destroy( void* _self, void* data );

void* _MixedStabiliserTerm_DefaultNew( Name name );

void _MixedStabiliserTerm_AssembleElement( void* _self,
					  StiffnessMatrix* stiffMat,
					  int elementIndex,
					  SystemLinearEquations* _sle,
					  FiniteElementContext* ctx,
					  double** elStiffMat );

#endif
