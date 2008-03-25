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
** $Id: types.h 1084 2008-03-25 00:04:14Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_SLE_LinearAlgrebra_types_h__
#define __StgFEM_SLE_LinearAlgrebra_types_h__
	
	/* Classes. */
	typedef struct Vector			Vector;
	typedef struct Matrix			Matrix;
	typedef struct MatrixSolver		MatrixSolver;
	typedef struct NonlinearSolver		NonlinearSolver;
#ifdef HAVE_PETSC
	typedef struct PETScVector		PETScVector;
	typedef struct PETScMatrix		PETScMatrix;
	typedef struct PETScMatrixSolver	PETScMatrixSolver;
	typedef struct PETScNonlinearSolver	PETScNonlinearSolver;
#endif

	typedef enum {
		MatrixSolver_Status_ConvergedRelative = 2, 
		MatrixSolver_Status_ConvergedAbsolute = 3, 
		MatrixSolver_Status_ConvergedIterations = 4, 
		MatrixSolver_Status_DivergedNull = -2, 
		MatrixSolver_Status_DivergedIterations = -3, 
		MatrixSolver_Status_DivergedTolerance = -4, 
		MatrixSolver_Status_Iterating = 0
	} MatrixSolver_Status;

	/* KSP enumerations. */
	typedef enum {
		PETScMatrixSolver_KSPType_Richardson, 
		PETScMatrixSolver_KSPType_GMRes, 
		PETScMatrixSolver_KSPType_FGMRes, 
		PETScMatrixSolver_KSPType_CG, 
		PETScMatrixSolver_KSPType_PreOnly
	} PETScMatrixSolver_KSPType;

	/* PC enumerations. */
	typedef enum {
		PETScMatrixSolver_PCType_Jacobi, 
		PETScMatrixSolver_PCType_BlockJacobi, 
		PETScMatrixSolver_PCType_SOR, 
		PETScMatrixSolver_PCType_ParallelSOR, 
		PETScMatrixSolver_PCType_LU, 
		PETScMatrixSolver_PCType_RedundantLU, 
		PETScMatrixSolver_PCType_ILU, 
		PETScMatrixSolver_PCType_Multigrid, 
		PETScMatrixSolver_PCType_None
	} PETScMatrixSolver_PCType;

	/* Norm enumerations. */
	typedef enum {
		PETScMatrixSolver_NormType_None = 0, 
		PETScMatrixSolver_NormType_Preconditioned = 1, 
		PETScMatrixSolver_NormType_Unpreconditioned = 2, 
		PETScMatrixSolver_NormType_Natural = 3
	} PETScMatrixSolver_NormType;


#endif /* __StgFEM_SLE_LinearAlgrebra_types_h__ */

