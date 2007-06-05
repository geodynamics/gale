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
**		Public types for this module.
**
** Assumptions:
**
** Comments:
**
** $Id: types.h 859 2007-06-05 06:55:16Z DavidLee $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_Discretisation_types_h__
#define __StgFEM_Discretisation_types_h__
	
	/* FE types/classes */
	typedef struct ElementType               ElementType;
	typedef struct ElementType_Register      ElementType_Register;
	typedef struct ConstantElementType       ConstantElementType;
	typedef struct BilinearElementType       BilinearElementType;
	typedef struct TrilinearElementType      TrilinearElementType;
	typedef struct RegularTrilinear			RegularTrilinear;
	typedef struct RegularBilinear			RegularBilinear;
	typedef struct Biquadratic		Biquadratic;
	typedef struct P1			P1;
	typedef struct LinearTriangleElementType LinearTriangleElementType;
	typedef struct FiniteElement_Element     FiniteElement_Element;
	typedef struct FeMesh			FeMesh;
	typedef struct FeMesh_ElementType	FeMesh_ElementType;
	typedef struct C0Generator		C0Generator;
	typedef struct C2Generator		C2Generator;
	typedef struct P1Generator		P1Generator;
	typedef struct LinkedDofInfo             LinkedDofInfo;
	typedef struct FeEquationNumber          FeEquationNumber;
	typedef struct FeVariable                FeVariable;
	typedef struct ShapeFeVariable           ShapeFeVariable;
	typedef struct OperatorFeVariable        OperatorFeVariable;
	typedef struct FeSwarmVariable           FeSwarmVariable;
	typedef struct AnalyticSolution          AnalyticSolution;
	typedef struct Triquadratic              Triquadratic;

	/* types for lists etc ... for readability */
	typedef FeVariable*				FeVariablePtr;
	typedef Stg_ObjectList				FeVariableList;

	/* Types for clarity */
	typedef Index                   Iteration_Index;

	/* Stuff for FeEquationNumber class */
	typedef Node_GlobalIndex Node_RemappedGlobalIndex;
	typedef Index RemappedNodeInfo_Index;
	typedef int					Dof_EquationNumber;
	typedef Dof_EquationNumber*			Dof_EquationNumbers;
	typedef Dof_EquationNumbers*			Dof_EquationNumbersList;

	/* basic types for ElementRegister */
	typedef Index					ElementType_Index;

	/* output streams: initialised in Init() */
	/* Let's give this library responsiblity for setting up the StgFEM stream, as it's the first to be compiled...  */
	extern Stream* StgFEM_Debug;
	extern Stream* StgFEM_Warning;

	extern Stream* StgFEM_Discretisation_Debug;

#endif /* __StgFEM_Discretisation_types_h__ */
