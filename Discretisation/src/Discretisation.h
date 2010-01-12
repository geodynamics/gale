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
**	External header file to this library/module
**
** Assumptions:
**	None so far.
**
** Comments:
**	None so far.
**
** $Id: Discretisation.h 1218 2008-09-04 06:18:44Z DavidLee $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_Discretisation_Discretisation_h__
#define __StgFEM_Discretisation_Discretisation_h__
	
	#include "units.h"
	#include "types.h"

	#include <petsc.h>
	#include <petscvec.h>
	#include <petscmat.h>
	#include <petscksp.h>
	#include <petscmg.h>
	#include <petscsnes.h>
	#include "PETScErrorChecking.h"

	#include "FeMesh_Algorithms.h"
	#include "FeMesh_ElementType.h"
	#include "ElementType.h"
	#include "ElementType_Register.h"
	#include "ConstantElementType.h"
	#include "BilinearElementType.h"
	#include "TrilinearElementType.h"
	#include "RegularTrilinear.h"
	#include "RegularBilinear.h"
	#include "Biquadratic.h"
	#include "Triquadratic.h"
	#include "P1.h"
	#include "LinearTriangleElementType.h"
	#include "BilinearInnerElType.h"
	#include "TrilinearInnerElType.h"

	#include "Element.h"
	#include "FeMesh.h"
	#include "C0Generator.h"
	#include "C2Generator.h"
	#include "P1Generator.h"
	#include "Inner2DGenerator.h"
	#include "LinkedDofInfo.h"
	#include "FeEquationNumber.h"
	#include "FeVariable.h"
	#include "ShapeFeVariable.h"
	#include "OperatorFeVariable.h"
	#include "FeSwarmVariable.h"
	#include "AnalyticSolution.h"
	#include "FieldTest.h"

	#include "FunctionSuite.h"

	#include "Init.h"
	#include "Finalise.h"

#endif /* __StgFEM_Discretisation_Discretisation_h__ */
