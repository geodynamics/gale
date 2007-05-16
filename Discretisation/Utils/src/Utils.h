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
**	External header file to this library.
**
** Assumptions:
**	None so far.
**
** Comments:
**	None so far.
**
** $Id: Utils.h 4103 2007-05-16 01:09:50Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Discretisation_Utils_h__
#define __Discretisation_Utils_h__

	#include <assert.h>
	#include "types.h"
	#include "RegularMeshUtils.h"
	#include "AllElementsVC.h"
	#include "AllNodesVC.h"
	#include "WallVC.h"
	#include "CornerVC.h"
	#include "InnerWallVC.h"
	#include "ShapeVC.h"
	#include "FrictionVC.h"
	#include "SplitFrictionWallVC.h"
	#include "DofLayout.h"
	#include "Operator.h"
	#include "FieldVariable_Register.h"
	#include "FieldVariable.h"
	#include "OperatorFieldVariable.h"
	#include "DiscretisationContext.h"
	#include "LinearRegression.h"
	#include "SobolGenerator.h"
	#include "Remesher.h"
/*
	#include "StripRemesher.h"
	#include "CellRemesher.h"
*/
	
	#include "TimeIntegratee.h"
	#include "TimeIntegrator.h"
	#include "ShapeAdvector.h"

/*
	#include "SemiRegDeform.h"
*/
	#include "Init.h"
	#include "Finalise.h"
	
#endif /* __Discretisation_Utils_h__ */
