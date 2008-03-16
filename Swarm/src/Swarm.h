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
**	External header file to libDomain_Swarm.
**
** Assumptions:
**	None so far.
**
** Comments:
**	None so far.
**
** $Id: Swarm.h 4171 2007-08-16 03:33:59Z DavidLee $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Domain_Swarm_h__
#define __Domain_Swarm_h__

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

	#include "types.h"
	#include "shortcuts.h"

	#include "ShadowInfo.h"
	#include "StandardParticle.h"
	#include "IntegrationPoint.h"
	#include "CellLayout.h"
	#include "ParticleLayout.h"
	#include "SwarmClass.h"
	#include "Swarm_Register.h"
	#include "ParticleCommHandler.h"
	#include "ParticleMovementHandler.h"
	#include "ParticleShadowSync.h"
	#include "SingleCellLayout.h"
	#include "ElementCellLayout.h"
	#include "SingleCellLayout.h"
	#include "TriSingleCellLayout.h"
	#include "Random.h"
	#include "PerCellParticleLayout.h"
	#include "RandomParticleLayout.h"
	#include "MeshParticleLayout.h"
	#include "GaussParticleLayout.h"
	#include "TriGaussParticleLayout.h"
	#include "GaussBorderParticleLayout.h"
	#include "GlobalParticleLayout.h"
	#include "LineParticleLayout.h"
	#include "FileParticleLayout.h"
	#include "SpaceFillerParticleLayout.h"
	#include "UnionParticleLayout.h"
	#include "WithinShapeParticleLayout.h"
	#include "PlaneParticleLayout.h"
	#include "ManualParticleLayout.h"
	#include "SwarmVariable_Register.h"
	#include "SwarmVariable.h"
	#include "OperatorSwarmVariable.h"
	#include "SwarmOutput.h"
	#include "SwarmDump.h"
	#include "SwarmShapeVC.h"
	#include "Init.h"
	#include "Finalise.h"

#endif /* __Domain_Swarm_h__ */
