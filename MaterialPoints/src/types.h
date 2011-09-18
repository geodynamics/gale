/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
** Copyright (c) 2005-2006, Monash Cluster Computing, Building 28, Monash University Clayton Campus,
**	Victoria, 3800, Australia
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**
** Contributors:
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Patrick D. Sunter, Software Engineer, VPAC. (patrick@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	David Stegman, Postdoctoral Fellow, Monash University. (david.stegman@sci.monash.edu.au)
**	Wendy Sharples, PhD Student, Monash University (wendy.sharples@sci.monash.edu.au)
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
**	Type definitions for this module.
**
** Assumptions:
**	None as yet.
**
** Comments:
**	None as yet.
**
** $Id: types.h 582 2008-08-12 03:15:11Z WendySharples $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __PICellerator_MaterialPoints_types_h__
#define __PICellerator_MaterialPoints_types_h__

	typedef struct MappedParticleLayout         MappedParticleLayout;
	typedef struct BackgroundParticleLayout     BackgroundParticleLayout;

	typedef struct IntegrationPointMapper       IntegrationPointMapper;
	typedef struct OneToOneMapper               OneToOneMapper;
        typedef struct OneToManyRef OneToManyRef;
        typedef struct OneToManyMapper OneToManyMapper;
        typedef struct NearestNeighborMapper NearestNeighborMapper;
	typedef struct CoincidentMapper             CoincidentMapper;
	typedef struct GaussCoincidentMapper        GaussCoincidentMapper;
	typedef struct GaussMapper                  GaussMapper;
	typedef struct ManyToOneMapper              ManyToOneMapper;

	typedef struct MaterialPointsSwarm           MaterialPointsSwarm;
	typedef struct IntegrationPointsSwarm       IntegrationPointsSwarm;
	typedef struct PICelleratorContext          PICelleratorContext;
	
	typedef struct Material                     Material;
	typedef struct Materials_Register           Materials_Register;
	typedef struct ParticleFeVariable           ParticleFeVariable;
	typedef struct MaterialFeVariable           MaterialFeVariable;
	typedef struct MaterialPoint                MaterialPoint;
	typedef struct MaterialPointRef             MaterialPointRef;
	typedef struct SwarmAdvector                SwarmAdvector;
	typedef struct SwarmAdvectionInAPlane                SwarmAdvectionInAPlane;

	typedef struct PeriodicBoundariesManager    PeriodicBoundariesManager;

	/* dave, 18.09.07 */
	typedef struct SwarmVariableField	SwarmVariableField;

	/* Types, for readability */
	typedef Index                             Material_Index;

#endif

