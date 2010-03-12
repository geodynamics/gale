/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Victorian Partnership for Advanced Computing (VPAC) Ltd, Australia
** (C) 2003-2005 All Rights Reserved
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	David May, PhD Student Monash University, VPAC. (david.may@sci.maths.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**
**
** Role:
**	Defines any header info, such as new structures, needed by this plugin
**
** Assumptions:
**
** Comments:
**
** $Id $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Underworld_plugins_HRS_Erosion_Context_h__
#define __Underworld_plugins_HRS_Erosion_Context_h__

	struct Underworld_HRS_Erosion_Context {
		Mesh*			mesh;
		double			K;
		double			DT;
		double			first_t_erosion;
		double			vT;
                FeVariable*             v;
		AbstractContext*        ctx;
                TimeIntegrator*         timeIntegrator;
	};

#endif
