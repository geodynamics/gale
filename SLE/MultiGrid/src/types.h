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
** $Id: types.h 200 2005-07-08 08:24:41Z AlanLo $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_SLE_MultiGrid_types_h__
#define __StgFEM_SLE_MultiGrid_types_h__
	
	
	/* 
	** Support stucts for MG.
	*/
	
	typedef struct {
		unsigned	nNodes;
		unsigned*	nodesTop;

		unsigned	nEls;
		unsigned*	nIncNodesLocal;
		unsigned**	incNodesLocal;
		unsigned*	nIncNodesPrev;
		unsigned**	incNodesPrev;
		unsigned*	incElTop;

		unsigned*	nIncElsLocal;
		unsigned**	incElsLocal;
		unsigned**	elNodeInds;

		unsigned	nOwnedEqNums;
		unsigned	nUCDofs;
		unsigned*	nDofs;
		unsigned**	lmTable;
	} MGMapping;
	
	
	typedef struct {
		FeVariable*	feVar;
		unsigned	maxLevels;
		MGMapping*	maps; /* x maxLevels */
	} MGGridMapping;
	
	
	typedef struct {
		StiffnessMatrix*	stiffMat;
		unsigned		nLevels;
		MGGridMapping*		mapping;
		
		unsigned*		nUpIts;
		unsigned*		nDownIts;
		unsigned*		nCycles;
		unsigned		nFinalIts;
		char**			downSmoothers;
		char**			upSmoothers;
		Matrix**		smoothers;
		Matrix**		rOps;
		Matrix**		pOps;
		Matrix**		iOps;
		
		unsigned*		nLocalRows; /* use this guy to break up coarse work vectors */
		unsigned*		nLocalCols; /* use this guy to break up coarse work vectors */
		Vector**		rhsVecs; /* work vector */
		Vector**		xVecs; /* work vector */
		Vector**		rVecs; /* work vector */
	} MGInfo;
	
	
	typedef struct {
		unsigned		nMappings;
		MGGridMapping*	mappings;
		unsigned		nInfos;
		MGInfo*		infos;
	} MGContext;
	
	
	/* 
	** MG classes.
	*/
	
	typedef struct MGCoarsener			MGCoarsener;
	typedef struct MGCoarsener_RegCartesian	MGCoarsener_RegCartesian;
	
	
#endif /* __StgFEM_SLE_MultiGrid_types_h__ */
