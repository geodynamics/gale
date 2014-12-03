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
** $Id: Init.c 1209 2008-08-25 01:16:22Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/Discretisation/Discretisation.h>
#include <StgFEM/SLE/SLE.h>

#include "Assembly.h"

#include <stdio.h>

Stream* StgFEM_Assembly_Debug = NULL;

/** Initialises the Linear Algebra package, then any init for this package
such as streams etc */
Bool StgFEM_Assembly_Init( int* argc, char** argv[] ) {
	Stg_ComponentRegister* componentRegister = Stg_ComponentRegister_Get_ComponentRegister();
	int tmp;
	
	Journal_Printf( Journal_Register( DebugStream_Type, (Name)"Context"  ), "In: %s\n", __func__ ); /* DO NOT CHANGE OR REMOVE */
	tmp = Stream_GetPrintingRank( Journal_Register( InfoStream_Type, (Name)"Context" )  );
	Stream_SetPrintingRank( Journal_Register( InfoStream_Type, (Name)"Context"  ), 0 );
	Journal_Printf( /* DO NOT CHANGE OR REMOVE */
		Journal_Register( InfoStream_Type, (Name)"Context"  ), 
		"StGermain FEM Assembly Library revision %s. Copyright (C) 2003-2005 VPAC.\n", VERSION );
	Stream_Flush( Journal_Register( InfoStream_Type, (Name)"Context" )  );
	Stream_SetPrintingRank( Journal_Register( InfoStream_Type, (Name)"Context"  ), tmp );
	
	/* initialise this level's streams */
	StgFEM_Assembly_Debug = Stream_RegisterChild( StgFEM_Debug, "Assembly" );
	
	Stg_ComponentRegister_Add( componentRegister, ThermalBuoyancyForceTerm_Type, (Name)"0", _ThermalBuoyancyForceTerm_DefaultNew  );
	Stg_ComponentRegister_Add( componentRegister, GradientStiffnessMatrixTerm_Type, (Name)"0", _GradientStiffnessMatrixTerm_DefaultNew  );
	Stg_ComponentRegister_Add( componentRegister, DivergenceMatrixTerm_Type, (Name)"0", _DivergenceMatrixTerm_DefaultNew  );
	Stg_ComponentRegister_Add( componentRegister, LaplacianStiffnessMatrixTerm_Type, (Name)"0", _LaplacianStiffnessMatrixTerm_DefaultNew  );
	Stg_ComponentRegister_Add( componentRegister, IsoviscousStressTensorTerm_Type, (Name)"0", _IsoviscousStressTensorTerm_DefaultNew  );
	Stg_ComponentRegister_Add( componentRegister, PressureGradMatrixTerm_Type, (Name)"0", _PressureGradMatrixTerm_DefaultNew  );
	Stg_ComponentRegister_Add( componentRegister, PressureGradForceTerm_Type, (Name)"0", _PressureGradForceTerm_DefaultNew  );
	Stg_ComponentRegister_Add( componentRegister, MassMatrixTerm_Type, (Name)"0", _MassMatrixTerm_DefaultNew  );

	RegisterParent( ThermalBuoyancyForceTerm_Type,     ForceTerm_Type );
	RegisterParent( GradientStiffnessMatrixTerm_Type,  StiffnessMatrixTerm_Type );
	RegisterParent( DivergenceMatrixTerm_Type,  StiffnessMatrixTerm_Type );
	RegisterParent( LaplacianStiffnessMatrixTerm_Type, StiffnessMatrixTerm_Type );
	RegisterParent( IsoviscousStressTensorTerm_Type,   StiffnessMatrixTerm_Type );
	RegisterParent( PressureGradMatrixTerm_Type,   StiffnessMatrixTerm_Type );
	RegisterParent( PressureGradForceTerm_Type,     ForceTerm_Type );
	RegisterParent( MassMatrixTerm_Type,   StiffnessMatrixTerm_Type );

	return True;
}


