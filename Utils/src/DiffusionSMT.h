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


#ifndef __PICellerator_Utils_DiffusionSMT_h__
#define __PICellerator_Utils_DiffusionSMT_h__

/** Textual name of this class */
extern const Type DiffusionSMT_Type;

typedef struct {
    double diffusion;
} DiffusionSMT_MaterialExt;

/** DiffusionSMT class contents */
#define __DiffusionSMT							\
    /* General info */							\
    __StiffnessMatrixTerm						\
									\
    /* Virtual info */							\
									\
    /* DiffusionSMT info */						\
									\
    Materials_Register*                                 materials_Register; \
    ExtensionInfo_Index                                 materialExtHandle; \
    MaterialSwarmVariable**                             diffusionSwarmVariables; \
    Index                                               materialSwarmCount;

struct DiffusionSMT { __DiffusionSMT };

DiffusionSMT* DiffusionSMT_New( 
    Name                                                name,
    StiffnessMatrix*                                    stiffnessMatrix,
    Swarm*                                              integrationSwarm );

DiffusionSMT* _DiffusionSMT_New( 
    SizeT                                               sizeOfSelf,  
    Type                                                type,
    Stg_Class_DeleteFunction*                           _delete,
    Stg_Class_PrintFunction*                            _print,
    Stg_Class_CopyFunction*                             _copy, 
    Stg_Component_DefaultConstructorFunction*           _defaultConstructor,
    Stg_Component_ConstructFunction*                    _construct,
    Stg_Component_BuildFunction*                        _build,
    Stg_Component_InitialiseFunction*                   _initialise,
    Stg_Component_ExecuteFunction*                      _execute,
    Stg_Component_DestroyFunction*                      _destroy,
    StiffnessMatrixTerm_AssembleElementFunction*        _assembleElement,		
    Name                                                name );
	
void DiffusionSMT_InitAll( 
    void*                                               matrixTerm,
    StiffnessMatrix*                                    stiffnessMatrix,
    Swarm*                                              integrationSwarm );

void _DiffusionSMT_Delete( void* matrixTerm );
void _DiffusionSMT_Print( void* matrixTerm, Stream* stream );

void* _DiffusionSMT_DefaultNew( Name name ) ;
void _DiffusionSMT_Construct( void* matrixTerm, Stg_ComponentFactory* cf, void* data ) ;
void _DiffusionSMT_Build( void* matrixTerm, void* data ) ;
void _DiffusionSMT_Initialise( void* matrixTerm, void* data ) ;
void _DiffusionSMT_Execute( void* matrixTerm, void* data ) ;
void _DiffusionSMT_Destroy( void* matrixTerm, void* data ) ;

void _DiffusionSMT_AssembleElement( 
    void*                                              matrixTerm,
    StiffnessMatrix*                                   stiffnessMatrix, 
    Element_LocalIndex                                 lElement_I, 
    SystemLinearEquations*                             sle,
    FiniteElementContext*                              context,
    double**                                           elStiffMat ) ;

#endif
