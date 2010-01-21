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
** Copyright (C) 2008, California Institute of Technology
** Modified for DivergenceForce by Walter Landry
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


#ifndef __Gale_Utils_DivergenceForce_h__
#define __Gale_Utils_DivergenceForce_h__

	/** Textual name of this class */
	extern const Type DivergenceForce_Type;

	/** DivergenceForce class contents */
	#define __DivergenceForce \
		/* General info */ \
		__ForceTerm \
		Stg_Shape*                                 domainShape;    \
                FeMesh*                                    geometryMesh; \
                StressBC_Entry                             force; \
                FiniteElementContext*                      context; \

	struct DivergenceForce { __DivergenceForce };

	DivergenceForce* DivergenceForce_New( 
		Name                                                name,
		ForceVector*                                        forceVector,
                Swarm* integrationSwarm,
                Stg_Shape* domainShape,
                FeMesh* geometryMesh,
                StressBC_Entry force,
                FiniteElementContext* context);

	DivergenceForce* _DivergenceForce_New( 
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
		ForceTerm_AssembleElementFunction*                  _assembleElement,		
		Name                                                name,
                FiniteElementContext*                               context);
	
	void DivergenceForce_InitAll( 
		void*                                               forceTerm,
		ForceVector*                                        forceVector,
		Swarm*                                              integrationSwarm,
                Stg_Shape* domainShape,
                FeMesh* geometryMesh,
                StressBC_Entry force,
                FiniteElementContext* context);

	void _DivergenceForce_Delete( void* forceTerm );
	void _DivergenceForce_Print( void* forceTerm, Stream* stream );

	void* _DivergenceForce_DefaultNew( Name name ) ;
void _DivergenceForce_Construct( void* forceTerm, Stg_ComponentFactory* cf, void* data ) ;
	void _DivergenceForce_Build( void* forceTerm, void* data ) ;
	void _DivergenceForce_Initialise( void* forceTerm, void* data ) ;
	void _DivergenceForce_Execute( void* forceTerm, void* data ) ;
	void _DivergenceForce_Destroy( void* forceTerm, void* data ) ;

	void _DivergenceForce_AssembleElement( void* forceTerm, ForceVector* forceVector, Element_LocalIndex lElement_I, double* elForceVec ) ;

#endif
