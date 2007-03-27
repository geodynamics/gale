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
**
** Assumptions:
**
** Invariants:
**
** Comments:
**
** $Id: LevelSet.h 2225 1970-01-02 13:48:23Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_plugins_LevelSet_LevelSet_h__
#define __StgFEM_plugins_LevelSet_LevelSet_h__
	

	/* Textual name of this class */
	extern const Type LevelSet_Type;

	/* Virtual function types */
	typedef Bool (VelocityFunc)( double* pos, double* vel, void* ctx );
	
	/* Class contents */
	#define __LevelSet \
		/* General info */ \
		__FieldVariable \
		\
		/* Virtual info */ \
		\
		/* LevelSet info ... */ \
		char*				velVarName; \
		Dictionary*			initDistDict; \
		\
		FieldVariable*			velVar; \
		unsigned			res; \
		unsigned			fac; \
		GRM				grm; \
		double*				min; \
		double*				max; \
		double*				step; \
		Sync*				nodeSync; \
		double				cutoff; \
		\
		unsigned			nEls; \
		unsigned			nElNodes; \
		unsigned**			elNodes; \
		unsigned*			baseNodeEl; \
		unsigned			nActEls; \
		unsigned*			actEls; \
		unsigned			nInactNodes; \
		unsigned*			inactNodes; \
		\
		double*				dists; \
		double*				vels; \
		double				minDist; \
		double				maxDist; \
		\
		/* FE info. */ \
		Quadrature*			qdr; \
		double				jacDet; \
		Function*			shapeFuncs; \
		double*				psiTbl; \
		double**			psiSqTbl; \
		double****			gradTbl; \
		double**			psiQdr; \
		double***			gradQdr; \
		\
		/* Solver info. */ \
		Matrix*				aMat; \
		Vector*				rhsVec; \
		Vector*				solVec; \
		MatrixSolver*			solver; \

	struct LevelSet { __LevelSet };


	/*-----------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	LevelSet* LevelSet_DefaultNew( Name name );

	/* Create a LevelSet */
	LevelSet* LevelSet_New( Name			name, 
				FieldVariable_Register*	fieldVariable_Register );

	/* Creation implementation */
	LevelSet* _LevelSet_New( SizeT						_sizeOfSelf, 
				 Type						type,
				 Stg_Class_DeleteFunction*			_delete,
				 Stg_Class_PrintFunction*			_print, 
				 Stg_Class_CopyFunction*			_copy, 
				 Stg_Component_DefaultConstructorFunction*	_defaultConstructor, 
				 Stg_Component_ConstructFunction*		_construct,
				 Stg_Component_BuildFunction*			_build,
				 Stg_Component_InitialiseFunction*		_initialise,
				 Stg_Component_ExecuteFunction*			_execute,
				 Stg_Component_DestroyFunction*			_destroy,
				 Name						name,
				 Bool						initFlag,
				 FieldVariable_InterpolateValueAtFunction*	_interpolateValueAt,
				 FieldVariable_GetValueFunction*		_getMinGlobalFieldMagnitude,
				 FieldVariable_GetValueFunction*		_getMaxGlobalFieldMagnitude,		
				 FieldVariable_GetCoordFunction*		_getMinAndMaxLocalCoords,
				 FieldVariable_GetCoordFunction*		_getMinAndMaxGlobalCoords,
				 Index						fieldComponentCount,
				 Dimension_Index				dim,
				 Bool                                     isCheckpointedAndReloaded,
				 MPI_Comm					communicator,
				 FieldVariable_Register*			fieldVariable_Register );


	/* Initialise a LevelSet */
	void LevelSet_Init( LevelSet*			self, 
			    FieldVariable_Register*	fV_Register );

	/* Initialisation implementation functions */
	void _LevelSet_Init( LevelSet* self );


	/*-----------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	/* Stg_Class_Delete implementation */
	void _LevelSet_Delete( void* levelSet );

	/* Print implementation */
	void _LevelSet_Print( void* levelSet, Stream* stream );

void _LevelSet_Construct( void* levelSet, Stg_ComponentFactory* cf, void* data ) ;
	void _LevelSet_Build( void* levelSet, void* data ) ;
	void _LevelSet_Execute( void* levelSet, void* data ) ;
	void _LevelSet_Destroy( void* levelSet, void* data ) ;
	void _LevelSet_Initialise( void* levelSet, void* data ) ;

	InterpolationResult _LevelSet_InterpValAt( void* levelSet, Coord gCrd, double* vecRes );

	double _LevelSet_GetMinField( void* levelSet );

	double _LevelSet_GetMaxField( void* levelSet );

	void _LevelSet_GetLocalCrds( void* levelSet, Coord min, Coord max );

	void _LevelSet_GetGlobalCrds( void* levelSet, Coord min, Coord max );


	/*-----------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void LevelSet_SetMesh( void* levelSet, 
			       GRM* grm, double* min, double* max, 
			       Sync* nodeSync );

	void LevelSet_ConvMesh( void* levelSet, 
				GRM* grm, Mesh* mesh );

	void LevelSet_SetVelocityField( void* levelSet, 
					FieldVariable* velVar );

	void LevelSet_SetVelocityFunc( void* levelSet, 
				       VelocityFunc* velFunc, void* ctx );

	void LevelSet_InitDistances( void* levelSet, 
				     double (func)( const unsigned* inds, const double* pos, void* ctx ), void* funcCtx );

	void LevelSet_Advance( void* levelSet, double dt );


	/*-----------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	void _LevelSet_Free( LevelSet* self );

	void _LevelSet_Solve( void* levelSet );

	void _LevelSet_BuildShapeFuncs( LevelSet* self, unsigned short nDims );

	void _LevelSet_Redistance( LevelSet* self, unsigned maxIts );

	void _LevelSet_InsertInactive( LevelSet* self );

	void _LevelSet_InitActiveEls( LevelSet* self );

	Bool _LevelSet_ValidateElNodes( GRM* grm, Sync* nodeSync, 
					unsigned* base, unsigned short curDim, 
					unsigned* elNodes, unsigned* curElNode );

	void _LevelSet_CalcMinMax( LevelSet* self );


#endif /* __StgFEM_plugins_LevelSet_LevelSet_h__ */
