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

#ifndef __Underworld_RecoveredFeVariable_h__
#define __Underworld_RecoveredFeVariable_h__

	/** Textual name of this class */
	extern const Type RecoveredFeVariable_Type;

	typedef void (REPPutElIntoProc)( RecoveredFeVariable* self, int lElement_I, double*** elHi_Mat, double** elFi_Mat );

	typedef void (REPCalcHi)( RecoveredFeVariable* self, double** GNx, double** pMatrix, double factor, double** Hi );

	typedef void (REPCalcFi)( RecoveredFeVariable* self, double** GNx, double* p_Stress, double factor, double* Fi );

	extern Stg_ObjectList* repRequiredRawFields_Reg;
	/* this functionPtr could be useful later, JG, 12May09
	typedef void (REPAssembleOnParticle)( RecoveredFeVariable*, ConstitutiveMatrix*, int, IntegrationPoint*, int, double*, double**, double, double*** elHi_Mat, double** elFi_Mat );
		REPAssembleOnParticle*        _assembleOnParticle; \
		*/

	/** RecoveredFeVariable class contents */ 
	#define __RecoveredFeVariable \
		/* parent class */ \
    __BaseRecoveryFeVar \
		/* RecoveredFeVariable info */ \
		/* virtual functions */ \
		REPPutElIntoProc*			_putElIntoProc; \
		REPCalcHi*					_calcHi; \
		REPCalcFi*					_calcFi; \
		/* class info */ \
		FeVariable*					rawPressureField; \
		FeVariable*					velGradField; \
		Bool							refreshMeshConnectivity; \
		Bool							recoverStrain; \
		Bool							nonLinearProblem; \
		double*						elementRep_H; /* [domainElNum][dofThatExist][nodesPerEl*dim][orderOfInterpolation] */ \
		double*						elementRep_F; /* [domainElNum][dofThatExist][nodesPerEl*dim] */ \
		int							nodesInPatch; \
		int							nodesPerEl; \
		int							coeffCount; \
		FILE*							dStream; \
		/* working space */ \
		double***					patch_H; \
		double**						patch_F; \
		double**						H_t; \
		double***					AMat; \
		double**						bVec; \
		double**						pMatrix;\
		double**						tmpEl_H; \
		double*						tmpEl_F; \
		double**						CPmat; \
		double**						tmpC; \
		
	struct RecoveredFeVariable { __RecoveredFeVariable };
	
	/* --- Contstructors / Destructors --- */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define RECOVEREDFEVARIABLE_DEFARGS \
                BASERECOVERYFEVAR_DEFARGS

	#define RECOVEREDFEVARIABLE_PASSARGS \
                BASERECOVERYFEVAR_PASSARGS

	RecoveredFeVariable* _RecoveredFeVariable_New(  RECOVEREDFEVARIABLE_DEFARGS  );
	
	void* _RecoveredFeVariable_DefaultNew( Name name );

   void _RecoveredFeVariable_Init( RecoveredFeVariable* self, StiffnessMatrix* stiffnessMatrix, FeVariable* rawPressureField, FeVariable* velGradField, Bool refreshMeshConnectivity, Bool recoverStrain );

	void _RecoveredFeVariable_Delete( void* recFeVariable );

	void _RecoveredFeVariable_Print( void* recFeVariable, Stream* stream );

	void* _RecoveredFeVariable_Copy( const void* recFeVariable, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

	void _RecoveredFeVariable_AssignFromXML( void* recFeVariable, Stg_ComponentFactory* cf, void* data );

	void _RecoveredFeVariable_Build( void* recFeVariable, void* data );

	void _RecoveredFeVariable_Initialise( void* recFeVariable, void* data );

	void _RecoveredFeVariable_Execute( void* recFeVariable, void* data );

	void _RecoveredFeVariable_Destroy( void* recFeVariable, void* data );

	/** Returns in particleValue the value represented on this particle */
	#define RecoveredFeVariable_ValueAtParticle( particleFeVariable, swarm, lElement_I, particle, particleValue ) \
		( (RecoveredFeVariable*)(particleFeVariable) )->_valueAtParticle( \
			(particleFeVariable), \
			(swarm), \
			(lElement_I), \
			(particle), \
			(particleValue) )

	void RecoveredFeVariable_Update( void* recFeVariable );

	void RecoveredFeVariable_CommunicateHF( RecoveredFeVariable* self );

	void RecoveredFeVariable_SolvePatch( RecoveredFeVariable* self, int node_I, int *elList, int nEl, LmStruct* lmStruct );

	void RecoveredFeVariable_SetupWorkSpace( RecoveredFeVariable* self );

	void RecoveredFeVariable_RemoveWorkSpace( RecoveredFeVariable* self );

	void _RecoveredFeVariable_AssembleAtParticle( RecoveredFeVariable* self, ConstitutiveMatrix* constitutiveMatrix,
		IntegrationPoint* particle, int lElement_I, 
		double* globalCoord, double** GNx, double detJac,
		double*** Hi_Mat, double** Fi_Mat	);

	void _RecoveredFeVariable_PutElIntoProc( RecoveredFeVariable* self, int lElement_I, double*** elHi_Mat, double** elFi_Mat );

	void _RecoveredFeVariable_UnpackElement( double* ptrH, double* ptrF, int rowsInH, int order, double **el_H, double *el_F );

	void _RecoveredFeVariable_ZeroHiFi( RecoveredFeVariable* self, double*** elHi_Mat, double** elFi_Mat );

	void _RecoveredFeVariable_CalcHi2D( RecoveredFeVariable* self, double** GNx, double** pMatrix, double factor, double** Hi );

	void _RecoveredFeVariable_CalcHi3D( RecoveredFeVariable* self, double** GNx, double** pMatrix, double factor, double** Hi );

	void _RecoveredFeVariable_CalcFi2D( RecoveredFeVariable* self, double** GNx, double* p_Stress, double factor, double* Fi );

	void _RecoveredFeVariable_CalcFi3D( RecoveredFeVariable* self, double** GNx, double* p_Stress, double factor, double* Fi );

	void _RecoveredFeVariable_Make_pVec( double* global, double* pVec );

	double _RecoveredFeVariable_ApplyCoeff( RecoveredFeVariable* self, double* coeff, double* coord, int order, double* pVec );

	void _RecoveredFeVariable_GetValueAtNode( void* feVariable, Node_DomainIndex dNode_I, double* value );

	void _RecoveredFeVariable_GetCoeffAtNode( void* feVariable, Node_DomainIndex dNode_I, double* coeff );

	void _RecoveredFeVariable_GetValueInElement( void* feVariable, Element_Index lEl_I, double* xi, double* value );

	void _TmpConstMat2D_CalculateStress( double** D, SymmetricTensor strainRate, SymmetricTensor stress, Bool isDiagonal );

	void _TmpConstMat3D_CalculateStress( double** D, SymmetricTensor strainRate, SymmetricTensor stress, Bool isDiagonal );

#endif 

