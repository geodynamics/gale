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


#ifndef __PICellerator_Utils_BuoyancyForceTermThermoChem_h__
#define __PICellerator_Utils_BuoyancyForceTermThermoChem_h__

	typedef double (BuoyancyForceTermThermoChem_CalcRaTFunction) (
				void* forceTerm, 
				Swarm* swarm, 
				Element_LocalIndex lElement_I, 
				void* particle );

	typedef double (BuoyancyForceTermThermoChem_CalcRaCFunction) (
				void* forceTerm, 
				Swarm* swarm, 
				Element_LocalIndex lElement_I, 
				void* particle );

	/** Textual name of this class */
	extern const Type BuoyancyForceTermThermoChem_Type;

	typedef struct {
		double density;
	} BuoyancyForceTermThermoChem_MaterialExt;

	/** BuoyancyForceTermThermoChem class contents */
	#define __BuoyancyForceTermThermoChem \
		/* General info */ \
		__ForceTerm \
		\
		/* Virtual info */ \
		BuoyancyForceTermThermoChem_CalcRaTFunction*              _calcRaT;                          \
		\
		BuoyancyForceTermThermoChem_CalcRaCFunction*              _calcRaC;                          \
		\
		/* BuoyancyForceTermThermoChem info */ \
		FeVariable*                                         temperatureField;                  \
		double                                              RaT;                           \
		double                                              RaC;                           \
		Bool                                                adjust;                            \
		Materials_Register*                                 materials_Register;                \
		ExtensionInfo_Index                                 materialExtHandle;                 \
		MaterialSwarmVariable**                             densitySwarmVariables;             \
		Index                                               materialSwarmCount;

	struct BuoyancyForceTermThermoChem { __BuoyancyForceTermThermoChem };

	BuoyancyForceTermThermoChem* BuoyancyForceTermThermoChem_New( 
		Name                                                name,
		ForceVector*                                        forceVector,
		Swarm*                                              integrationSwarm,
		FeVariable*                                         temperatureField,
		double                                              RaT,
		double                                              RaC,
		Bool                                                adjust,
		Materials_Register*                                 materials_Register );

	BuoyancyForceTermThermoChem* _BuoyancyForceTermThermoChem_New( 
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
		BuoyancyForceTermThermoChem_CalcRaTFunction*              _calcRaT,
		BuoyancyForceTermThermoChem_CalcRaCFunction*              _calcRaC,
		Name                                                name );
	
	void _BuoyancyForceTermThermoChem_Delete( void* forceTerm );
	void _BuoyancyForceTermThermoChem_Print( void* forceTerm, Stream* stream );

	void* _BuoyancyForceTermThermoChem_DefaultNew( Name name ) ;
void _BuoyancyForceTermThermoChem_Construct( void* forceTerm, Stg_ComponentFactory* cf, void* data ) ;
	void _BuoyancyForceTermThermoChem_Build( void* forceTerm, void* data ) ;
	void _BuoyancyForceTermThermoChem_Initialise( void* forceTerm, void* data ) ;
	void _BuoyancyForceTermThermoChem_Execute( void* forceTerm, void* data ) ;
	void _BuoyancyForceTermThermoChem_Destroy( void* forceTerm, void* data ) ;

	void _BuoyancyForceTermThermoChem_AssembleElement( void* forceTerm, ForceVector* forceVector, Element_LocalIndex lElement_I, double* elForceVec ) ;
	double _BuoyancyForceTermThermoChem_CalcRaT( void* forceTerm, Swarm* swarm, Element_LocalIndex lElement_I, void* particle ) ;
	#define BuoyancyForceTermThermoChem_CalcRaT( forceTerm, swarm, lElement_I, particle )\
		(( (BuoyancyForceTermThermoChem*) forceTerm )->_calcRaT( forceTerm, swarm, lElement_I, particle ) )

	double _BuoyancyForceTermThermoChem_CalcRaC( void* forceTerm, Swarm* swarm, Element_LocalIndex lElement_I, void* particle ) ;

	#define BuoyancyForceTermThermoChem_CalcRaC( forceTerm, swarm, lElement_I, particle )\
		(( (BuoyancyForceTermThermoChem*) forceTerm )->_calcRaC( forceTerm, swarm, lElement_I, particle ) )

#endif
