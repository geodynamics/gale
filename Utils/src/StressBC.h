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
**  Modified 2006 by Walter Landry (California Institute of Technology)
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


#ifndef __Gale_Utils_StressBC_h__
#define __Gale_Utils_StressBC_h__

	/** Textual name of this class */
	extern const Type StressBC_Type;

	/** StressBC class contents */
	#define __StressBC \
		/* General info */ \
		__ForceTerm \
		\
		/* StressBC info */ \
		Wall			         _wall; \
                StressBC_Entry                           _entryTbl[3]; \
                int                                      numEntries; \
		ConditionFunction_Register*		 conFunc_Register; \
                double                                   bottom_density; \

	struct StressBC { __StressBC };

	#ifndef ZERO
	#define ZERO 0
	#endif

	#define STRESSBC_DEFARGS \
                FORCETERM_DEFARGS

	#define STRESSBC_PASSARGS \
                FORCETERM_PASSARGS

        StressBC* _StressBC_New( STRESSBC_DEFARGS);

	StressBC* StressBC_New( 
		Name                                                name,
                FiniteElementContext*	context,
		ForceVector*                                        forceVector,
		Swarm*                                              integrationSwarm,
                ConditionFunction_Register*			    conFunc_Register);
	void StressBC_InitAll( 
		void*                                               forceTerm,
                FiniteElementContext*	context,
		ForceVector*                                        forceVector,
		Swarm*                                              integrationSwarm,
                ConditionFunction_Register*			    conFunc_Register);

	void _StressBC_Delete( void* forceTerm );
	void _StressBC_Print( void* forceTerm, Stream* stream );

	void* _StressBC_DefaultNew( Name name ) ;
        void _StressBC_AssignFromXML( void* forceTerm, Stg_ComponentFactory* cf, void* data ) ;
	void _StressBC_Build( void* forceTerm, void* data ) ;
	void _StressBC_Initialise( void* forceTerm, void* data ) ;
	void _StressBC_Execute( void* forceTerm, void* data ) ;
	void _StressBC_Destroy( void* forceTerm, void* data ) ;

	void _StressBC_AssembleElement( void* forceTerm, ForceVector* forceVector, Element_LocalIndex lElement_I, double* elForceVec ) ;
        unsigned StressBC_get_overcount(Dimension_Index dim, IJK ijk,
                                        unsigned sizes[]);

#endif
