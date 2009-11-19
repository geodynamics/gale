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


#ifndef __PICellerator_Utils_HydrostaticTerm_h__
#define __PICellerator_Utils_HydrostaticTerm_h__

	/** Textual name of this class */
	extern const Type HydrostaticTerm_Type;

	/** HydrostaticTerm class contents */
	#define __HydrostaticTerm \
		/* General info */ \
		__Stg_Component \
		\
		/* Virtual info */ \
                double upper_density; \
                double upper_alpha; \
                double lower_density; \
                double lower_alpha; \
                double height; \
                double material_boundary; \
                double T_0; \
                double linear_coefficient; \
                double exponential_coefficient1; \
                double exponential_coefficient2; \
                double gravity; \
                double v; \
                double width; \
                AbstractContext *context; \



	struct HydrostaticTerm { __HydrostaticTerm };

	HydrostaticTerm* HydrostaticTerm_New( 
		Name                                                name,
                double upper_density,
                double upper_alpha,
                double lower_density,
                double lower_alpha,
                double height,
                double material_boundary,
                double T_0,
                double linear_coefficient,
                double exponential_coefficient1,
                double exponential_coefficient2,
                double gravity,
                double v,
                double width,
                AbstractContext *context);

	HydrostaticTerm* _HydrostaticTerm_New( 
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
                double upper_density,
                double upper_alpha,
                double lower_density,
                double lower_alpha,
                double height,
                double material_boundary,
                double T_0,
                double linear_coefficient,
                double exponential_coefficient1,
                double exponential_coefficient2,
                double gravity,
                double v,
                double width,
                AbstractContext *context,
		Name                                                name );
	
        void HydrostaticTerm_InitAll(void* forceTerm,
                                     double upper_density,
                                     double upper_alpha,
                                     double lower_density,
                                     double lower_alpha,
                                     double height,
                                     double material_boundary,
                                     double T_0,
                                     double linear_coefficient,
                                     double exponential_coefficient1,
                                     double exponential_coefficient2,
                                     double gravity,
                                     double v,
                                     double width,
                                     AbstractContext *context);

	void _HydrostaticTerm_Delete( void* forceTerm );
	void _HydrostaticTerm_Print( void* forceTerm, Stream* stream );

	void* _HydrostaticTerm_DefaultNew( Name name ) ;
        void _HydrostaticTerm_Construct( void* forceTerm,
                                         Stg_ComponentFactory* cf,
                                         void* data ) ;
	void _HydrostaticTerm_Build( void* forceTerm, void* data ) ;
	void _HydrostaticTerm_Initialise( void* forceTerm, void* data ) ;
	void _HydrostaticTerm_Execute( void* forceTerm, void* data ) ;
	void _HydrostaticTerm_Destroy( void* forceTerm, void* data ) ;
        double HydrostaticTerm_Density( void* forceTerm, Coord coord);
        double HydrostaticTerm_Pressure( void* forceTerm, Coord coord);

#endif
