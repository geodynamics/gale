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
** $Id: /cig/src/PICellerator/Utils/src/HydrostaticTerm.c 2420 2008-12-12T20:12:30.951338Z boo  $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PopulationControl/PopulationControl.h>
#include <PICellerator/Weights/Weights.h>
#include <PICellerator/MaterialPoints/MaterialPoints.h>

#include "types.h"
#include "HydrostaticTerm.h"
#include "MaterialSwarmVariable.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class */
const Type HydrostaticTerm_Type = "HydrostaticTerm";

HydrostaticTerm* HydrostaticTerm_New(Name name,
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
                                     AbstractContext *context)
{
	HydrostaticTerm* self = (HydrostaticTerm*) _HydrostaticTerm_DefaultNew( name );

	HydrostaticTerm_InitAll(self,
                                upper_density,
                                upper_alpha,
                                lower_density,
                                lower_alpha,
                                height,
                                material_boundary,
                                T_0,
                                linear_coefficient,
                                exponential_coefficient1,
                                exponential_coefficient2,
                                gravity,
                                v,
                                width,
                                context);
	return self;
}

/* Creation implementation / Virtual constructor */
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
		Name                                                name )
{
  HydrostaticTerm* self;
  
  /* Allocate memory */
  assert( sizeOfSelf >= sizeof(HydrostaticTerm) );
  self = (HydrostaticTerm*)_Stg_Component_New(sizeOfSelf,
                                              type, 
                                              _delete,
                                              _print,
                                              _copy,
                                              _defaultConstructor,
                                              _construct,
                                              _build,
                                              _initialise,
                                              _execute,
                                              _destroy,
                                              name,
                                              NON_GLOBAL );

  /* Virtual info */
  
  self->upper_density=upper_density;
  self->upper_alpha=upper_alpha;
  self->lower_density=lower_density;
  self->lower_alpha=lower_alpha;
  self->height=height;
  self->material_boundary=material_boundary;
  self->T_0=T_0;
  self->linear_coefficient=linear_coefficient;
  self->exponential_coefficient1=exponential_coefficient1;
  self->exponential_coefficient2=exponential_coefficient2;
  self->gravity=gravity;
  self->v=v;
  self->width=width;
  self->context=context;
  return self;
}

void _HydrostaticTerm_Init(HydrostaticTerm* self,
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
                           AbstractContext *context)
{
  self->isConstructed    = True;

  self->upper_density=upper_density;
  self->upper_alpha=upper_alpha;
  self->lower_density=lower_density;
  self->lower_alpha=lower_alpha;
  self->height=height;
  self->material_boundary=material_boundary;
  self->T_0=T_0;
  self->linear_coefficient=linear_coefficient;
  self->exponential_coefficient1=exponential_coefficient1;
  self->exponential_coefficient2=exponential_coefficient2;
  self->gravity=gravity;
  self->v=v;
  self->width=width;
  self->context=context;
}

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
                             AbstractContext *context)
{
  HydrostaticTerm* self = (HydrostaticTerm*) forceTerm;
  _HydrostaticTerm_Init( self,
                         upper_density,
                         upper_alpha,
                         lower_density,
                         lower_alpha,
                         height,
                         material_boundary,
                         T_0,
                         linear_coefficient,
                         exponential_coefficient1,
                         exponential_coefficient2,
                         gravity,
                         v,
                         width,
                         context);
}

void _HydrostaticTerm_Delete( void* forceTerm ) {
}

void _HydrostaticTerm_Print( void* forceTerm, Stream* stream ) {
}

void* _HydrostaticTerm_DefaultNew( Name name ) {
	return (void*)_HydrostaticTerm_New( 
		sizeof(HydrostaticTerm), 
		HydrostaticTerm_Type,
		_HydrostaticTerm_Delete,
		_HydrostaticTerm_Print,
		NULL,
		_HydrostaticTerm_DefaultNew,
		_HydrostaticTerm_Construct,
		_HydrostaticTerm_Build,
		_HydrostaticTerm_Initialise,
		_HydrostaticTerm_Execute,
		_HydrostaticTerm_Destroy,
                0,0,0,0,0,0,0,0,0,0,0,0,0,NULL,
		name );
}

void _HydrostaticTerm_Construct( void* forceTerm, Stg_ComponentFactory* cf,
                                 void* data ) {
  HydrostaticTerm* self = (HydrostaticTerm*)forceTerm;
  double upper_density,upper_alpha,lower_density,lower_alpha,height,
    material_boundary,T_0,linear_coefficient,exponential_coefficient1,
    exponential_coefficient2,gravity,v,width;
  AbstractContext *context;

  context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data );

  upper_density=
    Stg_ComponentFactory_GetDouble(cf,self->name,"upperDensity",0.0);
  upper_alpha=
    Stg_ComponentFactory_GetDouble(cf,self->name,"upperAlpha",0.0);
  lower_density=
    Stg_ComponentFactory_GetDouble(cf,self->name,"lowerDensity",0.0);
  lower_alpha=Stg_ComponentFactory_GetDouble(cf,self->name,"lowerAlpha",0.0);
  height=Stg_ComponentFactory_GetDouble(cf,self->name,"height",0.0);
  material_boundary=
    Stg_ComponentFactory_GetDouble(cf,self->name,"materialBoundary",0.0);
  T_0=Stg_ComponentFactory_GetDouble(cf,self->name,"T_0",0.0);
  linear_coefficient=
    Stg_ComponentFactory_GetDouble(cf,self->name,"linearCoefficient",0.0);
  exponential_coefficient1=
    Stg_ComponentFactory_GetDouble(cf,self->name,"exponentialCoefficient1",0.0);
  exponential_coefficient2=
    Stg_ComponentFactory_GetDouble(cf,self->name,"exponentialCoefficient2",0.0);
  gravity=Stg_ComponentFactory_GetDouble(cf,self->name,"gravity",0.0);

  v=Stg_ComponentFactory_GetDouble(cf,self->name,"expansionVelocity",0.0);
  width=Stg_ComponentFactory_GetDouble(cf,self->name,"width",0.0);

  _HydrostaticTerm_Init( self, upper_density,
                         upper_alpha,
                         lower_density,
                         lower_alpha,
                         height,
                         material_boundary,
                         T_0,
                         linear_coefficient,
                         exponential_coefficient1,
                         exponential_coefficient2,
                         gravity,
                         v,
                         width,
                         context);
}

void _HydrostaticTerm_Build( void* forceTerm, void* data ) {
}

void _HydrostaticTerm_Initialise( void* forceTerm, void* data ) {
}

void _HydrostaticTerm_Execute( void* forceTerm, void* data ) {
}

void _HydrostaticTerm_Destroy( void* forceTerm, void* data ) {
}

void HydrostaticTerm_Current_Heights(void* forceTerm, double *current_height,
                                     double *current_boundary)
{
  HydrostaticTerm *self=(HydrostaticTerm *)forceTerm;
  double current_thickness;
  double t=self->context->currentTime;

  if(self->v!=0)
    {
      current_thickness=(self->height - self->material_boundary)
        *self->width / (self->width + self->v * t);
      *current_boundary=self->material_boundary
        + (self->upper_density / self->lower_density)
        * (self->height - self->material_boundary - current_thickness);
    }
  else
    {
      current_thickness=self->height - self->material_boundary;
      *current_boundary=self->material_boundary;
    }
  *current_height=current_thickness + *current_boundary;
}

double HydrostaticTerm_Density( void* forceTerm, Coord coord)
{
  HydrostaticTerm *self=(HydrostaticTerm *)forceTerm;
  double h=self->height-coord[1];
  double T, density, alpha;
  double current_height, current_boundary;

  HydrostaticTerm_Current_Heights(self,&current_height,&current_boundary);
  /* printf("Time %g %g %g\n",current_height,current_boundary,self->context->currentTime); */

  /* First, get the temperature */
  if(h<0)
    {
      T=self->T_0;
    }
  else
    {
      T=self->T_0 + self->linear_coefficient*h
        + self->exponential_coefficient1
        *(1-exp(-self->exponential_coefficient2*h));
    }

  /* Then, get the density */
  if(coord[1]>current_height)
    {
      density=alpha=0;
    }
  else if(coord[1]>current_boundary)
    {
      density=self->upper_density;;
      alpha=self->upper_alpha;
    }
  else
    {
      density=self->lower_density;
      alpha=self->lower_alpha;
    }

  return density*(1-alpha*T);
}

double HydrostaticTerm_Pressure_Analytic(double density, double gravity,
                                         double h, double alpha,
                                         double T_0, double A,
                                         double B, double C)
{
  if(C==0)
    return density*gravity*h*(1-alpha*(T_0 + A*h/2));
  else
    return density*gravity*h*(1-alpha*(T_0 + A*h/2 + B*(1+(exp(-C*h)-1)/(C*h))));
}

double HydrostaticTerm_Pressure( void* forceTerm, Coord coord)
{
  HydrostaticTerm *self=(HydrostaticTerm *)forceTerm;
  double h=self->height-coord[1];
  double p;

  if(coord[1]>=self->height)
    {
      p=0;
    }
  else if(coord[1]>self->material_boundary)
    {
      p=HydrostaticTerm_Pressure_Analytic(self->upper_density, self->gravity, h,
                                          self->upper_alpha,
                                          self->T_0,
                                          self->linear_coefficient,
                                          self->exponential_coefficient1,
                                          self->exponential_coefficient2);
    }
  else
    {
      p=HydrostaticTerm_Pressure_Analytic(self->upper_density,
                                          self->gravity,
                                          self->height-self->material_boundary,
                                          self->upper_alpha,
                                          self->T_0,
                                          self->linear_coefficient,
                                          self->exponential_coefficient1,
                                          self->exponential_coefficient2)
        - HydrostaticTerm_Pressure_Analytic(self->lower_density, self->gravity,
                                            self->height-self->material_boundary,
                                            self->lower_alpha,
                                            self->T_0,
                                            self->linear_coefficient,
                                            self->exponential_coefficient1,
                                            self->exponential_coefficient2)
        + HydrostaticTerm_Pressure_Analytic(self->lower_density, self->gravity,
                                            h,
                                            self->lower_alpha,
                                            self->T_0,
                                            self->linear_coefficient,
                                            self->exponential_coefficient1,
                                            self->exponential_coefficient2);
    }
  return p;
}

