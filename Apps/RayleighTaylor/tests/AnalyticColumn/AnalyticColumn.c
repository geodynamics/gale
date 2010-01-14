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
** $Id: AnalyticColumn.c 518 2007-10-11 08:07:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#define MAX_FOURIER_TERMS 220

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include <assert.h>

const Type AnalyticColumn_Type = "AnalyticColumn";

typedef struct {
	__AnalyticSolution
	FeVariable*		velocityField;
	FeVariable*		pressureField;
	Dimension_Index         dim;
	double                  sigma;
	double                  viscosity;
	double                  startColumnX;
	double                  endColumnX;
	double                  startColumnZ;
	double                  endColumnZ;
} AnalyticColumn;

void AnalyticColumn_Constants( void* analyticSolution, double n, double* C1, double* C2, double* C3, double *C4 ) {
	AnalyticColumn*          self           = (AnalyticColumn*)analyticSolution;
	double                  viscosity      = self->viscosity;
	double                  deltaRho       = self->sigma;
	double                  x0             = 0.5 * ( self->startColumnX + self->endColumnX );
	double                  dx             = self->endColumnX - self->startColumnX;
	double                  factor = M_PI * n;

	*C1 = 0.2e1 * deltaRho * sin(factor * dx / 0.2e1) * cos(factor * x0) * exp(-factor) * (0.2e1 * exp(-factor) + 0.2e1 + factor) * pow(exp(-factor) + 0.1e1, -0.2e1) / viscosity * pow(M_PI, -0.3e1) * pow(n, -0.3e1); 
		
	*C2 = -0.2e1 * (-0.2e1 * exp(-factor) + factor * exp(-factor) - 0.2e1) * deltaRho * sin(factor * dx / 0.2e1) * cos(factor * x0) * pow(exp(-factor) + 0.1e1, -0.2e1) / viscosity * pow(M_PI, -0.3e1) * pow(n, -0.3e1); 
		
	*C3 = -0.2e1 * deltaRho * sin(factor * dx / 0.2e1) * cos(factor * x0) * exp(-factor) / viscosity * pow(n, -0.2e1) * pow(M_PI, -0.2e1) / (exp(-factor) + 0.1e1); 
		
	*C4 = 0.2e1 * deltaRho * sin(factor * dx / 0.2e1) * cos(factor * x0) / viscosity * pow(n, -0.2e1) * pow(M_PI, -0.2e1) / (exp(-factor) + 0.1e1);
}


void _AnalyticColumn_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity ) {
	AnalyticColumn*         self           = (AnalyticColumn*)analyticSolution;
	double                  n, m;
	double                  x, y, z;
	double                  C1, C2, C3, C4;
	double                  viscosity      = self->viscosity;
	double                  sigma          = self->sigma;
	double                  deltaRho;
	double                  x0;
	double                  z0;
	double                  dx;
	double                  dz;
	double                  u1, u2, u3;
	double                  factor;
	
	/* Find coordinate of node */
	x = coord[ I_AXIS ];
	y = coord[ J_AXIS ];
	z = coord[ K_AXIS ];

	x0 = 0.5 * ( self->startColumnX + self->endColumnX );
	dx = self->endColumnX - self->startColumnX;
	z0 = 0.5 * ( self->startColumnZ + self->endColumnZ );
	dz = self->endColumnZ - self->startColumnZ;

	/* Initialise */
	if ( self->dim == 2 ) {
		velocity[ I_AXIS ] = 0.0;
		velocity[ J_AXIS ] = 0.0;

		for( n = 1.0 ; n < MAX_FOURIER_TERMS ; n++ ) {
			factor = M_PI * n;

			AnalyticColumn_Constants( self, n, &C1, &C2, &C3, &C4 );

			u1 = (-0.2e1 * sigma * sin(factor * x0 + factor * dx / 0.2e1) + 0.2e1 * sigma * sin(factor * x0 - factor * dx / 0.2e1) + C1 * exp(factor * y) * pow(n, 0.3e1) * pow(M_PI, 0.3e1) * viscosity + C2 * exp(-factor * y) * pow(n, 0.3e1) * pow(M_PI, 0.3e1) * viscosity + C3 * exp(factor * y) * y * pow(n, 0.3e1) * pow(M_PI, 0.3e1) * viscosity + C4 * exp(-factor * y) * y * pow(n, 0.3e1) * pow(M_PI, 0.3e1) * viscosity) * pow(n, -0.3e1) * pow(M_PI, -0.3e1) / viscosity;

			u2 = (-C1 * factor * exp(factor * y) + C2 * factor * exp(-factor * y) - C3 * factor * exp(factor * y) * y - C3 * exp(factor * y) + C4 * factor * exp(-factor * y) * y - C4 * exp(-factor * y)) / n / M_PI;
			
			u2 *= sin( factor * x );
			velocity[ I_AXIS ] += u2;
			u1 *= cos( factor * x );
			velocity[ J_AXIS ] += u1;
		}
	}
	else {
		double L1,kn,km; 
		double Am,Ap,Bm,Bp,C,D,E;

		velocity[ I_AXIS ] = 0.0;
		velocity[ J_AXIS ] = 0.0;
		velocity[ K_AXIS ] = 0.0;

		for ( n = 0 ; n < 45 ; n++ ) {
			for ( m = 0 ; m < 45 ; m++ ) {

				if ( n != 0 && m!=0 ){
					deltaRho = 4.0*sigma*sin((double)n*M_PI*dx)*sin((double)m*M_PI*dz)/(double)n/M_PI/(double)m/M_PI; 
				} 
				else if ( n == 0 && m != 0 ) { 
					deltaRho = 2.0*sigma*dx*sin((double)m*M_PI*dz)/(double)m/M_PI; 
				}
				else if  ( n != 0 && m == 0 ) { 
					deltaRho = 2.0*sigma*dz*sin((double)n*M_PI*dx)/(double)n/M_PI; 
				} 
				else { 
					deltaRho = sigma*dx*dz; 
				} 
				kn = M_PI * (double)n; 
				km = M_PI * (double)m;
				
				L1 = M_PI * sqrt( (double)(n*n + m*m)); 
				
				Am = exp((y-2.0)*L1)-exp(-y*L1);
				Ap = exp((y-2.0)*L1)+exp(-y*L1); 
				Bm = exp((y-1.0)*L1)-exp(-(y+1.0)*L1); 
				Bp = exp((y-1.0)*L1)+exp(-(y+1.0)*L1);

				C = (exp(-y*L1)-1.0)*(exp((y-1.0)*L1)-1.0); 
				D = exp((y-1.0)*L1)-exp(-y*L1); 
				E = (1.0+exp(-L1)); 
				
				u1 = (n!=0 || m!=0) 
					?  -( y*Am+(y-1.0)*Bm )*deltaRho/( 2*viscosity*L1*E*E ) - C*deltaRho/(viscosity*L1*L1*E) 
					:  0.0; 
				if ( m != 0 ){ 
					u2 = ( y*Ap+(y-1.0)*Bp )*deltaRho*km/( 2*viscosity*L1*L1*E*E ) - D*deltaRho*km/(2*viscosity*L1*L1*L1*E);
				} 
				else { 
					u2 = 0.0; 
				} 
				
				if ( n != 0 ){ 
					u3 = ( y*Ap+(y-1.0)*Bp )*deltaRho*kn/( 2*viscosity*L1*L1*E*E ) - D*deltaRho*kn/(2*viscosity*L1*L1*L1*E);
				} 
				else { 
					u3 = 0.0; 
				} 

				u1 *= cos(n*M_PI*x)*cos(m*M_PI*z); 
				velocity[ J_AXIS ] += u1;
				u2 *= cos(n*M_PI*x)*sin(m*M_PI*z);
				velocity[ K_AXIS ] += u2;
				u3 *= sin(n*M_PI*x)*cos(m*M_PI*z);
				velocity[ I_AXIS ] += u3;
			}
		}
	}
}

void _AnalyticColumn_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure ) {
	AnalyticColumn*          self           = (AnalyticColumn*)analyticSolution;
	double                  viscosity      = self->viscosity;
	double                  sigma          = self->sigma;
	double                  n, m;
	double                  x, y, z;
	double                  C1, C2, C3, C4;
	double                  deltaRho;
	double                  x0, z0;
	double                  dx, dz;
	double                  u2, u3, u4;
	double                  pp;
	double                  factor;
	
	/* Find coordinate of node */
	x = coord[ I_AXIS ];
	y = coord[ J_AXIS ];
	z = coord[ K_AXIS ];

	x0 = 0.5 * ( self->startColumnX + self->endColumnX );
	dx = self->endColumnX - self->startColumnX;
	z0 = 0.5 * ( self->startColumnZ + self->endColumnZ );
	dz = self->endColumnZ - self->startColumnZ;

	/* Initialise */
	*pressure = 0.0;

	if ( self->dim == 2 ) {
		for( n = 1.0 ; n < MAX_FOURIER_TERMS ; n++ ) {
			factor = M_PI * n;

			AnalyticColumn_Constants( self, n, &C1, &C2, &C3, &C4 );

			u2 = (-C1 * factor * exp(factor * y) + C2 * factor * exp(-factor * y) - C3 * factor * exp(factor * y) * y - C3 * exp(factor * y) + C4 * factor * exp(-factor * y) * y - C4 * exp(-factor * y)) / n / M_PI;
			
			u3 = -viscosity * (-C1 * exp(factor * y) + C2 * exp(-factor * y) - C3 * exp(factor * y) * y + C4 * exp(-factor * y) * y);
			*pressure += (double)( -2.0* viscosity * factor * u2 - u3*2.0*factor )*cos( factor * x );
		}
	}
	else {
		double L1,kn,km; 
		double Am,Ap,Bm,Bp,C,D,E;

		for ( n = 0 ; n < 45 ; n++ ) {
			for ( m = 0 ; m < 45 ; m++ ) {

				if ( n!=0 && m!=0 ){
					deltaRho = 4.0*sigma*sin((double)n*M_PI*dx)*sin((double)m*M_PI*dz)/(double)n/M_PI/(double)m/M_PI; 
				} 
				else if ( n==0 && m !=0) { 
					deltaRho = 2.0*sigma*dx*sin((double)m*M_PI*dz)/(double)m/M_PI; 
				}
				else if  ( n!=0 && m ==0) { 
					deltaRho = 2.0*sigma*dz*sin((double)n*M_PI*dx)/(double)n/M_PI; 
				} 
				else { 
					deltaRho = sigma*dx*dz; 
				} 
				kn = M_PI * (double)n; 
				km = M_PI * (double)m;
				
				L1 = M_PI * sqrt( (double)(n*n + m*m)); 
				
				Am = exp((y-2.0)*L1)-exp(-y*L1);
				Ap = exp((y-2.0)*L1)+exp(-y*L1); 
				Bm = exp((y-1.0)*L1)-exp(-(y+1.0)*L1); 
				Bp = exp((y-1.0)*L1)+exp(-(y+1.0)*L1);

				C = (exp(-y*L1)-1.0)*(exp((y-1.0)*L1)-1.0); 
				D = exp((y-1.0)*L1)-exp(-y*L1); 
				E = (1.0+exp(-L1)); 
				
				if ( m != 0 ){ 
					u2 = ( y*Ap+(y-1.0)*Bp )*deltaRho*km/( 2*viscosity*L1*L1*E*E ) - D*deltaRho*km/(2*viscosity*L1*L1*L1*E);
				} 
				else { 
					u2 = 0.0; 
				} 
				
				if ( n != 0 ){ 
					u3 = ( y*Ap+(y-1.0)*Bp )*deltaRho*kn/( 2*viscosity*L1*L1*E*E ) - D*deltaRho*kn/(2*viscosity*L1*L1*L1*E);
				} 
				else { 
					u3 = 0.0; 
				} 
				
				u4 = (n == 0 && m == 0 ) 
					?  deltaRho*(y-0.5) 
					:  -( y*Ap+(y-1.0)*Bp )*deltaRho/( E*E ) + 2.0*D*deltaRho/(L1*E) ; 
				
				pp = -u4-2.0*viscosity*(kn*u3+km*u2); 
				pp *= cos(n*M_PI*x)*cos(m*M_PI*z); 
				*pressure += pp; /* total pressure */
			}
		}
	}
}

	
void _AnalyticColumn_AssignFromXML( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	AnalyticColumn*          self           = (AnalyticColumn*)analyticSolution;

	/* Construct Parent */
	_AnalyticSolution_AssignFromXML( self, cf, data );

	/* Create Analytic Fields */
	self->velocityField = Stg_ComponentFactory_ConstructByName( cf, (Name)"VelocityField", FeVariable, True, data  ); 
	self->pressureField = Stg_ComponentFactory_ConstructByName( cf, (Name)"PressureField", FeVariable, True, data  ); 

	self->dim          = Stg_ComponentFactory_GetRootDictUnsignedInt( cf, (Dictionary_Entry_Key)"dim", 0  );
	self->startColumnX = Stg_ComponentFactory_GetRootDictDouble( cf, (Dictionary_Entry_Key)"startColumnX", 0.0  );
	self->endColumnX   = Stg_ComponentFactory_GetRootDictDouble( cf, (Dictionary_Entry_Key)"endColumnX", 0.0  );
	self->startColumnZ = Stg_ComponentFactory_GetRootDictDouble( cf, (Dictionary_Entry_Key)"startColumnZ", 0.0  );
	self->endColumnZ   = Stg_ComponentFactory_GetRootDictDouble( cf, (Dictionary_Entry_Key)"endColumnZ", 0.0  );
	self->viscosity    = Stg_ComponentFactory_GetRootDictDouble( cf, (Dictionary_Entry_Key)"viscosity", 1.0  );
	self->sigma        = Stg_ComponentFactory_GetRootDictDouble( cf, (Dictionary_Entry_Key)"sigma", 1.0  );
}

void _AnalyticColumn_Build( void* analyticSolution, void* data ) {
	AnalyticColumn*	self = (AnalyticColumn*)analyticSolution;

	assert( self && Stg_CheckType( self, AnalyticColumn ) );

	Build( self->velocityField, data, False );
	Build( self->pressureField, data, False );
	AnalyticSolution_CreateAnalyticField( self, self->velocityField, _AnalyticColumn_VelocityFunction );
	AnalyticSolution_CreateAnalyticField( self, self->pressureField, _AnalyticColumn_PressureFunction );

	_AnalyticSolution_Build( self, data );
}

void* _AnalyticColumn_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(AnalyticColumn);
	Type                                                      type = AnalyticColumn_Type;
	Stg_Class_DeleteFunction*                              _delete = _AnalyticSolution_Delete;
	Stg_Class_PrintFunction*                                _print = _AnalyticSolution_Print;
	Stg_Class_CopyFunction*                                  _copy = _AnalyticSolution_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _AnalyticColumn_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _AnalyticColumn_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _AnalyticColumn_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _AnalyticSolution_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _AnalyticSolution_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _AnalyticSolution_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return _AnalyticSolution_New(  ANALYTICSOLUTION_PASSARGS  );
}

Index _PICellerator_AnalyticColumn_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, AnalyticColumn_Type, (Name)"0", _AnalyticColumn_DefaultNew  );
}


