/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
** $Id: Simplex.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include "StGermain/StGermain.h"
#include "Geometry.h"


Bool Simplex_Search3D( double** verts, unsigned* inc, 
		       unsigned nSimplices, unsigned** inds, 
		       double* point, double* bc, unsigned* inside )
{
	unsigned	s_i;

	assert( inds );

	for( s_i = 0; s_i < nSimplices; s_i++ ) {
		unsigned	ind_i;

		Simplex_Barycenter3D( verts, inc, inds[s_i], point, bc );
		for( ind_i = 0; ind_i < 4; ind_i++ ) {
			if( bc[ind_i] < 0.0 || bc[ind_i] > 1.0 )
				break;
		}
		if( ind_i == 4 ) {
			*inside = s_i;
			return True;
		}
	}

	return False;
}

Bool Simplex_Search2D( double** verts, unsigned* inc, 
		       unsigned nSimplices, unsigned** inds, 
		       double* point, double* bc, unsigned* inside )
{
	unsigned	s_i;

	assert( inds );

	for( s_i = 0; s_i < nSimplices; s_i++ ) {
		unsigned	ind_i;

		Simplex_Barycenter2D( verts, inc, inds[s_i], point, bc );
		for( ind_i = 0; ind_i < 3; ind_i++ ) {
			if( bc[ind_i] < 0.0 || bc[ind_i] > 1.0 )
				break;
		}
		if( ind_i == 3 ) {
			*inside = s_i;
			return True;
		}
	}

	return False;
}

void Simplex_Barycenter3D( double** verts, unsigned* inc, unsigned* inds, double* point, double* bc ) {
	double*	tet[4];
	double	x0, y0, z0;
	double	x1, y1, z1;
	double	x2, y2, z2;
	double	x3, y3, z3;
	double	px, py, pz;
	double	den;

	assert( verts );
	assert( inc );
	assert( inds );
	assert( point );
	assert( bc );
	assert( verts[inc[inds[0]]] );
	assert( verts[inc[inds[1]]] );
	assert( verts[inc[inds[2]]] );
	assert( verts[inc[inds[3]]] );

	tet[0] = verts[inc[inds[0]]];
	tet[1] = verts[inc[inds[1]]];
	tet[2] = verts[inc[inds[2]]];
	tet[3] = verts[inc[inds[3]]];
	x0 = tet[0][0]; x1 = tet[1][0]; x2 = tet[2][0]; x3 = tet[3][0];
	y0 = tet[0][1]; y1 = tet[1][1]; y2 = tet[2][1]; y3 = tet[3][1];
	z0 = tet[0][2]; z1 = tet[1][2]; z2 = tet[2][2]; z3 = tet[3][2];
	px = point[0]; py = point[1]; pz = point[2];
	den = 1.0 / (x1*(y0*(z3 - z2) + y2*(z0 - z3) + y3*(z2 - z0)) + 
			     x0*(y2*(z3 - z1) + y1*(z2 - z3) + y3*(z1 - z2)) + 
			     x2*(y1*(z3 - z0) + y0*(z1 - z3) + y3*(z0 - z1)) + 
			     x3*(y0*(z2 - z1) + y1*(z0 - z2) + y2*(z1 - z0)));

	bc[1] = -(x0*(py*(z3 - z2) + y2*(pz - z3) + y3*(z2 - pz)) + 
		   px*(y2*(z3 - z0) + y0*(z2 - z3) + y3*(z0 - z2)) + 
		   x2*(y0*(z3 - pz) + py*(z0 - z3) + y3*(pz - z0)) + 
		   x3*(py*(z2 - z0) + y0*(pz - z2) + y2*(z0 - pz))) * den;
	if( Num_Approx( bc[1], 0.0 ) ) bc[1] = 0.0;
	else if( Num_Approx( bc[1], 1.0 ) ) bc[1] = 1.0;

	bc[2] = (x0*(py*(z3 - z1) + y1*(pz - z3) + y3*(z1 - pz)) + 
		  px*(y1*(z3 - z0) + y0*(z1 - z3) + y3*(z0 - z1)) + 
		  x1*(y0*(z3 - pz) + py*(z0 - z3) + y3*(pz - z0)) + 
		  x3*(py*(z1 - z0) + y0*(pz - z1) + y1*(z0 - pz))) * den;
	if( Num_Approx( bc[2], 0.0 ) ) bc[2] = 0.0;
	else if( Num_Approx( bc[2], 1.0 ) ) bc[2] = 1.0;

	bc[3] = -(x0*(py*(z2 - z1) + y1*(pz - z2) + y2*(z1 - pz)) + 
		   px*(y1*(z2 - z0) + y0*(z1 - z2) + y2*(z0 - z1)) + 
		   x1*(y0*(z2 - pz) + py*(z0 - z2) + y2*(pz - z0)) + 
		   x2*(py*(z1 - z0) + y0*(pz - z1) + y1*(z0 - pz))) * den;
	if( Num_Approx( bc[3], 0.0 ) ) bc[3] = 0.0;
	else if( Num_Approx( bc[3], 1.0 ) ) bc[3] = 1.0;

	bc[0] = 1.0 - bc[1] - bc[2] - bc[3];
	if( Num_Approx( bc[0], 0.0 ) ) bc[0] = 0.0;
	else if( Num_Approx( bc[0], 1.0 ) ) bc[0] = 1.0;
}

void Simplex_Barycenter2D( double** verts, unsigned* inc, unsigned* inds, double* point, double* bc ) {
	double*	tri[3];
	double	a;
	double	b;
	double	c;
	double	d;
	double	e;
	double	f;

	assert( verts );
	assert( inc );
	assert( inds );
	assert( point );
	assert( bc );
	assert( verts[inc[inds[0]]] );
	assert( verts[inc[inds[1]]] );
	assert( verts[inc[inds[2]]] );

	tri[0] = verts[inc[inds[0]]];
	tri[1] = verts[inc[inds[1]]];
	tri[2] = verts[inc[inds[2]]];
	a = tri[0][0] - tri[2][0];
	b = tri[1][0] - tri[2][0];
	c = tri[2][0] - point[0];
	d = tri[0][1] - tri[2][1];
	e = tri[1][1] - tri[2][1];
	f = tri[2][1] - point[1];

	bc[0] = (b * f - c * e) / (a * e - b * d);
	if( Num_Approx( bc[0], 0.0 ) ) bc[0] = 0.0;
	else if( Num_Approx( bc[0], 1.0 ) ) bc[0] = 1.0;
	bc[1] = (a * f - c * d) / (b * d - a * e);
	if( Num_Approx( bc[1], 0.0 ) ) bc[1] = 0.0;
	else if( Num_Approx( bc[1], 1.0 ) ) bc[1] = 1.0;
	bc[2] = 1.0 - bc[0] - bc[1];
	if( Num_Approx( bc[2], 0.0 ) ) bc[2] = 0.0;
	else if( Num_Approx( bc[2], 1.0 ) ) bc[2] = 1.0;
}

double Simplex_Volume( double** verts, unsigned* inc, unsigned* inds ) {
	static const double	fac = 1.0 / 6.0;
	double			da[3], db[3], dc[3];

	assert( verts );
	assert( inc );
	assert( inds );

	Vec_Sub3D( da, verts[inc[inds[2]]], verts[inc[inds[0]]] );
	Vec_Sub3D( db, verts[inc[inds[2]]], verts[inc[inds[1]]] );
	Vec_Sub3D( dc, verts[inc[inds[2]]], verts[inc[inds[2]]] );
	Vec_Cross3D( db, db, dc );

	return fac * fabs( Vec_Dot3D( da, db ) );
}

double Simplex_Area( double** verts, unsigned* inc, unsigned* inds ) {
	unsigned	a = inc[inds[0]];
	unsigned	b = inc[inds[1]];
	unsigned	c = inc[inds[2]];

	assert( verts );
	assert( inc );
	assert( inds );

	return 0.5 * fabs( verts[a][0] * (verts[c][1] - verts[b][1]) + 
			   verts[b][0] * (verts[a][1] - verts[c][1]) + 
			   verts[c][0] * (verts[b][1] - verts[a][1]) );
}
