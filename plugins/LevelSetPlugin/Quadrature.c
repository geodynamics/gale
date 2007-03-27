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
** $Id:  $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <StGermain/StGermain.h>
#include "Polynomial.h"
#include "Quadrature.h"


void Quadrature_Init( Quadrature* self, unsigned nDims, unsigned nPnts1D ) {
	assert( nDims > 0 );
	assert( nPnts1D > 0 );

	self->nDims = nDims;
	self->nPnts = (unsigned)pow( (double)nPnts1D, (double)nDims );
	self->nPnts1D = nPnts1D;
	self->basis = Memory_Alloc_Array_Unnamed( unsigned, nDims );
	self->x = Memory_Alloc_2DArray_Unnamed( double, self->nPnts, nDims );
	self->w = Memory_Alloc_Array_Unnamed( double, self->nPnts );
}


void Quadrature_Delete( Quadrature* self ) {
	FreeArray( self->basis );
	FreeArray( self->x );
	FreeArray( self->w );
}


/* Adapted from numerical recipes. */
void Quadrature_GaussLegendre( Quadrature* self, double low, double upp ) {
	const double	relPrec = 3.0e-11;
	double		xm, xl;
	double		root, oldRoot;
	double		p1, p2, p3, pp;
	unsigned	hn;
	unsigned	root_i, root_j;
	unsigned	n = self->nPnts1D;
	double*		x;
	double*		w;
	unsigned*	inds;
	unsigned	pnt_i, dim_i;

	assert( self );


	/*
	** Calculate the Gauss-Legendre set of abscissas and weight points for evaluating a quadrature
	** of the form W(x) = 1.
	*/

	x = Memory_Alloc_Array_Unnamed( double, n );
	w = Memory_Alloc_Array_Unnamed( double, n );

	/* Only need half the roots due to symmetry. */
	hn = (n + 1) / 2;

	/* Some useful values. */
	xm = 0.5 * (upp + low);
	xl = 0.5 * (upp - low);

	for( root_i = 1; root_i <= hn; root_i++ ) {
		/* Our initial approximation to the root. */
		root = cos( M_PI * (root_i - 0.25) / (n + 0.5) );
		do {
			p1 = 1.0;
			p2 = 0.0;
			for( root_j = 1; root_j <= n; root_j++ ) {
				p3 = p2;
				p2 = p1;
				p1 = ((2.0 * root_j - 1.0) * root * p2 - (root_j - 1.0) * p3) / root_j;
			}

			/* Now for the derivative. */
			pp = n * (root * p1 - p2) / (root * root - 1.0);

			/* Update the root. */
			oldRoot = root;
			root = oldRoot - p1 / pp; /* Newton's method. */
		} while( fabs( root - oldRoot ) > relPrec );

		x[root_i - 1] = xm - xl * root;
		x[n - root_i] = xm + xl * root;
		w[root_i - 1] = 2.0 * xl / ((1.0 - root * root) * pp * pp);
		w[n - root_i] = w[root_i - 1];
	}


	/*
	** Extend to multiple dimensions.
	*/

	/* Calculate a basis for the quadrature points. */
	self->basis[0] = 1;
	for( dim_i = 1; dim_i < self->nDims; dim_i++ ) {
		self->basis[dim_i] = self->basis[dim_i - 1] * self->nPnts1D;
	}

	/* Calculate multi-dimensional quadrature weights and abscissas. */
	inds = Memory_Alloc_Array_Unnamed( unsigned, self->nDims );
	for( pnt_i = 0; pnt_i < self->nPnts; pnt_i++ ) {
		/* Calc multi-dimensional indices. */
		inds[self->nDims - 1] = pnt_i / self->basis[self->nDims - 1];
		for( dim_i = self->nDims - 1; dim_i > 0; dim_i-- ) {
			inds[dim_i - 1] = (pnt_i - inds[dim_i] * self->basis[dim_i]) / self->basis[dim_i - 1];
		}

		self->x[pnt_i][dim_i] = x[inds[0]];
		self->w[pnt_i] = w[inds[0]];
		for( dim_i = 1; dim_i < self->nDims; dim_i++ ) {
			self->x[pnt_i][dim_i] = x[inds[dim_i]];
			self->w[pnt_i] *= w[inds[dim_i]];
		}
	}
	FreeArray( inds );
}


double Quadrature_Eval( Quadrature* self, Function* f ) {
	double		res = 0.0;
	unsigned	pnt_i;

	assert( self );
	assert( self->nPnts > 0 );
	assert( self->x );
	assert( self->w );
	assert( f );
	assert( self->nDims == f->nDims );

	for( pnt_i = 0; pnt_i < self->nPnts; pnt_i++ ) {
		res += self->w[pnt_i] * Function_Eval( f, self->x[pnt_i] );
	}

	return res;
}
