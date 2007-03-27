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
#include <string.h>
#include <assert.h>
#include <StGermain/StGermain.h>
#include "Polynomial.h"


void Polynomial_Clone( Polynomial* p, Polynomial* dst ) {
	assert( p );
	assert( dst );

	dst->order = p->order;
	dst->coefs = Memory_Realloc_Array( dst->coefs, double, dst->order + 1 );
	memcpy( dst->coefs, p->coefs, sizeof(double) * (dst->order + 1) );
}


void Polynomial_Mult( Polynomial* p1, Polynomial* p2, Polynomial* dst ) {
	unsigned	dstOrder;
	unsigned	o_i, o_j;

	assert( p1 );
	assert( p2 );
	assert( dst );

	dst->order = p1->order + p2->order;
	dst->coefs = Memory_Realloc_Array( dst->coefs, double, dst->order + 1 );
	memset( dst->coefs, 0, sizeof(double) * (dst->order + 1) );
	for( o_i = 0; o_i < p1->order + 1; o_i++ ) {
		for( o_j = 0; o_j < p2->order + 1; o_j++ ) {
			dstOrder = o_i + o_j;
			dst->coefs[dstOrder] += p1->coefs[o_i] * p2->coefs[o_j];
		}
	}
}


void Polynomial_ScalarMult( Polynomial* p, double s, Polynomial* dst ) {
	unsigned	o_i;

	assert( p );
	assert( dst );

	dst->order = p->order;
	dst->coefs = Memory_Realloc_Array( dst->coefs, double, dst->order );
	for( o_i = 0; o_i <= p->order; o_i++ ) {
		dst->coefs[o_i] = s * p->coefs[o_i];
	}
}


void Polynomial_Deriv( Polynomial* p, Polynomial* dst ) {
	assert( p );
	assert( dst );

	if( p->order > 0 ) {
		unsigned	o_i;

		dst->order = p->order - 1;
		dst->coefs = Memory_Realloc_Array( dst->coefs, double, dst->order + 1 );
		for( o_i = 1; o_i <= p->order; o_i++ ) {
			dst->coefs[o_i - 1] = p->coefs[o_i] * (double)o_i;
		}
	}
	else {
		dst->order = 0;
		dst->coefs = Memory_Realloc_Array( dst->coefs, double, dst->order + 1 );
		dst->coefs[0] = 0.0;
	}
}


void Polynomial_Integrate( Polynomial* p, Polynomial* dst ) {
	unsigned	o_i;

	assert( p );
	assert( dst );

	dst->order = p->order + 1;
	dst->coefs = Memory_Realloc_Array( dst->coefs, double, dst->order + 1 );
	dst->coefs[0] = 0.0;
	for( o_i = 0; o_i <= p->order; o_i++ ) {
		dst->coefs[o_i + 1] = p->coefs[o_i] / (double)(o_i + 1);
	}
}


double Polynomial_Eval( Polynomial* p, double x ) {
	double		val;
	unsigned	o_i;

	assert( p );
	assert( p->coefs );

	val = p->coefs[0];
	for( o_i = 1; o_i <= p->order; o_i++ ) {
		val += p->coefs[o_i] * pow( x, (double)o_i );
	}

	return val;
}


void Polynomial_Delete( Polynomial* p ) {
	assert( p );
	KillArray( p->coefs );
	p->order = 0;
}


void Function_Mult( Function* f1, Function* f2, Function* dst ) {
	assert( f1 );
	assert( f2 );
	assert( f1->nDims == f2->nDims );
	assert( dst );

	dst->nTerms = f1->nTerms * f2->nTerms;
	dst->nDims = f1->nDims;
	if( dst->nTerms ) {
		unsigned	curTerm;
		unsigned	t_i, t_j, d_i;

		dst->polys = Memory_Realloc_2DArray( dst->polys, Polynomial, dst->nTerms, dst->nDims );
		curTerm = 0;
		for( t_i = 0; t_i < f1->nTerms; t_i++ ) {
			for( t_j = 0; t_j < f2->nTerms; t_j++ ) {
				memset( dst->polys[curTerm], 0, sizeof(Polynomial) * dst->nDims );
				for( d_i = 0; d_i < dst->nDims; d_i++ ) {
					Polynomial_Mult( f1->polys[t_i] + d_i, f2->polys[t_j] + d_i, 
							 dst->polys[curTerm] + d_i );
				}
				curTerm++;
			}
		}
	}
	else {
		dst->polys = NULL;
	}
}


void Function_Partial( Function* f, unsigned dim, Function* dst ) {
	unsigned	t_i, d_i;

	assert( f );
	assert( dim < f->nDims );
	assert( dst );

	dst->nTerms = f->nTerms;
	dst->nDims = f->nDims;
	dst->polys = Memory_Realloc_2DArray( dst->polys, Polynomial, dst->nTerms, dst->nDims );
	for( t_i = 0; t_i < f->nTerms; t_i++ ) {
		memset( dst->polys[t_i], 0, sizeof(Polynomial) * dst->nDims );
		for( d_i = 0; d_i < dst->nDims; d_i++ ) {
			if( d_i == dim ) {
				Polynomial_Deriv( f->polys[t_i] + d_i, dst->polys[t_i] + d_i );
			}
			else {
				Polynomial_Clone( f->polys[t_i] + d_i, dst->polys[t_i] + d_i );
			}
		}
	}
}


void Function_Grad( Function* f, Function* dstVec ) {
	unsigned	d_i;

	assert( f );
	assert( dstVec );

	for( d_i = 0; d_i < f->nDims; d_i++ ) {
		Function_Partial( f, d_i, dstVec + d_i );
	}
}


void Function_Dot( Function* fvec, double* svec, Function* dst ) {
	unsigned	curTerm;
	unsigned	d_i, d_j, t_i;

	assert( fvec );
	assert( svec );
	assert( dst );

	dst->nTerms = 0;
	dst->nDims = fvec[0].nDims;
	for( d_i = 0; d_i < dst->nDims; d_i++ ) {
		dst->nTerms += fvec[d_i].nTerms;
	}
	dst->polys = Memory_Realloc_2DArray( dst->polys, Polynomial, dst->nTerms, dst->nDims );

	curTerm = 0;
	for( d_i = 0; d_i < dst->nDims; d_i++ ) {
		for( t_i = 0; t_i < fvec->nTerms; t_i++ ) {
			memset( dst->polys[curTerm], 0, sizeof(Polynomial) * dst->nDims );
			Polynomial_ScalarMult( fvec[d_i].polys[t_i], svec[d_i], dst->polys[curTerm] );
			for( d_j = 1; d_j < dst->nDims; d_j++ ) {
				Polynomial_Clone( fvec[d_i].polys[t_i] + d_j, dst->polys[curTerm] + d_j );
			}
			curTerm++;
		}
	}
}


void Function_Integrate( Function* f, Function* dst ) {
	unsigned	t_i;

	assert( f );
	assert( dst );

	dst->nTerms = f->nTerms;
	dst->nDims = f->nDims;
	dst->polys = Memory_Realloc_2DArray( dst->polys, Polynomial, dst->nTerms, dst->nDims );
	for( t_i = 0; t_i < f->nTerms; t_i++ ) {
		unsigned	d_i;

		memset( dst->polys[t_i], 0, sizeof(Polynomial) * dst->nDims );
		for( d_i = 0; d_i < dst->nDims; d_i++ ) {
			Polynomial_Integrate( f->polys[t_i] + d_i, dst->polys[t_i] + d_i );
		}
	}
}


double Function_Eval( Function* f, double* x ) {
	double		val, pVal;
	unsigned	t_i;

	assert( f );
	assert( !f->nTerms || f->polys );
	assert( x );

	val = 0.0;
	for( t_i = 0; t_i < f->nTerms; t_i++ ) {
		unsigned	d_i;

		pVal = 1.0;
		for( d_i = 0; d_i < f->nDims; d_i++ ) {
			pVal *= Polynomial_Eval( f->polys[t_i] + d_i, x[d_i] );
		}
		val += pVal;
	}

	return val;
}


void Function_Delete( Function* f ) {
	unsigned	term_i;

	for( term_i = 0; term_i < f->nTerms; term_i++ ) {
		unsigned	dim_i;

		for( dim_i = 0; dim_i < f->nDims; dim_i++ ) {
			Polynomial_Delete( f->polys[term_i] + dim_i );
		}
	}
	KillArray( f->polys );
	f->nTerms = 0;
	f->nDims = 0;
}
