/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2005, Monash Cluster Computing 
** All rights reserved.
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
** 		* Redistributions of source code must retain the above copyright notice, 
** 			this list of conditions and the following disclaimer.
** 		* Redistributions in binary form must reproduce the above copyright 
**			notice, this list of conditions and the following disclaimer in the 
**			documentation and/or other materials provided with the distribution.
** 		* Neither the name of the Monash University nor the names of its contributors 
**			may be used to endorse or promote products derived from this software 
**			without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
** THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
** PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS 
** BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
** CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
** SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
** HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
** LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
** OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
**
** Contact:
*%		Louis Moresi - Louis.Moresi@sci.monash.edu.au
*%
** Contributors:
*+		Mirko Velic
*+		Julian Giordani
*+		Robert Turnbull
*+		Vincent Lemiale
*+		Louis Moresi
*+		David May
*+		David Stegman
*+		Patrick Sunter
** $Id: solA.c 565 2006-05-19 02:33:01Z RobertTurnbull $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
#include <assert.h>
#include "solKx.h"

const Type Velic_solKx_Type = "Velic_solKx";

double Calculate_AA( Velic_solKx* self ) {
	double t1, t4, t6, t7, t9, t12, t14, t16, t17; 
	
	double    sigma    = self->sigma;
	double    kn       = M_PI * self->kn;
	double    km       = M_PI * self->km;
	double    B        = self->B;

	t1 = kn * kn;
	t4 = km * km;
	t6 = t4 * t4;
	t7 = B * B;
	t9 = 0.4e1 * t7 * t4;
	t12 = 0.8e1 * t7 * kn * km;
	t14 = 0.4e1 * t7 * t1;
	t16 = 0.2e1 * t4 * t1;
	t17 = t1 * t1;
	return -0.4e1 * B * t1 * sigma * (t4 + t1) / (t6 + t9 + t12 + t14 + t16 + t17) / (t6 + t9 - t12 + t14 + t16 + t17);
}

double Calculate_BB( Velic_solKx* self ) {
	double t2, t3, t4, t6, t7, t9, t10, t12, t16; 

	double    sigma    = self->sigma;
	double    kn       = M_PI * self->kn;
	double    km       = M_PI * self->km;
	double    B        = self->B;

	t2 = kn * kn;
	t3 = t2 * t2;
	t4 = B * B;
	t6 = 0.4e1 * t4 * t2;
	t7 = km * km;
	t9 = 0.4e1 * t7 * t4;
	t10 = t7 * t7;
	t12 = 0.2e1 * t7 * t2;
	t16 = 0.8e1 * t4 * kn * km;
	return sigma * kn * (t3 - t6 + t9 + t10 + t12) / (t10 + t9 + t16 + t6 + t12 + t3) / (t10 + t9 - t16 + t6 + t12 + t3);
}

double Calculate__C1( Velic_solKx* self, double x, double z, double y ) {
	double t1, t2, t3, t4, t6, t7, t11, t13, t14, t15, t16, t17, t18, t20, t22, t23, t25, t30, t31, t34, t39, t45, t50, t52, t53, t54, t58, t60, t63, t70, t83, t85, t86, t88, t93, t101, t103, t104, t117, t130, t135, t146, t150;

	double    kn       = M_PI * self->kn;
	double    km       = M_PI * self->km;
	double    B        = self->B;

	double a  = B*B + km*km;
	double b  = 2.0*km*B;
	double r  = sqrt(a*a + b*b);
	double Rp = sqrt( (r+a)/2.0 );
	double Rm = sqrt( (r-a)/2.0 );

	t1 = -B + Rp;
	t2 = Rm * t1;
	t3 = Rm * Rp;
	t4 = Rm * B;
	t6 = 0.2e1 * B * kn;
	t7 = t3 + t4 + t6;
	t11 = Rp * Rp;
	t13 = 0.2e1 * B * Rp;
	t14 = B * B;
	t15 = 0.3e1 * t14;
	t16 = kn * kn;
	t17 = Rm * Rm;
	t18 = t11 + t13 - t15 + t16 - t17;
	t20 = t2 * t18 * Calculate_BB( self );
	t22 = -Rm + kn;
	t23 = cos(t22);
	t25 = t3 + t4 - t6;
	t30 = Rm + kn;
	t31 = cos(t30);
	t34 = t2 * t18 * Calculate_AA( self );
	t39 = sin(t30);
	t45 = sin(t22);
	t50 = exp(-0.3e1 * Rp - B);
	t52 = Rp + B;
	t53 = Rm * t52;
	t54 = t3 - t4 - t6;
	t58 = t11 - t13 - t17 - t15 + t16;
	t60 = t53 * t58 * Calculate_BB( self );
	t63 = t3 - t4 + t6;
	t70 = t53 * t58 * Calculate_AA( self );
	t83 = exp(-t52);
	t85 = t11 * B;
	t86 = t14 * Rp;
	t88 = t17 * Rp;
	t93 = t17 * B;
	t101 = 0.8e1 * t14 * Calculate_BB( self ) * kn * Rp;
	t103 = 0.2e1 * Rm;
	t104 = cos(t103);
	t117 = sin(t103);
	t130 = exp(-0.2e1 * Rp);
	t135 = exp(-0.4e1 * Rp);
	t146 = t17 * t14;
	t150 = t17 * t11;
	return (((-0.2e1 * t2 * t7 * Calculate_AA( self ) + t20) * t23 + (-0.2e1 * t2 * t25 * Calculate_AA( self ) - t20) * t31 + (t34 - 0.2e1 * t2 * t25 * Calculate_BB( self )) * t39 + (-t34 - 0.2e1 * t2 * t7 * Calculate_BB( self )) * t45) * t50 + ((0.2e1 * t53 * t54 * Calculate_AA( self ) + t60) * t23 + (0.2e1 * t53 * t63 * Calculate_AA( self ) - t60) * t31 + (t70 + 0.2e1 * t53 * t63 * Calculate_BB( self )) * t39 + (-t70 + 0.2e1 * t53 * t54 * Calculate_BB( self )) * t45) * t83 + ((-0.2e1 * Rp * (t85 + 0.2e1 * t86 + 0.2e1 * t88 - 0.3e1 * t14 * B + t16 * B + t93) * Calculate_AA( self ) - t101) * t104 + (-0.2e1 * t3 * (t11 - 0.5e1 * t14 + t16 - t17) * Calculate_AA( self ) - 0.8e1 * B * Calculate_BB( self ) * kn * Rm * Rp) * t117 + 0.2e1 * B * (t11 * Rp + 0.2e1 * t85 - 0.3e1 * t86 + t88 + t16 * Rp + 0.2e1 * t93) * Calculate_AA( self ) + t101) * t130 + 0.4e1 * t17 * t1 * t52 * Calculate_AA( self ) * t135) / (((0.8e1 * t14 + 0.8e1 * t17) * t11 * t104 - 0.8e1 * t11 * t14 - 0.8e1 * t146) * t130 + (-0.4e1 * t150 + 0.4e1 * t146) * t135 - 0.4e1 * t150 + 0.4e1 * t146);
}

double Calculate__C2( Velic_solKx* self, double x, double z, double y ) {
	double t1, t2, t4, t5, t6, t7, t8, t10, t11, t12, t13, t14, t15, t16, t17, t18, t20, t21, t22, t24, t25, t28, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t41, t42, t43, t47, t48, t50, t52, t56, t57, t63, t69, t74, t76, t77, t80, t81, t83, t84, t85, t86, t92, t102, t112, t120, t125, t126, t137, t142, t152, t160; 

	double    kn       = M_PI * self->kn;
	double    km       = M_PI * self->km;
	double    B        = self->B;

	double a  = B*B + km*km;
	double b  = 2.0*km*B;
	double r  = sqrt(a*a + b*b);
	double Rp = sqrt( (r+a)/2.0 );
	double Rm = sqrt( (r-a)/2.0 );

	t1 = B * B;
	t2 = t1 * kn;
	t4 = 0.8e1 * t2 * Rp;
	t5 = Rp * Rp;
	t6 = t5 * Rp;
	t7 = Rm * t6;
	t8 = Rp * t1;
	t10 = 0.5e1 * t8 * Rm;
	t11 = kn * kn;
	t12 = t11 * Rp;
	t13 = t12 * Rm;
	t14 = t11 * B;
	t15 = t14 * Rm;
	t16 = Rm * Rm;
	t17 = t16 * Rm;
	t18 = t17 * Rp;
	t20 = Rm * t5 * B;
	t21 = t17 * B;
	t22 = t1 * B;
	t24 = 0.3e1 * t22 * Rm;
	t25 = -t4 + t7 - t10 + t13 + t15 - t18 - t20 - t21 - t24;
	t28 = 0.3e1 * t22 * Rp;
	t30 = 0.2e1 * t2 * Rm;
	t31 = t16 * t1;
	t32 = t16 * t5;
	t33 = t6 * B;
	t34 = t16 * Rp;
	t35 = t34 * B;
	t36 = t5 * t1;
	t37 = 0.2e1 * t36;
	t38 = B * kn;
	t39 = Rm * Rp;
	t41 = 0.2e1 * t38 * t39;
	t42 = t12 * B;
	t43 = -t28 + t30 + t31 + t32 + t33 + t35 + t37 + t41 + t42;
	t47 = -Rm + kn;
	t48 = cos(t47);
	t50 = t13 + t15 - t18 - t24 - t20 - t21 + t7 - t10 + t4;
	t52 = -t30 + t37 + t33 - t28 + t31 + t35 - t41 + t32 + t42;
	t56 = Rm + kn;
	t57 = cos(t56);
	t63 = sin(t56);
	t69 = sin(t47);
	t74 = exp(-0.3e1 * Rp - B);
	t76 = Rp + B;
	t77 = Rm * t76;
	t80 = 0.3e1 * t1;
	t81 = t5 - 0.2e1 * B * Rp - t16 - t80 + t11;
	t83 = t77 * t81 * Calculate_AA( self );
	t84 = Rm * B;
	t85 = 0.2e1 * t38;
	t86 = t39 - t84 - t85;
	t92 = t39 - t84 + t85;
	t102 = t77 * t81 * Calculate_BB( self );
	t112 = exp(-t76);
	t120 = kn * Rm;
	t125 = 0.2e1 * Rm;
	t126 = cos(t125);
	t137 = t1 * Calculate_BB( self );
	t142 = sin(t125);
	t152 = exp(-0.2e1 * Rp);
	t160 = exp(-0.4e1 * Rp);
	return (((t25 * Calculate_AA( self ) + 0.2e1 * t43 * Calculate_BB( self )) * t48 + (t50 * Calculate_AA( self ) - 0.2e1 * t52 * Calculate_BB( self )) * t57 + (0.2e1 * t52 * Calculate_AA( self ) + t50 * Calculate_BB( self )) * t63 + (-0.2e1 * t43 * Calculate_AA( self ) + t25 * Calculate_BB( self )) * t69) * t74 + ((-t83 + 0.2e1 * t77 * t86 * Calculate_BB( self )) * t48 + (-t83 - 0.2e1 * t77 * t92 * Calculate_BB( self )) * t57 + (0.2e1 * t77 * t92 * Calculate_AA( self ) - t102) * t63 + (-0.2e1 * t77 * t86 * Calculate_AA( self ) - t102) * t69) * t112 + ((0.2e1 * t39 * (t5 - 0.5e1 * t1 + t11 - t16) * Calculate_AA( self ) + 0.8e1 * B * Calculate_BB( self ) * t120 * Rp) * t126 + (-0.2e1 * Rp * (t5 * B + 0.2e1 * t8 + 0.2e1 * t34 - 0.3e1 * t22 + t14 + t16 * B) * Calculate_AA( self ) - 0.8e1 * t137 * kn * Rp) * t142 - 0.2e1 * t84 * (t5 + t80 + t16 - t11) * Calculate_AA( self ) + 0.8e1 * t137 * t120) * t152 + (-0.2e1 * t83 - 0.8e1 * t38 * t77 * Calculate_BB( self )) * t160) / (((0.8e1 * t1 + 0.8e1 * t16) * t5 * t126 - 0.8e1 * t36 - 0.8e1 * t31) * t152 + (-0.4e1 * t32 + 0.4e1 * t31) * t160 - 0.4e1 * t32 + 0.4e1 * t31);	
}

double Calculate__C3( Velic_solKx* self, double x, double z, double y ) {
	double t1, t2, t3, t4, t6, t7, t11, t13, t14, t15, t16, t17, t18, t20, t22, t23, t25, t30, t31, t34, t39, t45, t50, t52, t53, t54, t58, t60, t63, t70, t83, t85, t86, t91, t92, t101, t103, t104, t117, t130, t143, t147, t151; 
	
	double    kn       = M_PI * self->kn;
	double    km       = M_PI * self->km;
	double    B        = self->B;

	double a  = B*B + km*km;
	double b  = 2.0*km*B;
	double r  = sqrt(a*a + b*b);
	double Rp = sqrt( (r+a)/2.0 );
	double Rm = sqrt( (r-a)/2.0 );

	t1 = -B + Rp;
	t2 = Rm * t1;
	t3 = Rm * Rp;
	t4 = Rm * B;
	t6 = 0.2e1 * B * kn;
	t7 = t3 + t4 + t6;
	t11 = Rp * Rp;
	t13 = 0.2e1 * B * Rp;
	t14 = B * B;
	t15 = 0.3e1 * t14;
	t16 = kn * kn;
	t17 = Rm * Rm;
	t18 = t11 + t13 - t15 + t16 - t17;
	t20 = t2 * t18 * Calculate_BB( self );
	t22 = -Rm + kn;
	t23 = cos(t22);
	t25 = t3 + t4 - t6;
	t30 = Rm + kn;
	t31 = cos(t30);
	t34 = t2 * t18 * Calculate_AA( self );
	t39 = sin(t30);
	t45 = sin(t22);
	t50 = exp(-0.3e1 * Rp - B);
	t52 = Rp + B;
	t53 = Rm * t52;
	t54 = t3 - t4 - t6;
	t58 = t11 - t13 - t17 - t15 + t16;
	t60 = t53 * t58 * Calculate_BB( self );
	t63 = t3 - t4 + t6;
	t70 = t53 * t58 * Calculate_AA( self );
	t83 = exp(-t52);
	t85 = t17 * B;
	t86 = t17 * Rp;
	t91 = t11 * B;
	t92 = t14 * Rp;
	t101 = 0.8e1 * t14 * Calculate_BB( self ) * kn * Rp;
	t103 = 0.2e1 * Rm;
	t104 = cos(t103);
	t117 = sin(t103);
	t130 = exp(-0.2e1 * Rp);
	t143 = t17 * t14;
	t147 = t17 * t11;
	t151 = exp(-0.4e1 * Rp);
	return (((0.2e1 * t2 * t7 * Calculate_AA( self ) - t20) * t23 + (0.2e1 * t2 * t25 * Calculate_AA( self ) + t20) * t31 + (-t34 + 0.2e1 * t2 * t25 * Calculate_BB( self )) * t39 + (t34 + 0.2e1 * t2 * t7 * Calculate_BB( self )) * t45) * t50 + ((-0.2e1 * t53 * t54 * Calculate_AA( self ) - t60) * t23 + (-0.2e1 * t53 * t63 * Calculate_AA( self ) + t60) * t31 + (-t70 - 0.2e1 * t53 * t63 * Calculate_BB( self )) * t39 + (t70 - 0.2e1 * t53 * t54 * Calculate_BB( self )) * t45) * t83 + ((0.2e1 * Rp * (t85 - 0.2e1 * t86 - 0.3e1 * t14 * B + t16 * B + t91 - 0.2e1 * t92) * Calculate_AA( self ) + t101) * t104 + (0.2e1 * t3 * (t11 - 0.5e1 * t14 + t16 - t17) * Calculate_AA( self ) + 0.8e1 * B * Calculate_BB( self ) * kn * Rm * Rp) * t117 - 0.2e1 * B * (t11 * Rp - 0.2e1 * t91 - 0.3e1 * t92 + t86 + t16 * Rp - 0.2e1 * t85) * Calculate_AA( self ) - t101) * t130 + 0.4e1 * t17 * t1 * t52 * Calculate_AA( self )) / (((0.8e1 * t14 + 0.8e1 * t17) * t11 * t104 - 0.8e1 * t11 * t14 - 0.8e1 * t143) * t130 + (-0.4e1 * t147 + 0.4e1 * t143) * t151 - 0.4e1 * t147 + 0.4e1 * t143);
			
}

double Calculate__C4( Velic_solKx* self, double x, double z, double y ) {
	double t2, t3, t6, t7, t8, t9, t10, t12, t13, t14, t15, t16, t17, t22, t23, t25, t30, t31, t37, t39, t45, t50, t52, t54, t56, t57, t59, t60, t62, t63, t64, t65, t66, t67, t68, t69, t70, t71, t72, t74, t75, t76, t77, t79, t81, t82, t83, t84, t85, t87, t88, t93, t95, t112, t120, t125, t126, t137, t142, t152, t170;
	
	double    kn       = M_PI * self->kn;
	double    km       = M_PI * self->km;
	double    B        = self->B;

	double a  = B*B + km*km;
	double b  = 2.0*km*B;
	double r  = sqrt(a*a + b*b);
	double Rp = sqrt( (r+a)/2.0 );
	double Rm = sqrt( (r-a)/2.0 );

	t2 = Rm * (-B + Rp);
	t3 = Rp * Rp;
	t6 = B * B;
	t7 = 0.3e1 * t6;
	t8 = kn * kn;
	t9 = Rm * Rm;
	t10 = t3 + 0.2e1 * B * Rp - t7 + t8 - t9;
	t12 = t2 * t10 * Calculate_AA( self );
	t13 = Rm * Rp;
	t14 = Rm * B;
	t15 = B * kn;
	t16 = 0.2e1 * t15;
	t17 = t13 + t14 + t16;
	t22 = -Rm + kn;
	t23 = cos(t22);
	t25 = t13 + t14 - t16;
	t30 = Rm + kn;
	t31 = cos(t30);
	t37 = t2 * t10 * Calculate_BB( self );
	t39 = sin(t30);
	t45 = sin(t22);
	t50 = exp(-0.3e1 * Rp - B);
	t52 = Rp * t6;
	t54 = 0.5e1 * Rm * t52;
	t56 = Rm * t3 * B;
	t57 = t6 * B;
	t59 = 0.3e1 * t57 * Rm;
	t60 = t6 * kn;
	t62 = 0.8e1 * t60 * Rp;
	t63 = t9 * Rm;
	t64 = Rp * t63;
	t65 = B * t8;
	t66 = t65 * Rm;
	t67 = B * t63;
	t68 = Rp * t8;
	t69 = t68 * Rm;
	t70 = t3 * Rp;
	t71 = Rm * t70;
	t72 = -t54 + t56 + t59 - t62 - t64 - t66 + t67 + t69 + t71;
	t74 = t70 * B;
	t75 = t9 * t3;
	t76 = t3 * t6;
	t77 = 0.2e1 * t76;
	t79 = 0.3e1 * t57 * Rp;
	t81 = 0.2e1 * t15 * t13;
	t82 = t9 * Rp;
	t83 = t82 * B;
	t84 = t68 * B;
	t85 = t9 * t6;
	t87 = 0.2e1 * t60 * Rm;
	t88 = t74 - t75 - t77 - t79 + t81 + t83 + t84 - t85 - t87;
	t93 = -t66 + t69 - t64 + t59 - t54 + t62 + t56 + t71 + t67;
	t95 = -t77 + t74 + t87 - t85 - t79 + t84 - t75 + t83 - t81;
	t112 = exp(-Rp - B);
	t120 = kn * Rm;
	t125 = 0.2e1 * Rm;
	t126 = cos(t125);
	t137 = t6 * Calculate_BB( self );
	t142 = sin(t125);
	t152 = exp(-0.2e1 * Rp);
	t170 = exp(-0.4e1 * Rp);
	return (((t12 + 0.2e1 * t2 * t17 * Calculate_BB( self )) * t23 + (t12 - 0.2e1 * t2 * t25 * Calculate_BB( self )) * t31 + (0.2e1 * t2 * t25 * Calculate_AA( self ) + t37) * t39 + (-0.2e1 * t2 * t17 * Calculate_AA( self ) + t37) * t45) * t50 + ((-t72 * Calculate_AA( self ) - 0.2e1 * t88 * Calculate_BB( self )) * t23 + (-t93 * Calculate_AA( self ) + 0.2e1 * t95 * Calculate_BB( self )) * t31 + (-0.2e1 * t95 * Calculate_AA( self ) - t93 * Calculate_BB( self )) * t39 + (0.2e1 * t88 * Calculate_AA( self ) - t72 * Calculate_BB( self )) * t45) * t112 + ((-0.2e1 * t13 * (t3 - 0.5e1 * t6 + t8 - t9) * Calculate_AA( self ) - 0.8e1 * B * Calculate_BB( self ) * t120 * Rp) * t126 + (0.2e1 * Rp * (t9 * B - 0.2e1 * t82 - 0.3e1 * t57 + t65 + t3 * B - 0.2e1 * t52) * Calculate_AA( self ) + 0.8e1 * t137 * kn * Rp) * t142 - 0.2e1 * t14 * (t3 + t7 + t9 - t8) * Calculate_AA( self ) + 0.8e1 * t137 * t120) * t152 + 0.2e1 * t12 + 0.8e1 * t15 * t2 * Calculate_BB( self )) / (((0.8e1 * t6 + 0.8e1 * t9) * t3 * t126 - 0.8e1 * t76 - 0.8e1 * t85) * t152 + (-0.4e1 * t75 + 0.4e1 * t85) * t170 - 0.4e1 * t75 + 0.4e1 * t85);
	
}


double Calculate_u1( Velic_solKx* self, double x, double z, double y ) {
	double t2, t3, t4, t6, t11, t18, t19, t20, t22;
	
	double    kn       = M_PI * self->kn;
	double    km       = M_PI * self->km;
	double    B        = self->B;

	double a  = B*B + km*km;
	double b  = 2.0*km*B;
	double r  = sqrt(a*a + b*b);
	double Rp = sqrt( (r+a)/2.0 );
	double Rm = sqrt( (r-a)/2.0 );
	double UU = Rp - B;
	double VV = Rp + B;

	t2 = exp(UU * x);
	t3 = Rm * x;
	t4 = cos(t3);
	t6 = sin(t3);
	t11 = exp(-VV * x);
	t18 = exp(-0.2e1 * x * B);
	t19 = kn * x;
	t20 = cos(t19);
	t22 = sin(t19);
	return -km * (t2 * (Calculate__C1( self, x, y, z ) * t4 + Calculate__C2( self, x, y, z ) * t6) + t11 * (Calculate__C3( self, x, y, z ) * t4 + Calculate__C4( self, x, y, z ) * t6) + t18 * (Calculate_AA( self ) * t20 + Calculate_BB( self ) * t22));

}

double Calculate_u2( Velic_solKx* self, double x, double z, double y ) {
	double t2, t4, t5, t7, t18, t32, t34, t35, t37;
	
	double    kn       = M_PI * self->kn;
	double    km       = M_PI * self->km;
	double    B        = self->B;

	double a  = B*B + km*km;
	double b  = 2.0*km*B;
	double r  = sqrt(a*a + b*b);
	double Rp = sqrt( (r+a)/2.0 );
	double Rm = sqrt( (r-a)/2.0 );
	double UU = Rp - B;
	double VV = Rp + B;

	t2 = exp(UU * x);
	t4 = Rm * x;
	t5 = cos(t4);
	t7 = sin(t4);
	t18 = exp(-VV * x);
	t32 = exp(-0.2e1 * x * B);
	t34 = kn * x;
	t35 = cos(t34);
	t37 = sin(t34);
	return UU * t2 * (Calculate__C1( self, x, y, z ) * t5 + Calculate__C2( self, x, y, z ) * t7) + t2 * (-Calculate__C1( self, x, y, z ) * t7 * Rm + Calculate__C2( self, x, y, z ) * t5 * Rm) - VV * t18 * (Calculate__C3( self, x, y, z ) * t5 + Calculate__C4( self, x, y, z ) * t7) + t18 * (-Calculate__C3( self, x, y, z ) * t7 * Rm + Calculate__C4( self, x, y, z ) * t5 * Rm) - 0.2e1 * B * t32 * (Calculate_AA( self ) * t35 + Calculate_BB( self ) * t37) + t32 * (-Calculate_AA( self ) * t37 * kn + Calculate_BB( self ) * t35 * kn);
	
}

double Calculate_u3( Velic_solKx* self, double x, double z, double y ) {
	double t2, t3, t6, t8, t9, t11, t22, t34, t36, t37, t39;

	double    kn       = M_PI * self->kn;
	double    km       = M_PI * self->km;
	double    B        = self->B;

	double a  = B*B + km*km;
	double b  = 2.0*km*B;
	double r  = sqrt(a*a + b*b);
	double Rp = sqrt( (r+a)/2.0 );
	double Rm = sqrt( (r-a)/2.0 );
	double UU = Rp - B;
	double VV = Rp + B;

	t2 = 0.2e1 * x * B;
	t3 = exp(t2);
	t6 = exp(UU * x);
	t8 = Rm * x;
	t9 = cos(t8);
	t11 = sin(t8);
	t22 = exp(-VV * x);
	t34 = exp(-t2);
	t36 = kn * x;
	t37 = cos(t36);
	t39 = sin(t36);
	return -0.2e1 * t3 * km * (UU * t6 * (Calculate__C1( self, x, y, z ) * t9 + Calculate__C2( self, x, y, z ) * t11) + t6 * (-Calculate__C1( self, x, y, z ) * t11 * Rm + Calculate__C2( self, x, y, z ) * t9 * Rm) - VV * t22 * (Calculate__C3( self, x, y, z ) * t9 + Calculate__C4( self, x, y, z ) * t11) + t22 * (-Calculate__C3( self, x, y, z ) * t11 * Rm + Calculate__C4( self, x, y, z ) * t9 * Rm) - 0.2e1 * B * t34 * (Calculate_AA( self ) * t37 + Calculate_BB( self ) * t39) + t34 * (-Calculate_AA( self ) * t39 * kn + Calculate_BB( self ) * t37 * kn));
}

double Calculate_u4( Velic_solKx* self, double x, double z, double y ) {
	double t2, t3, t4, t6, t7, t8, t9, t10, t11, t12, t15, t16, t17, t18, t20, t21, t22, t23, t24, t25, t26, t30, t41, t46, t61, t73;
	
	double    kn       = M_PI * self->kn;
	double    km       = M_PI * self->km;
	double    B        = self->B;

	double a  = B*B + km*km;
	double b  = 2.0*km*B;
	double r  = sqrt(a*a + b*b);
	double Rp = sqrt( (r+a)/2.0 );
	double Rm = sqrt( (r-a)/2.0 );
	double UU = Rp - B;
	double VV = Rp + B;

	t2 = 0.2e1 * x * B;
	t3 = exp(t2);
	t4 = km * km;
	t6 = exp(UU * x);
	t7 = Rm * x;
	t8 = cos(t7);
	t9 = Calculate__C1( self, x, y, z ) * t8;
	t10 = sin(t7);
	t11 = Calculate__C2( self, x, y, z ) * t10;
	t12 = t9 + t11;
	t15 = exp(-VV * x);
	t16 = Calculate__C3( self, x, y, z ) * t8;
	t17 = Calculate__C4( self, x, y, z ) * t10;
	t18 = t16 + t17;
	t20 = exp(-t2);
	t21 = kn * x;
	t22 = cos(t21);
	t23 = Calculate_AA( self ) * t22;
	t24 = sin(t21);
	t25 = Calculate_BB( self ) * t24;
	t26 = t23 + t25;
	t30 = UU * UU;
	t41 = Rm * Rm;
	t46 = VV * VV;
	t61 = B * B;
	t73 = kn * kn;
	return t3 * (t4 * (t6 * t12 + t15 * t18 + t20 * t26) + t30 * t6 * t12 + 0.2e1 * UU * t6 * (-Calculate__C1( self, x, y, z ) * t10 * Rm + Calculate__C2( self, x, y, z ) * t8 * Rm) + t6 * (-t9 * t41 - t11 * t41) + t46 * t15 * t18 - 0.2e1 * VV * t15 * (-Calculate__C3( self, x, y, z ) * t10 * Rm + Calculate__C4( self, x, y, z ) * t8 * Rm) + t15 * (-t16 * t41 - t17 * t41) + 0.4e1 * t61 * t20 * t26 - 0.4e1 * B * t20 * (-Calculate_AA( self ) * t24 * kn + Calculate_BB( self ) * t22 * kn) + t20 * (-t23 * t73 - t25 * t73));
}

double Calculate_u5( Velic_solKx* self, double x, double z, double y ) {
	double t2, t3, t4, t7, t9, t10, t11, t12, t13, t14, t16, t17, t19, t21, t24, t25, t28, t31, t36, t39, t41, t42, t43, t45, t46, t48, t50, t53, t56, t63, t65, t67, t68, t69, t70, t71, t72, t75, t76, t78, t80, t83, t84, t87, t90, t111, t128; 

	double    kn       = M_PI * self->kn;
	double    km       = M_PI * self->km;
	double    B        = self->B;

	double a  = B*B + km*km;
	double b  = 2.0*km*B;
	double r  = sqrt(a*a + b*b);
	double Rp = sqrt( (r+a)/2.0 );
	double Rm = sqrt( (r-a)/2.0 );
	double UU = Rp - B;
	double VV = Rp + B;

	t2 = 0.2e1 * x * B;
	t3 = exp(t2);
	t4 = UU * UU;
	t7 = exp(UU * x);
	t9 = Rm * x;
	t10 = cos(t9);
	t11 = Calculate__C1( self, x, y, z ) * t10;
	t12 = sin(t9);
	t13 = Calculate__C2( self, x, y, z ) * t12;
	t14 = t11 + t13;
	t16 = t4 * t7;
	t17 = Calculate__C1( self, x, y, z ) * t12;
	t19 = Calculate__C2( self, x, y, z ) * t10;
	t21 = -t17 * Rm + t19 * Rm;
	t24 = UU * t7;
	t25 = Rm * Rm;
	t28 = -t11 * t25 - t13 * t25;
	t31 = t25 * Rm;
	t36 = VV * VV;
	t39 = exp(-VV * x);
	t41 = Calculate__C3( self, x, y, z ) * t10;
	t42 = Calculate__C4( self, x, y, z ) * t12;
	t43 = t41 + t42;
	t45 = t36 * t39;
	t46 = Calculate__C3( self, x, y, z ) * t12;
	t48 = Calculate__C4( self, x, y, z ) * t10;
	t50 = -t46 * Rm + t48 * Rm;
	t53 = VV * t39;
	t56 = -t41 * t25 - t42 * t25;
	t63 = B * B;
	t65 = exp(-t2);
	t67 = kn * x;
	t68 = cos(t67);
	t69 = Calculate_AA( self ) * t68;
	t70 = sin(t67);
	t71 = Calculate_BB( self ) * t70;
	t72 = t69 + t71;
	t75 = t63 * t65;
	t76 = Calculate_AA( self ) * t70;
	t78 = Calculate_BB( self ) * t68;
	t80 = -t76 * kn + t78 * kn;
	t83 = B * t65;
	t84 = kn * kn;
	t87 = -t69 * t84 - t71 * t84;
	t90 = t84 * kn;
	t111 = km * km;
	t128 = t4 * UU * t7 * t14 + 0.3e1 * t16 * t21 + 0.3e1 * t24 * t28 + t7 * (t17 * t31 - t19 * t31) - t36 * VV * t39 * t43 + 0.3e1 * t45 * t50 - 0.3e1 * t53 * t56 + t39 * (t46 * t31 - t48 * t31) - 0.8e1 * t63 * B * t65 * t72 + 0.12e2 * t75 * t80 - 0.6e1 * t83 * t87 + t65 * (t76 * t90 - t78 * t90) + 0.2e1 * B * (t16 * t14 + 0.2e1 * t24 * t21 + t7 * t28 + t45 * t43 - 0.2e1 * t53 * t50 + t39 * t56 + 0.4e1 * t75 * t72 - 0.4e1 * t83 * t80 + t65 * t87) - t111 * (t24 * t14 + t7 * t21 - t53 * t43 + t39 * t50 - 0.2e1 * t83 * t72 + t65 * t80) + 0.2e1 * B * t111 * (t7 * t14 + t39 * t43 + t65 * t72);
	 return -t3 * t128 / km;
}

double Calculate_u6( Velic_solKx* self, double x, double z, double y ) {
	double t2, t3, t6, t8, t9, t11, t22, t34, t36, t37, t39;
	
	double    kn       = M_PI * self->kn;
	double    km       = M_PI * self->km;
	double    B        = self->B;

	double a  = B*B + km*km;
	double b  = 2.0*km*B;
	double r  = sqrt(a*a + b*b);
	double Rp = sqrt( (r+a)/2.0 );
	double Rm = sqrt( (r-a)/2.0 );
	double UU = Rp - B;
	double VV = Rp + B;

	t2 = 0.2e1 * x * B;
	t3 = exp(t2);
	t6 = exp(UU * x);
	t8 = Rm * x;
	t9 = cos(t8);
	t11 = sin(t8);
	t22 = exp(-VV * x);
	t34 = exp(-t2);
	t36 = kn * x;
	t37 = cos(t36);
	t39 = sin(t36);
	return 0.2e1 * t3 * km * (UU * t6 * (Calculate__C1( self, x, y, z ) * t9 + Calculate__C2( self, x, y, z ) * t11) + t6 * (-Calculate__C1( self, x, y, z ) * t11 * Rm + Calculate__C2( self, x, y, z ) * t9 * Rm) - VV * t22 * (Calculate__C3( self, x, y, z ) * t9 + Calculate__C4( self, x, y, z ) * t11) + t22 * (-Calculate__C3( self, x, y, z ) * t11 * Rm + Calculate__C4( self, x, y, z ) * t9 * Rm) - 0.2e1 * B * t34 * (Calculate_AA( self ) * t37 + Calculate_BB( self ) * t39) + t34 * (-Calculate_AA( self ) * t39 * kn + Calculate_BB( self ) * t37 * kn));
}

double Calculate_SS( Velic_solKx* self, double x, double z, double y ) {
	
	double    kn       = M_PI * self->kn;
	double    km       = M_PI * self->km;
	double    B        = self->B;

	double a  = B*B + km*km;
	double b  = 2.0*km*B;
	double r  = sqrt(a*a + b*b);
	double Rp = sqrt( (r+a)/2.0 );
	double Rm = sqrt( (r-a)/2.0 );
	double UU = Rp - B;
	double VV = Rp + B;

	return sin(km*z)*(exp(UU*x)*(Calculate__C1( self, x, y, z )*cos(Rm*x)+Calculate__C2( self, x, y, z )*sin(Rm*x)) + exp(-VV*x)*(Calculate__C3( self, x, y, z )*cos(Rm*x)+Calculate__C4( self, x, y, z )*sin(Rm*x)) + exp(-2*x*B)*(Calculate_AA( self )*cos(kn*x)+Calculate_BB( self )*sin(kn*x)));

}
void Velic_solKx_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure ) {
	Velic_solKx* self = (Velic_solKx*) analyticSolution;
	double x    = coord[ I_AXIS ];
	double z    = coord[ J_AXIS ];
	double y    = coord[ K_AXIS ];

	double    km       = M_PI * self->km;
	
	*pressure = Calculate_u5( self, x, y, z )*cos(km*z);
}

void Velic_solKx_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity ) {
	Velic_solKx* self = (Velic_solKx*) analyticSolution;
	double x    = coord[ I_AXIS ];
	double z    = coord[ J_AXIS ];
	double y    = coord[ K_AXIS ];
	
	double    km       = M_PI * self->km;

	velocity[ I_AXIS ] = Calculate_u1( self, x, y, z ) * cos(km*z);
	velocity[ J_AXIS ] = Calculate_u2( self, x, y, z ) * sin(km*z);
	if ( analyticFeVariable->dim == 3 ) 
		velocity[ K_AXIS ] = 0.0;
	printf("%0.7f %0.7f %0.7f %0.7f \n",x,z,velocity[ I_AXIS ],velocity[ J_AXIS ]);
}

void Velic_solKx_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress ) {
	Velic_solKx* self = (Velic_solKx*) analyticSolution;
	double x    = coord[ I_AXIS ];
	double z    = coord[ J_AXIS ];
	double y    = coord[ K_AXIS ];

	double    km       = M_PI * self->km;
	double tmp;
	double U5   = Calculate_u5( self, x, y, z );
	tmp = Calculate_u6( self, x, y, z ) - U5; 
	stress[0] = tmp*cos(km*z); /* xx stress */
	tmp = Calculate_u3( self, x, y, z ) - U5;
	stress[1] = Calculate_u3( self, x, y, z ) * cos(km*z); /* zz stress */
	stress[2] = Calculate_u4( self, x, y, z ) * sin(km*z); /* zx stress */
	
}
	
void _Velic_solKx_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	Velic_solKx* self = (Velic_solKx*) analyticSolution;
	FeVariable*              velocityField;
	FeVariable*              pressureField;
	FeVariable*              stressField;

	/* Construct Parent */
	_AnalyticSolution_Construct( self, cf, data );

	/* Create Analytic Fields */
	velocityField = Stg_ComponentFactory_ConstructByName( cf, "VelocityField", FeVariable, True, data ); 
	AnalyticSolution_CreateAnalyticVectorField( self, velocityField, Velic_solKx_VelocityFunction );

	pressureField = Stg_ComponentFactory_ConstructByName( cf, "PressureField", FeVariable, True, data ); 
	AnalyticSolution_CreateAnalyticField( self, pressureField, Velic_solKx_PressureFunction );

	stressField = Stg_ComponentFactory_ConstructByName( cf, "StressField", FeVariable, False, data ); 
	if ( stressField )
		AnalyticSolution_CreateAnalyticSymmetricTensorField( self, stressField, Velic_solKx_StressFunction );

	
	self->sigma = Stg_ComponentFactory_GetRootDictDouble( cf, "sigma", 1.0 );
	self->B = Stg_ComponentFactory_GetRootDictDouble( cf, "B", log(100.0)/2.0 );
	self->kn = Stg_ComponentFactory_GetRootDictDouble( cf, "kn", 1.0 );
	self->km = Stg_ComponentFactory_GetRootDictDouble( cf, "km", 1.0 );
}

void* _Velic_solKx_DefaultNew( Name name ) {
	return _AnalyticSolution_New(
			sizeof(Velic_solKx),
			Velic_solKx_Type,
			_AnalyticSolution_Delete,
			_AnalyticSolution_Print,
			_AnalyticSolution_Copy,
			_Velic_solKx_DefaultNew,
			_Velic_solKx_Construct,
			_AnalyticSolution_Build,
			_AnalyticSolution_Initialise,
			_AnalyticSolution_Execute,
			_AnalyticSolution_Destroy,
			name );
}

Index _Velic_solKx_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Velic_solKx_Type, "0", _Velic_solKx_DefaultNew );
}
