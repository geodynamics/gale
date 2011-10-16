#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

void _Velic_solKx( 
		const double pos[],
		double _sigma,
		double _m, int _n, 
		double _B,
		double vel[], double* presssure, 
		double total_stress[], double strain_rate[], double* viscosity );


int main( int argc, char **argv )
{
	int i,j;
	double pos[2], vel[2], pressure, total_stress[3], strain_rate[3];
	double x,z;
	
	for (i=0;i<101;i++){
		for(j=0;j<101;j++){
			x = i/100.0;
			z = j/100.0;
			
			pos[0] = x;
			pos[1] = z;
			_Velic_solKx(
					pos,
					1.0,
					(double)1, 1,
					log(100.0)/2.0,
					vel, &pressure, total_stress, strain_rate, NULL );
			printf("%0.7f %0.7f %0.7f %0.7f %0.7f %0.7f %0.7f %0.7f %0.7f %0.7f %0.7f \n",
					pos[0],pos[1],
					vel[0],vel[1], pressure, 
					total_stress[0], total_stress[1], total_stress[2], 
					strain_rate[0], strain_rate[1], strain_rate[2] );
		}
		printf("\n");
	}
	
	return 0;
}




/* this solution matches the stream function values in Zhong's paper for the few given viscosity constrasts */
void _Velic_solKx( 
		const double pos[],
		double _sigma, /* density */
		double _m, int _n, /* wavelength in z, wavenumber in x */
		double _B, /* viscosity parameter */
		double vel[], double* presssure, 
		double total_stress[], double strain_rate[], double* viscosity )
{
	double Z;
	double u1,u2,u3,u4,u5,u6;
	double sum1,sum2,sum3,sum4,sum5,sum6,mag,sum7,x,z;
	double sigma;
	int n;
	double kn;
	double _C1,_C2,_C3,_C4;
	double B, Rp, UU, VV;
	double rho,a,b,r,_aa,_bb,AA,BB,Rm,km,SS;
	
	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;
	double t11,t12,t13,t14,t15,t16,t17,t18,t19,t20;
	double t21,t22,t23,t24,t25,t26,t28,t30,t31,t32;
	double t33,t34,t35,t36,t37,t38,t39,t41,t42,t43;
	double t45,t46,t47,t48,t50,t52,t53,t54,t56,t57;
	double t58,t59,t60,t61,t62,t63,t64,t65,t66,t67;
	double t68,t69,t70,t71,t72,t73,t74,t75,t76,t77;
	double t78,t79,t80,t81,t82,t83,t84,t85,t86,t87;
	double t88,t90,t91,t92,t93,t95,t101,t102,t103,t104;
	double t111,t112,t117,t120,t125,t126,t128,t130,t135,t137;
	double t142,t143,t146,t147,t150,t151,t152,t160,t170;
	
	
	
	/*************************************************************************/
	/*************************************************************************/
	/* rho = -sin(km*z)*cos(kn*x) */ 
	Z = exp( 2.0 * B * x );
	B = _B; /* viscosity parameter must be non-zero*/
	km = _m*M_PI; /* solution valid for km not zero -- should get trivial solution if km=0 */
	n = _n; /* solution valid for n not zero */
	kn = (double)_n*M_PI;
	sigma = _sigma;
	/*************************************************************************/
	/*************************************************************************/
	a = B*B + km*km;
	b = 2.0*km*B;
	r = sqrt(a*a + b*b);
	Rp = sqrt( (r+a)/2.0 );
	Rm  = sqrt( (r-a)/2.0 );
	UU  = Rp - B;
	VV = Rp + B;
	
	x = pos[0];
	z = pos[1];
	
	sum1=0.0;
	sum2=0.0;
	sum3=0.0;
	sum4=0.0;
	sum5=0.0;
	sum6=0.0;
	sum7=0.0;
	
	
	/*******************************************/
	/*         calculate the constants         */
	/*******************************************/
	
	
	t1 = kn * kn;
	t4 = km * km;
	t6 = t4 * t4;
	t7 = B * B;
	t9 = 0.4e1 * t7 * t4;
	t12 = 0.8e1 * t7 * kn * km;
	t14 = 0.4e1 * t7 * t1;
	t16 = 0.2e1 * t4 * t1;
	t17 = t1 * t1;
	_aa = -0.4e1 * B * t1 * sigma * (t4 + t1) / (t6 + t9 + t12 + t14 + t16 + t17) / (t6 + t9 - t12 + t14 + t16 + t17);
	
	t2 = kn * kn;
	t3 = t2 * t2;
	t4 = B * B;
	t6 = 0.4e1 * t4 * t2;
	t7 = km * km;
	t9 = 0.4e1 * t7 * t4;
	t10 = t7 * t7;
	t12 = 0.2e1 * t7 * t2;
	t16 = 0.8e1 * t4 * kn * km;
	_bb = sigma * kn * (t3 - t6 + t9 + t10 + t12) / (t10 + t9 + t16 + t6 + t12 + t3) / (t10 + t9 - t16 + t6 + t12 + t3);
	
	AA = _aa;
	BB = _bb;
	
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
	t20 = t2 * t18 * BB;
	t22 = -Rm + kn;
	t23 = cos(t22);
	t25 = t3 + t4 - t6;
	t30 = Rm + kn;
	t31 = cos(t30);
	t34 = t2 * t18 * AA;
	t39 = sin(t30);
	t45 = sin(t22);
	t50 = exp(-0.3e1 * Rp - B);
	t52 = Rp + B;
	t53 = Rm * t52;
	t54 = t3 - t4 - t6;
	t58 = t11 - t13 - t17 - t15 + t16;
	t60 = t53 * t58 * BB;
	t63 = t3 - t4 + t6;
	t70 = t53 * t58 * AA;
	t83 = exp(-t52);
	t85 = t11 * B;
	t86 = t14 * Rp;
	t88 = t17 * Rp;
	t93 = t17 * B;
	t101 = 0.8e1 * t14 * BB * kn * Rp;
	t103 = 0.2e1 * Rm;
	t104 = cos(t103);
	t117 = sin(t103);
	t130 = exp(-0.2e1 * Rp);
	t135 = exp(-0.4e1 * Rp);
	t146 = t17 * t14;
	t150 = t17 * t11;
	_C1 = (((-0.2e1 * t2 * t7 * AA + t20) * t23 + (-0.2e1 * t2 * t25 * AA - t20) * t31 + (t34 - 0.2e1 * t2 * t25 * BB) * t39 + (-t34 - 0.2e1 * t2 * t7 * BB) * t45) * t50 + ((0.2e1 * t53 * t54 * AA + t60) * t23 + (0.2e1 * t53 * t63 * AA - t60) * t31 + (t70 + 0.2e1 * t53 * t63 * BB) * t39 + (-t70 + 0.2e1 * t53 * t54 * BB) * t45) * t83 + ((-0.2e1 * Rp * (t85 + 0.2e1 * t86 + 0.2e1 * t88 - 0.3e1 * t14 * B + t16 * B + t93) * AA - t101) * t104 + (-0.2e1 * t3 * (t11 - 0.5e1 * t14 + t16 - t17) * AA - 0.8e1 * B * BB * kn * Rm * Rp) * t117 + 0.2e1 * B * (t11 * Rp + 0.2e1 * t85 - 0.3e1 * t86 + t88 + t16 * Rp + 0.2e1 * t93) * AA + t101) * t130 + 0.4e1 * t17 * t1 * t52 * AA * t135) / (((0.8e1 * t14 + 0.8e1 * t17) * t11 * t104 - 0.8e1 * t11 * t14 - 0.8e1 * t146) * t130 + (-0.4e1 * t150 + 0.4e1 * t146) * t135 - 0.4e1 * t150 + 0.4e1 * t146);
	
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
	t83 = t77 * t81 * AA;
	t84 = Rm * B;
	t85 = 0.2e1 * t38;
	t86 = t39 - t84 - t85;
	t92 = t39 - t84 + t85;
	t102 = t77 * t81 * BB;
	t112 = exp(-t76);
	t120 = kn * Rm;
	t125 = 0.2e1 * Rm;
	t126 = cos(t125);
	t137 = t1 * BB;
	t142 = sin(t125);
	t152 = exp(-0.2e1 * Rp);
	t160 = exp(-0.4e1 * Rp);
	_C2 = (((t25 * AA + 0.2e1 * t43 * BB) * t48 + (t50 * AA - 0.2e1 * t52 * BB) * t57 + (0.2e1 * t52 * AA + t50 * BB) * t63 + (-0.2e1 * t43 * AA + t25 * BB) * t69) * t74 + ((-t83 + 0.2e1 * t77 * t86 * BB) * t48 + (-t83 - 0.2e1 * t77 * t92 * BB) * t57 + (0.2e1 * t77 * t92 * AA - t102) * t63 + (-0.2e1 * t77 * t86 * AA - t102) * t69) * t112 + ((0.2e1 * t39 * (t5 - 0.5e1 * t1 + t11 - t16) * AA + 0.8e1 * B * BB * t120 * Rp) * t126 + (-0.2e1 * Rp * (t5 * B + 0.2e1 * t8 + 0.2e1 * t34 - 0.3e1 * t22 + t14 + t16 * B) * AA - 0.8e1 * t137 * kn * Rp) * t142 - 0.2e1 * t84 * (t5 + t80 + t16 - t11) * AA + 0.8e1 * t137 * t120) * t152 + (-0.2e1 * t83 - 0.8e1 * t38 * t77 * BB) * t160) / (((0.8e1 * t1 + 0.8e1 * t16) * t5 * t126 - 0.8e1 * t36 - 0.8e1 * t31) * t152 + (-0.4e1 * t32 + 0.4e1 * t31) * t160 - 0.4e1 * t32 + 0.4e1 * t31);
	
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
	t20 = t2 * t18 * BB;
	t22 = -Rm + kn;
	t23 = cos(t22);
	t25 = t3 + t4 - t6;
	t30 = Rm + kn;
	t31 = cos(t30);
	t34 = t2 * t18 * AA;
	t39 = sin(t30);
	t45 = sin(t22);
	t50 = exp(-0.3e1 * Rp - B);
	t52 = Rp + B;
	t53 = Rm * t52;
	t54 = t3 - t4 - t6;
	t58 = t11 - t13 - t17 - t15 + t16;
	t60 = t53 * t58 * BB;
	t63 = t3 - t4 + t6;
	t70 = t53 * t58 * AA;
	t83 = exp(-t52);
	t85 = t17 * B;
	t86 = t17 * Rp;
	t91 = t11 * B;
	t92 = t14 * Rp;
	t101 = 0.8e1 * t14 * BB * kn * Rp;
	t103 = 0.2e1 * Rm;
	t104 = cos(t103);
	t117 = sin(t103);
	t130 = exp(-0.2e1 * Rp);
	t143 = t17 * t14;
	t147 = t17 * t11;
	t151 = exp(-0.4e1 * Rp);
	_C3 = (((0.2e1 * t2 * t7 * AA - t20) * t23 + (0.2e1 * t2 * t25 * AA + t20) * t31 + (-t34 + 0.2e1 * t2 * t25 * BB) * t39 + (t34 + 0.2e1 * t2 * t7 * BB) * t45) * t50 + ((-0.2e1 * t53 * t54 * AA - t60) * t23 + (-0.2e1 * t53 * t63 * AA + t60) * t31 + (-t70 - 0.2e1 * t53 * t63 * BB) * t39 + (t70 - 0.2e1 * t53 * t54 * BB) * t45) * t83 + ((0.2e1 * Rp * (t85 - 0.2e1 * t86 - 0.3e1 * t14 * B + t16 * B + t91 - 0.2e1 * t92) * AA + t101) * t104 + (0.2e1 * t3 * (t11 - 0.5e1 * t14 + t16 - t17) * AA + 0.8e1 * B * BB * kn * Rm * Rp) * t117 - 0.2e1 * B * (t11 * Rp - 0.2e1 * t91 - 0.3e1 * t92 + t86 + t16 * Rp - 0.2e1 * t85) * AA - t101) * t130 + 0.4e1 * t17 * t1 * t52 * AA) / (((0.8e1 * t14 + 0.8e1 * t17) * t11 * t104 - 0.8e1 * t11 * t14 - 0.8e1 * t143) * t130 + (-0.4e1 * t147 + 0.4e1 * t143) * t151 - 0.4e1 * t147 + 0.4e1 * t143);
	
	t2 = Rm * (-B + Rp);
	t3 = Rp * Rp;
	t6 = B * B;
	t7 = 0.3e1 * t6;
	t8 = kn * kn;
	t9 = Rm * Rm;
	t10 = t3 + 0.2e1 * B * Rp - t7 + t8 - t9;
	t12 = t2 * t10 * AA;
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
	t37 = t2 * t10 * BB;
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
	t137 = t6 * BB;
	t142 = sin(t125);
	t152 = exp(-0.2e1 * Rp);
	t170 = exp(-0.4e1 * Rp);
	_C4 = (((t12 + 0.2e1 * t2 * t17 * BB) * t23 + (t12 - 0.2e1 * t2 * t25 * BB) * t31 + (0.2e1 * t2 * t25 * AA + t37) * t39 + (-0.2e1 * t2 * t17 * AA + t37) * t45) * t50 + ((-t72 * AA - 0.2e1 * t88 * BB) * t23 + (-t93 * AA + 0.2e1 * t95 * BB) * t31 + (-0.2e1 * t95 * AA - t93 * BB) * t39 + (0.2e1 * t88 * AA - t72 * BB) * t45) * t112 + ((-0.2e1 * t13 * (t3 - 0.5e1 * t6 + t8 - t9) * AA - 0.8e1 * B * BB * t120 * Rp) * t126 + (0.2e1 * Rp * (t9 * B - 0.2e1 * t82 - 0.3e1 * t57 + t65 + t3 * B - 0.2e1 * t52) * AA + 0.8e1 * t137 * kn * Rp) * t142 - 0.2e1 * t14 * (t3 + t7 + t9 - t8) * AA + 0.8e1 * t137 * t120) * t152 + 0.2e1 * t12 + 0.8e1 * t15 * t2 * BB) / (((0.8e1 * t6 + 0.8e1 * t9) * t3 * t126 - 0.8e1 * t76 - 0.8e1 * t85) * t152 + (-0.4e1 * t75 + 0.4e1 * t85) * t170 - 0.4e1 * t75 + 0.4e1 * t85);
	
	
	/*******************************************/
	/*       calculate the velocities etc      */
	/*******************************************/
	
	
	t2 = exp(UU * x);
	t3 = Rm * x;
	t4 = cos(t3);
	t6 = sin(t3);
	t11 = exp(-VV * x);
	t18 = exp(-0.2e1 * x * B);
	t19 = kn * x;
	t20 = cos(t19);
	t22 = sin(t19);
	u1 = -km * (t2 * (_C1 * t4 + _C2 * t6) + t11 * (_C3 * t4 + _C4 * t6) + t18 * (AA * t20 + BB * t22));
	
	t2 = exp(UU * x);
	t4 = Rm * x;
	t5 = cos(t4);
	t7 = sin(t4);
	t18 = exp(-VV * x);
	t32 = exp(-0.2e1 * x * B);
	t34 = kn * x;
	t35 = cos(t34);
	t37 = sin(t34);
	u2 = UU * t2 * (_C1 * t5 + _C2 * t7) + t2 * (-_C1 * t7 * Rm + _C2 * t5 * Rm) - VV * t18 * (_C3 * t5 + _C4 * t7) + t18 * (-_C3 * t7 * Rm + _C4 * t5 * Rm) - 0.2e1 * B * t32 * (AA * t35 + BB * t37) + t32 * (-AA * t37 * kn + BB * t35 * kn);
	
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
	u3 = -0.2e1 * t3 * km * (UU * t6 * (_C1 * t9 + _C2 * t11) + t6 * (-_C1 * t11 * Rm + _C2 * t9 * Rm) - VV * t22 * (_C3 * t9 + _C4 * t11) + t22 * (-_C3 * t11 * Rm + _C4 * t9 * Rm) - 0.2e1 * B * t34 * (AA * t37 + BB * t39) + t34 * (-AA * t39 * kn + BB * t37 * kn));
	
	t2 = 0.2e1 * x * B;
	t3 = exp(t2);
	t4 = km * km;
	t6 = exp(UU * x);
	t7 = Rm * x;
	t8 = cos(t7);
	t9 = _C1 * t8;
	t10 = sin(t7);
	t11 = _C2 * t10;
	t12 = t9 + t11;
	t15 = exp(-VV * x);
	t16 = _C3 * t8;
	t17 = _C4 * t10;
	t18 = t16 + t17;
	t20 = exp(-t2);
	t21 = kn * x;
	t22 = cos(t21);
	t23 = AA * t22;
	t24 = sin(t21);
	t25 = BB * t24;
	t26 = t23 + t25;
	t30 = UU * UU;
	t41 = Rm * Rm;
	t46 = VV * VV;
	t61 = B * B;
	t73 = kn * kn;
	u4 = t3 * (t4 * (t6 * t12 + t15 * t18 + t20 * t26) + t30 * t6 * t12 + 0.2e1 * UU * t6 * (-_C1 * t10 * Rm + _C2 * t8 * Rm) + t6 * (-t9 * t41 - t11 * t41) + t46 * t15 * t18 - 0.2e1 * VV * t15 * (-_C3 * t10 * Rm + _C4 * t8 * Rm) + t15 * (-t16 * t41 - t17 * t41) + 0.4e1 * t61 * t20 * t26 - 0.4e1 * B * t20 * (-AA * t24 * kn + BB * t22 * kn) + t20 * (-t23 * t73 - t25 * t73));
	
	t2 = 0.2e1 * x * B;
	t3 = exp(t2);
	t4 = UU * UU;
	t7 = exp(UU * x);
	t9 = Rm * x;
	t10 = cos(t9);
	t11 = _C1 * t10;
	t12 = sin(t9);
	t13 = _C2 * t12;
	t14 = t11 + t13;
	t16 = t4 * t7;
	t17 = _C1 * t12;
	t19 = _C2 * t10;
	t21 = -t17 * Rm + t19 * Rm;
	t24 = UU * t7;
	t25 = Rm * Rm;
	t28 = -t11 * t25 - t13 * t25;
	t31 = t25 * Rm;
	t36 = VV * VV;
	t39 = exp(-VV * x);
	t41 = _C3 * t10;
	t42 = _C4 * t12;
	t43 = t41 + t42;
	t45 = t36 * t39;
	t46 = _C3 * t12;
	t48 = _C4 * t10;
	t50 = -t46 * Rm + t48 * Rm;
	t53 = VV * t39;
	t56 = -t41 * t25 - t42 * t25;
	t63 = B * B;
	t65 = exp(-t2);
	t67 = kn * x;
	t68 = cos(t67);
	t69 = AA * t68;
	t70 = sin(t67);
	t71 = BB * t70;
	t72 = t69 + t71;
	t75 = t63 * t65;
	t76 = AA * t70;
	t78 = BB * t68;
	t80 = -t76 * kn + t78 * kn;
	t83 = B * t65;
	t84 = kn * kn;
	t87 = -t69 * t84 - t71 * t84;
	t90 = t84 * kn;
	t111 = km * km;
	t128 = t4 * UU * t7 * t14 + 0.3e1 * t16 * t21 + 0.3e1 * t24 * t28 + t7 * (t17 * t31 - t19 * t31) - t36 * VV * t39 * t43 + 0.3e1 * t45 * t50 - 0.3e1 * t53 * t56 + t39 * (t46 * t31 - t48 * t31) - 0.8e1 * t63 * B * t65 * t72 + 0.12e2 * t75 * t80 - 0.6e1 * t83 * t87 + t65 * (t76 * t90 - t78 * t90) + 0.2e1 * B * (t16 * t14 + 0.2e1 * t24 * t21 + t7 * t28 + t45 * t43 - 0.2e1 * t53 * t50 + t39 * t56 + 0.4e1 * t75 * t72 - 0.4e1 * t83 * t80 + t65 * t87) - t111 * (t24 * t14 + t7 * t21 - t53 * t43 + t39 * t50 - 0.2e1 * t83 * t72 + t65 * t80) + 0.2e1 * B * t111 * (t7 * t14 + t39 * t43 + t65 * t72);
	u5 = -t3 * t128 / km;
	
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
	u6 = 0.2e1 * t3 * km * (UU * t6 * (_C1 * t9 + _C2 * t11) + t6 * (-_C1 * t11 * Rm + _C2 * t9 * Rm) - VV * t22 * (_C3 * t9 + _C4 * t11) + t22 * (-_C3 * t11 * Rm + _C4 * t9 * Rm) - 0.2e1 * B * t34 * (AA * t37 + BB * t39) + t34 * (-AA * t39 * kn + BB * t37 * kn));
	
	
	
	
	
	SS = sin(km*z)*(exp(UU*x)*(_C1*cos(Rm*x)+_C2*sin(Rm*x)) + exp(-VV*x)*(_C3*cos(Rm*x)+_C4*sin(Rm*x)) + exp(-2*x*B)*(AA*cos(kn*x)+BB*sin(kn*x)));
	
	/* u1 = Vx, u2 = Vz, u3 = txx, u4 = tzx, u5 = pressure, u6 = tzz */
	
	
	sum5 += u5*cos(km*z);  /* pressure */
	u6 -= u5; /* get total stress */
	sum6 += u6*cos(km*z);  /* xx stress */
	
	u1 *= cos(km*z); /* x velocity */
	sum1 += u1;
	u2 *= sin(km*z); /* z velocity */
	sum2 += u2;
	u3 -= u5; /* get total stress */
	u3 *= cos(km*z); /* zz stress */
	sum3 += u3;
	u4 *= sin(km*z); /* zx stress */
	sum4 += u4;
	
	rho = -sigma*sin(km*z)*cos(kn*x); /* density */
	sum7 += rho;
	
	
	mag=sqrt(u1*u1+u2*u2);
	//printf("%0.7f %0.7f %0.7f %0.7f %0.7f %0.7f\n",x,z,u1,u2,mag,SS);
	/*printf("%0.7f %0.7f %0.7f %0.7f %0.7f %0.7f %0.7f %0.7f %0.7f %0.7f %0.7f\n",x,z,sum1,sum2,sum3,sum4,sum5,sum6,mag,sum7,SS);*/
	
	if ( viscosity != NULL ) {
		*viscosity = Z;
	}
	
	/* Output */
	if( vel != NULL ) {
		vel[0] = sum2;
		vel[1] = sum1;
	}
	if( presssure != NULL ) {
		(*presssure) = sum5;
	}
	if( total_stress != NULL ) {
		total_stress[0] = sum6;
		total_stress[1] = sum3;
		total_stress[2] = sum4;
	}
	if( strain_rate != NULL ) {
		/* sigma = tau - p, tau = sigma + p, tau[] = 2*eta*strain_rate[] */
		strain_rate[0] = (sum6+sum5)/(2.0*Z);
		strain_rate[1] = (sum3+sum5)/(2.0*Z);
		strain_rate[2] = (sum4)/(2.0*Z);
	}
	/* Value checks, could be cleaned up if needed. Julian Giordani 9-Oct-2006*/
        if( fabs( sum5 - ( -0.5*(sum6+sum3) ) ) > 1e-5 ) {
                assert(0);
        }
}


