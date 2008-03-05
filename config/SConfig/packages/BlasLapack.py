import os
import SConfig

class BlasLapack(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.setup_search_defaults()
        self.setup_options()
        self.libraries = [['blas', 'lapack'],
                          ['cblas', 'clapack'],
                          ['mkl']]
        self.frameworks = ['Accelerate']
        self.symbols = [(['dgeev'], 'FORTRAN_NORMAL'),
                        (['dgeev_'], 'FORTRAN_SINGLE_TRAILINGBAR'),
                        (['dgeev__'], 'FORTRAN_DOUBLE_TRAILINGBAR'),
                        (['DGEEV'], 'FORTRAN_UPPERCASE')]
        self.symbol_setup = '''char jobVecLeft='N';
char jobVecRight='N';
int dim=1;
double* arrayA=NULL;
double* outputReal=NULL;
double* outputImag=NULL;
double* leftEigenVec=NULL;
double* rightEigenVec=NULL;
int leadDimVL=1;
int leadDimVR=1;
double* workSpace=NULL;
int dimWorkSpace;
int INFO=0;

dimWorkSpace=10*dim;
arrayA=malloc(dim*dim*sizeof(double));
memset(arrayA, 0, dim*dim*sizeof(double));
outputReal=malloc(dim*sizeof(double));
outputImag=malloc(dim*sizeof(double));
memset(outputReal, 0, dim*sizeof(double));
memset(outputImag, 0, dim*sizeof(double));
workSpace=malloc(dimWorkSpace*sizeof(double));
leftEigenVec=malloc(leadDimVL*dim*sizeof(double));
rightEigenVec=malloc(leadDimVR*dim*sizeof(double));
memset(leftEigenVec, 0, leadDimVL*dim*sizeof(double));
memset(rightEigenVec, 0, leadDimVR*dim*sizeof(double));
'''
        self.symbol_teardown = '''free(arrayA);
free(outputReal);
free(outputImag);
free(workSpace);
free(leftEigenVec);
free(rightEigenVec);
'''
        self.symbol_prototypes = ['void %s(char*,char*,int*,double*,int*,double*,double*,double*,int*,double*,int*,double*,int*,int*);']
        self.symbol_calls = ['%s(&jobVecLeft, &jobVecRight, &dim, arrayA, &dim, outputReal, outputImag, leftEigenVec, &leadDimVL, rightEigenVec, &leadDimVR, workSpace, &dimWorkSpace, &INFO );']
