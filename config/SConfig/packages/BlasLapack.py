import os
import SConfig

class BlasLapack(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.cmath = self.dependency(SConfig.packages.cmath)
        self.libraries = [['blas', 'lapack'],
                          ['cblas', 'clapack'],
                          ['mkl', 'mkl_lapack']]
        self.frameworks = [['Accelerate']]
        self.use_rpath = True
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

    # Thanks to there not being a C version of Blas/Lapack on edda, we need
    # to be able to search for that installation specifically.
    def generate_locations(self):
        lib_dir = ['/usr/local/IBM_compilers/xlf/9.1/lib64',
                   '/opt/ibmcmp/lib64']
        use_dir = True
        for d in lib_dir:
            if not os.path.exists(d):
                use_dir = False
                break
        for loc in SConfig.Package.generate_locations(self):
            yield loc
            if use_dir:
                yield [loc[0], loc[1], loc[2] + lib_dir, loc[3]]

    def generate_libraries(self, location):
        lib_dir = '/usr/local/IBM_compilers/xlf/9.1/lib64'
        if lib_dir in location[2]:
            yield ['blas', 'lapack', 'xlf90', 'xlfmath', 'xl']
        else:
            for libs in SConfig.Package.generate_libraries(self, location):
                yield libs
