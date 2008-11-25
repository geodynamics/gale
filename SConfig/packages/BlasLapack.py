import os
import SConfig

class BlasLapack(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Package.__init__(self, scons_env, scons_opts, required)
        self.cmath = self.dependency(SConfig.packages.cmath)
        self.libraries = [['lapack', 'blas'],
                          ['flapack', 'fblas'],
                          ['flapack', 'fblas', 'gfortran'],
                          ['clapack', 'cblas'],
                          ['f2clapack', 'f2cblas'],
                          ['mkl_lapack', 'mkl']]
        self.shared_libraries = ['lapack', 'blas']
        self.extra_libraries = ['gfortran']
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

    # Thanks to there not being a C version of Blas/Lapack on edda, we need
    # to be able to search for that installation specifically.
    def process_installations(self, inst):
        lib_dir = ['/usr/local/IBM_compilers/xlf/9.1/lib64',
                   '/opt/ibmcmp/lib64']
        use_dir = True
        for d in lib_dir:
            if not os.path.exists(d):
                use_dir = False
                break
        if use_dir:
            inst.add_lib_dirs(lib_dir)

    def generate_libraries(self, inst):
        lib_dir = '/usr/local/IBM_compilers/xlf/9.1/lib64'
        if lib_dir in inst.lib_dirs:
            old_libs = list(inst.extra_libs)
            inst.extra_libs += ['xlf90', 'xlfmath', 'xl']
            yield ['blas', 'lapack', 'xlf90', 'xlfmath', 'xl']
            inst.extra_libs = old_libs
        else:
            for libs in SConfig.Package.generate_libraries(self, inst):
                yield libs
