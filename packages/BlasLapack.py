import os, platform
import config
import config.utils as utils

class BlasLapack(config.Package):

    def __init__(self, ctx, **kw):
        config.Package.__init__(self, ctx, **kw)
        self.lib_exts = ["em64t"]
	if platform.system() == "Darwin":
            self.base_dirs.append("/System/Library/Frameworks/vecLib.framework")


    def _setup_dependencies(self):
        config.Package._setup_dependencies(self)
        self.syml = self.add_dependency(config.tools.SymLister)


    def setup_libraries(self):
        self.add_library_set([], ["lapack", "blas"])
        self.add_library_set([], ["clapack", "cblas"])
        self.add_library_set([], ["flapack", "fblas"])
        self.add_library_set([], ["f2clapack", "f2cblas"])
        self.add_library_set([], ["mkl_lapack", "mkl"])
        self.add_library_set([], ["guide", "mkl_lapack", "mkl"])
        self.add_library_set([], ["scs"])
        self.add_auxilliary_libs(["c", "f77"], ["m"])
        self.add_auxilliary_libs(["c", "f77"], ["gfortran", "m"])
        self.add_auxilliary_libs(["c", "f77"], ["pthread", "m"])

	if platform.system() == "Darwin":
            self.add_library_set([], ["vecLib"])


    def gen_sources(self, comp, lang):
        c_set = set(["c", "cpp"])
        fort_set = set(["f77", "f90", "f95"])
        comp_set = set(comp["languages"])
        syms = zip(["dgeev", "dgeev_", "dgeev__", "DGEEV"],
                   ["NORMAL", "SINGLE_UNDERBAR", "DOUBLE_UNDERBAR", "CAPS"])

        if lang in fort_set:
            if comp_set.intersection(fort_set):
                # Fortran libraries with fortran compiler.
                yield self.f77_source_code
            elif comp_set.intersection(c_set):
                # Fortran libraries with C compiler.
                for sym, id in syms:
                    yield (self.c_source_code%(sym, sym), "BLASLAPACK_" + id)

        elif lang in c_set:
            if comp_set.intersection(c_set):
                # C libraries with C compiler.
                for sym, id in syms:
                    yield (self.c_source_code%(sym, sym), "BLASLAPACK_" + id)

    c_source_code = """#include <stdlib.h>
#include <stdio.h>
#include <string.h>
void %s(char*,char*,int*,double*,int*,double*,double*,
  double*,int*,double*,int*,double*,int*,int*);
int main( int argc, char** argv ) {
  char jobVecLeft='N';
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
  %s(&jobVecLeft, &jobVecRight, &dim, arrayA, &dim,
    outputReal, outputImag, leftEigenVec, &leadDimVL,
    rightEigenVec, &leadDimVR, workSpace, &dimWorkSpace, &INFO );
  free(arrayA);
  free(outputReal);
  free(outputImag);
  free(workSpace);
  free(leftEigenVec);
  free(rightEigenVec);
  return EXIT_SUCCESS;
}
"""

    f77_source_code = """      program blaslapack
      implicit none
      character jobvl, jobvr
      real a(1,1), wr, wi, vl(1,1), vr(1,1), work(10)
      integer info
      jobvl = 'N'
      jobvr = 'N'
      call dgeev(jobvl, jobvr, 1, a, 1, wr, wi, vl, 1, vr, 1, work,
     % 10, info)
      stop
      end
"""
