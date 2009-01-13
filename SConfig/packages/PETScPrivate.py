import os
import SConfig

class PETScPrivate( SConfig.Package ):
    def __init__( self, scons_env, scons_opts, required = False, **kw ):
        SConfig.Package.__init__( self, scons_env, scons_opts, required, **kw )
        self.petsc = self.dependency( SConfig.packages.PETSc )
        self.checks += [ self.testPetscPrivateSrc ]

    def testPetscPrivateSrc( self ):
	p = os.path.join( self.petsc.base_dir, 'src/mat/impls/aij/mpi/mpiaij.h' )
	if os.path.exists( p ):
            self.env['CPPPATH'].append( self.petsc.base_dir )
            return True
        else:
            return False

