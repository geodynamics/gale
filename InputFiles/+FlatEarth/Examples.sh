# Example runs using the 2D and 3D convection apps in this directory

# UWPATH points to a working build of Underworld (with experimental to light up multigrid)

export UWPATH=/home/louismoresi/stgUnderworld-2008.04.01
export UWEXEC=$UWPATH/build-optimised/bin/Underworld


## Basic 2D PIC convection setup using U/W multigrid solver

OUTPUTDIR=output/ConvectionPIC-Ra1e6-Arrh1e6-4proc-64x64; mkdir $OUTPUTDIR >& /dev/null; \
nohup mpirun -np 4 $UWEXEC $UWPATH/Underworld/InputFiles/+FlatEarth/MantleConvection2D+PIC+PLASTICITY+REP.xml \
$UWPATH/Experimental/InputFiles/Experimental_Components/MultigridForRegular.xml \
--outputPath=$OUTPUTDIR --checkpointPath=$OUTPUTDIR --dumpEvery=10 --checkpointEvery=50 \
--nonLinearTolerance=1.0e-3 --nonLinearMaxIterations=25 \
--mgLevels=5 --elementResI=64 --elementResJ=64  \
--maxTimeSteps=25000 \
< /dev/null >& $OUTPUTDIR/output.txt &

## Basic 3D PIC convection setup using U/W multigrid solver

OUTPUTDIR=output/ConvectionPIC-Ra1e6-Arrh1e6-4proc-64x64; mkdir $OUTPUTDIR >& /dev/null; \
nohup mpirun -np 4 $UWEXEC $UWPATH/Underworld/InputFiles/+FlatEarth/MantleConvection3D+PIC+PLASTICITY+REP.xml \
$UWPATH/Experimental/InputFiles/Experimental_Components/MultigridForRegular.xml \
--outputPath=$OUTPUTDIR --checkpointPath=$OUTPUTDIR --dumpEvery=10 --checkpointEvery=250 \
--nonLinearTolerance=1.0e-3 --nonLinearMaxIterations=25 \
--mgLevels=5 --elementResI=64 --elementResJ=32 --elementResJ=64 --maxX=2.0 --maxZ=2.0 \
--maxTimeSteps=25000 \
< /dev/null >& $OUTPUTDIR/output.txt &


# Rerun from checkpoint using PETSCext from Dave May (david.may@sci.monash.edu.au to obtain these solvers)
# $PETSCEXT_OPTS_PARALLEL_MG contains the runtime flags to set up the solver environment

OUTPUTDIR=output/ConvectionPIC-Ra1e6-Arrh1e6-4proc-64x64; mkdir $OUTPUTDIR >& /dev/null; \
nohup mpirun -np 4 $UWEXEC $UWPATH/Underworld/InputFiles/+FlatEarth/MantleConvection2D+PIC+PLASTICITY+REP.xml \
$UWPATH/Experimental/InputFiles/PetscExt-Solvers+mg.xml $PETSCEXT_OPTS_PARALLEL_MG \
--outputPath=$OUTPUTDIR --checkpointPath=$OUTPUTDIR --dumpEvery=10 --checkpointEvery=50 --restartTimestep=50 \
--nonLinearTolerance=1.0e-3 --nonLinearMaxIterations=25 \
--mgLevels=5 --elementResI=64 --elementResJ=64  \
--maxTimeSteps=25000 \
< /dev/null >& $OUTPUTDIR/output3.txt &

