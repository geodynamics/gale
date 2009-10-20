./checkpointTest.pl Anisotropic.xml
./checkpointTest.pl Arrhenius.xml
./checkpointTest.pl ArrheniusPIC.xml
./checkpointTest.pl DepthDependentViscosity.xml
./checkpointTest.pl Extension.xml
./checkpointTest.pl FrankKamenetskii.xml
./checkpointTest.pl InternalHeating.xml
./checkpointTest.pl MobileLid.xml
./checkpointTest.pl MultiComponent.xml
./checkpointTest.pl NonNewtonian.xml
#./checkpointTest.pl PeriodicBenchmark.xml
./checkpointTest.pl RayleighTaylorBenchmark.xml

./analyticTest.pl NonNewtonianShear.xml -optionsFile np-1.dat
./analyticTest.pl NonNewtonianShear.xml -optionsFile np-2.dat
./analyticTest.pl NonNewtonianShear.xml -optionsFile np-4.dat
./analyticTest.pl NonNewtonianShear.xml -optionsFile np-8.dat
./analyticTest.pl Trubitsyn2006Isoviscous.xml -optionsFile np-1.dat
./analyticTest.pl Trubitsyn2006Isoviscous.xml -optionsFile np-2.dat
./analyticTest.pl Trubitsyn2006Isoviscous.xml -optionsFile np-4.dat
./analyticTest.pl Trubitsyn2006Isoviscous.xml -optionsFile np-8.dat
