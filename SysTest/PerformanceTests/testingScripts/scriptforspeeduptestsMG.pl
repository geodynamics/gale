#! /usr/bin/perl
$namevar = 'ThermalConvectionMG';
$filenamevar = $namevar . '.xml';
$procvar = 1;
$walltimevar = 10;
$resIvar = 32;
$resJvar = 32;
$mgLevelsvar = 1;
$stoptime = 5000;
$countervar = 1;
print "before the first while loop\n";
#going up to res 128x128
while($resIvar < 129)
{
	print "before the 2nd while loop\n";
	#going up to 32 procs
	while($procvar < 33)
	{
		#reset walltime var:
		$walltimevar = 10;
		$pbsfilename = $namevar . "_I" . $resIvar . "_J" . $resJvar . "procs_" . $procvar . ".pbs";
		print "$pbsfilename\n";
		$jobnamevar = $namevar . $resIvar . $resJvar . $procvar;
		#want to do whole number multiplication- can't have decimals in the wall time
		$walltimevar = $walltimevar*$countervar; 
		#outputpathvar
		$outputPathvar = "./output/testingsuite/" . $namevar . "_I" . $resIvar . "_J" . $resJvar . "procs_" . $procvar;
		print "$outputPathvar\n";
		#writing the output paths to a text file so we can look em up later to get the cpu time etc etc.
		$outputsfilevar = "outputPathsMG.txt";
		open(outputpathsfile, '>>',$outputsfilevar);
		print outputpathsfile "$outputPathvar ";
		print outputpathsfile "$namevar ";
		print outputpathsfile "$procvar ";
		print outputpathsfile "$resIvar ";
		print outputpathsfile "$resJvar\n";
		close (outputpathsfile);
		#to get mglevels if luke's multigrid works out what the levels are based on how many procs is NOT in the stable build
		#if($procvar == 4){
		#	$resIonprocvar = $resIvar/2;
		#}
		#if($procvar == 8){
		#	$resIonprocvar = $resIvar/4;
		#}
		#if($procvar == 16){
		#	$resIonprocvar = $resIvar/4;
		#}
		#if($procvar == 32){
		#	$resIonprocvar = $resIvar/8;
		#}
		#else{
		#	$resIonprocvar = $resIvar;
		#}
		#if it is in the stable build- works out the levels irrespective of procs
		$resIonprocvar = $resIvar;
		print "$resIonprocvar\n";
		#then count how many times $resIonprocvar can be divided by 2 with no remainders and then you have your # of mg levels
		$counterlevels = 0;
		$rem = 0;
		while($rem == 0){
			$rem = $resIonprocvar%2;
			$resIonprocvar = $resIonprocvar/2;
			$counterlevels = $counterlevels + 1;
		}
		$mgLevelsvar = $counterlevels - 1;
		print "number of MG levels!! $mgLevelsvar\n";
		open(pbsscript, '>',$pbsfilename);
		print pbsscript "#!/bin/bash\n";	
		print pbsscript "#PBS -N " . $filenamevar . $resIvar . $resJvar . $procvar . "MG\n";	
		print pbsscript "#PBS -l nodes=" . $procvar . "\n";
		print pbsscript "#PBS -l walltime=" . $walltimevar . ":00:00" . "\n";	
		print pbsscript "cd " . '$';
		print pbsscript "PBS_O_WORKDIR\n";
		print pbsscript "#PBS -V\n";
		print pbsscript "mpiexec ./Underworld ./InputFiles/" . $filenamevar;
		print pbsscript " --elementResI=" . $resIvar;
		print pbsscript " --elementResJ=" . $resJvar;
		print pbsscript " --outputPath=" . $outputPathvar;
		print pbsscript " --mgLevels=" . $mgLevelsvar;
		print pbsscript " --stopTime=" . $stoptime;
		print pbsscript "\n";
		close (pbsscript);
		#qsub the pbsscript
		$command = "qsub " . $pbsfilename;
		print STDERR "$command\n";
		system($command);
		if($procvar == 1){	
			$procvar = $procvar*4;
		}
		else{
			$procvar = $procvar*2;
		}
	}
	#reset procvar
	$procvar = 1;
	#increment res vars
	$resIvar = $resIvar + $resIvar; 
	$resJvar = $resJvar + $resJvar;
	$counter = $counter + 1;
}

print 'hello this file,';
print filename; 
print 'ran successfully';

