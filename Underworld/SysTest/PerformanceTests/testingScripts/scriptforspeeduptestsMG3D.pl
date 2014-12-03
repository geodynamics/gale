#! /usr/bin/perl
$namevar = 'RayleighTaylorBenchmarkMG';
$filenamevar = $namevar . '.xml';
$procvar = 1;
$walltimevar = 10;
$resIvar = 8;
$resJvar = 8;
$resKvar = 8;
$mgLevelsvar = 1;
$stoptime = 1000;
$countervar = 1;
print "before the first while loop\n";
#going up to res 32x32x32
while($resIvar < 33)
{
	print "before the 2nd while loop\n";
	#going up to 32 procs
	while($procvar < 33)
	{
		#reset walltime var:
		$walltimevar = 10;
		$pbsfilename = $namevar . "_I" . $resIvar . "_J" . $resJvar . "_K".$resKvar."procs_" . $procvar . ".pbs";
		print "$pbsfilename\n";
		$jobnamevar = $namevar . $resIvar . $resJvar . $resKvar . $procvar;
		#want to do whole number multiplication- can't have decimals in the wall time
		$walltimevar = $walltimevar*$countervar; 
		#outputpathvar
		$outputPathvar = "./output/testingsuite/" . $namevar . "_I" . $resIvar . "_J" . $resJvar . "_K".$resKvar."procs_" . $procvar;
		print "$outputPathvar\n";
		#writing the output paths to a text file so we can look em up later to get the cpu time etc etc.
		$outputsfilevar = "outputPathsMG3D.txt";
		open(outputpathsfile, '>>',$outputsfilevar);
		print outputpathsfile "$outputPathvar ";
		print outputpathsfile "$namevar ";
		print outputpathsfile "$procvar ";
		print outputpathsfile "$resIvar ";
		print outputpathsfile "$resKvar ";
		print outputpathsfile "$resKvar\n";
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
		$resJonprocvar = $resJvar;
		print "$resJonprocvar\n";
		#then count how many times $resIonprocvar can be divided by 2 with no remainders and then you have your # of mg levels
		$counterlevels = 0;
		$rem = 0;
		while($rem == 0){
			$rem = $resJonprocvar%2;
			$resJonprocvar = $resJonprocvar/2;
			$counterlevels = $counterlevels + 1;
		}
		$mgLevelsvar = $counterlevels - 1;
		print "number of MG levels!! $mgLevelsvar\n";
		open(pbsscript, '>',$pbsfilename);
		print pbsscript "#!/bin/bash\n";	
		print pbsscript "#PBS -N " . $filenamevar . $resIvar . $resJvar . $resKvar . $procvar . "MG\n";	
		print pbsscript "#PBS -l nodes=" . $procvar . "\n";
		print pbsscript "#PBS -l walltime=" . $walltimevar . ":00:00" . "\n";	
		print pbsscript "cd " . '$';
		print pbsscript "PBS_O_WORKDIR\n";
		print pbsscript "#PBS -V\n";
		print pbsscript "mpiexec ./Underworld ./InputFiles/" . $filenamevar;
		print pbsscript " --elementResI=" . $resIvar;
		print pbsscript " --elementResJ=" . $resJvar;
		print pbsscript " --elementResK=" . $resKvar;
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
	$resKvar = $resKvar + $resKvar;
	$counter = $counter + 1;
}

print 'hello this file,';
print filename; 
print 'ran successfully';

