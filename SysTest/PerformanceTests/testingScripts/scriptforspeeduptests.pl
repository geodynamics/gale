#! /usr/bin/perl
$namevar = 'ThermalConvection';
$filenamevar = $namevar . '.xml';
$procvar = 1;
$walltimevar = 10;
$resIvar = 32;
$resJvar = 32;
$stoptime = 5000;
$countervar = 1;
print "before the first while loop\n";
#going up to res 128x128
while($resIvar < 129)
{
	print "before the 2nd while loop\n";
	#going up to 32 procs
	while($procvar <= 33)
	{
		#reset walltimevar
		$walltimevar = 10;
		$pbsfilename = $namevar . "_I" . $resIvar . "_J" . $resJvar . "procs_" . $procvar . ".pbs";
		print "$pbsfilename\n";
		$jobnamevar = $namevar . $resIvar . $resJvar . $procvar;
		#want to do whole number division because we can't have decimals in the wall time
		$walltimevar = $walltimevar*$countervar; 
		print "walltime is $walltimevar procs are $procvar counter is $countervar";
		#outputpathvar
		$outputPathvar = "./output/testingsuite/" . $namevar . "_I" . $resIvar . "_J" . $resJvar . "procs_" . $procvar;
		print "$outputPathvar\n";
		#writing the output paths to a text file so we can look em up later to get the cpu time etc etc.
		$outputsfilevar = "outputPaths.txt";
		open(outputpathsfile, '>>',$outputsfilevar);
		print outputpathsfile "$outputPathvar ";
		print outputpathsfile "$namevar ";
		print outputpathsfile "$procvar ";
		print outputpathsfile "$resIvar ";
		print outputpathsfile "$resJvar\n";
		close (outputpathsfile);

		open(pbsscript, '>',$pbsfilename);
		print pbsscript "#!/bin/bash\n";	
		print pbsscript "#PBS -N " . $filenamevar . $resIvar . $resJvar . $procvar . "\n";
		print pbsscript "#PBS -l nodes=" . $procvar . "\n";
		print pbsscript "#PBS -l walltime=" . $walltimevar . ":00:00" . "\n";	
		print pbsscript "cd " . '$';
		print pbsscript "PBS_O_WORKDIR\n";
		print pbsscript "#PBS -V\n";
		print pbsscript "mpiexec ./Underworld ./InputFiles/$filenamevar";
		print pbsscript " --elementResI=$resIvar";
		print pbsscript " --elementResJ=$resJvar";
		print pbsscript " --outputPath=$outputPathvar";
		print pbsscript " --stopTime=$stoptime\n";
		close (pbsscript);
		#qsub the pbsscript
		$command = "qsub $pbsfilename";
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
	$countervar = $countervar + 1;
}

print 'hello this file,';
print filename; 
print 'ran successfully';

