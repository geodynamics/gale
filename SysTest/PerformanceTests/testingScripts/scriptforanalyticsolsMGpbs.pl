#! /usr/bin/perl
@namevar = ("testVelicSolAMG", "testVelicSolCxMG", "testVelicSolKzMG");
@mgLevelsvar = ("4", "5", "6");
$procvar = 1;
$walltimevar = 2;
$resIvar = 32;
$resJvar = 32;
$countervar = 1;
print "doing the for loop\n";
#going up to res 128x128
for($i = 0; $i < 3; $i++)
{
	print "before the while loop\n";
	while($resIvar < 129){
		while($procvar < 17){
			$filenamevar = "$namevar[$i].xml";
			$mg = $mgLevelsvar[$i];
			$pbsfilename = $namevar[$i] . "_I" . $resIvar . "_J" . $resJvar . "_K".$resKvar."procs_" . $procvar . ".pbs";
			print "$pbsfilename\n";
			$jobnamevar = $namevar[$i] . $resIvar . $resJvar . $resKvar . $procvar;
			#outputpathvar	
			$outputPathvar = "./output/MG/" . $namevar[$i] . "_I" . $resIvar . "_J" . $resJvar . "_P" . $procvar;
			print "$outputPathvar\n";
			#writing the output paths to a text file so we can look em up later to get the cpu time etc etc.
			$outputsfilevar = "outputPathsMG.txt";
			open(outputpathsfile, '>>',$outputsfilevar);
			print outputpathsfile "$outputPathvar ";
			print outputpathsfile "$namevar[$i] ";
		  	print outputpathsfile "$procvar ";	
			print outputpathsfile "$resIvar ";
			print outputpathsfile "$resJvar\n";
			close (outputpathsfile);
			open(pbsscript, '>',$pbsfilename);
			print pbsscript "#!/bin/bash\n";	
			print pbsscript "#PBS -N " . $namevar[$i] . $resIvar . $resJvar . $procvar . "\n";
			print pbsscript "#PBS -l ncpus=" . $procvar . "\n";
			print pbsscript "#PBS -l walltime=" . $walltimevar . ":00:00\n";
			print pbsscript "cd " . '$';
			print pbsscript "PBS_O_WORKDIR\n";
			print pbsscript "export MPI_DIR=" . '$' . "PBS_O_WORKDIR\n";
			print pbsscript "#PBS -V\n";
			print pbsscript "mpiexec ./Underworld ./InputFiles/$filenamevar";
			print pbsscript " --elementResI=$resIvar";
			print pbsscript " --elementResJ=$resJvar";
			print pbsscript " --mgLevels=$mg";
			print pbsscript " --outputPath=$outputPathvar\n";
			close (pbsscript);
			
			#qsub the pbsscript	
			$command = "qsub $pbsfilename";
			print STDERR "$command\n";
			system($command);	
			$procvar = $procvar + $procvar;		
		}
		#increment res vars
		$resIvar = $resIvar + $resIvar; 
		$resJvar = $resJvar + $resJvar;
		#reset procvar
		$procvar = 1;
		$countervar = $countervar + 1;
	}
	#reset resvars
	$resIvar = 32;
	$resJvar = 32;
}

print 'hello this file,';
print filename; 
print 'ran successfully';

