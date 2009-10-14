#! /usr/bin/perl
$outputfile = "outputPathsMGOPTS.txt";
$spreadsheetfile = "AllModelsSpreadSheetMGOPTS2D.txt";
$counter = 0;
@procs = (1,2,4,8,16);

for($i=0;$i<3;$i++)
{
	open(outpathfile,$outputfile);
	foreach $line (<outpathfile>)	
	{
		
		($outputpath, $name, $procvar, $resI, $resJ) = split(' ',$line);
		print $resI;
		#open the FrequentOutput.dat file from the outputpath:
		$outputvar = $outputpath . "/FrequentOutput.dat";
		print "$outputvar $name\n";
		if($counter == 0){
			$resIfirstvar = $resI;
			$rescounter = 0;
			$firstnamevar = $name;
			print "first set of models are $firstnamevar and first resolution is $resIfirstvar\n";
		}
		
		#print "$firstnamevar $name\n";
		if($firstnamevar ne $name){
			$counter = 0;
			print "setting counter back to zero and resetting rescounter and also resIfirst\n";
			$firstnamevar = $name;
			$resIfirstvar = $resI;
			$rescounter = 0;
			
		}
		open(freqout, $outputvar);
		foreach $line (<freqout>)
		{
			#read each line of data and split it into variables
			if($name ne "testVelicSolKzMG"){
				($timestep, $time, $vrms, $CPUtime) = split(' ',$line);
				print "$resI $CPUtime\n";
			}
			else{
				($timestep, $time, $nusselt, $vrms, $CPUtime) = split(' ',$line);
				print "$resI $CPUtime\n";

			}
			
		}
		print "$counter\n";
		#do them as seperate files first then bring them all into 1 file by storing the names of the separate files 
		$filename = "datafromFrequentOutputsForSpreadSheet_" . "$name" . ".txt";
		
				
		open (datafile, '>>', $filename);
		#doing CPU time first
		if($i == 0){
			if($counter == 0){
				print datafile "$name\n";
				print datafile "Res	CPUProc$procs[0]	CPUProc$procs[1]	CPUProc$procs[2]	CPUProc$procs[3]	CPUProc$procs[4]\n";
				print datafile "$resI"."x"."$resJ";
				print datafile "	$CPUtime";
				#only need to store this once
				open(storednamesfile, '>>', 'storednames.txt');
				print storednamesfile "$filename\n";

			}
			if(($resIfirstvar == $resI) && ($firstnamevar eq $name) && ($counter != 0)){
				print datafile "	$CPUtime";
			}
			#do this after as we are setting resIfirst to be resI
			if($resIfirstvar != $resI){
				print datafile "\n$resI"."x"."$resJ";
				
				print datafile "	$CPUtime";
				$resIfirstvar = $resI;
			}
		}
		#then speedup = time taken 1 proc/ time taken x procs
		if($i == 1){
			print $vrms;
			if($counter == 0){
				#setting CPU1_proc
				$CPU_1proc = $CPUtime;
				$speedup = $CPU_1proc/$CPUtime;
				print datafile "\n\n$name\n";
				print datafile "Res	SPEEDUPPROC$procs[0]	SPEEDUPPROC$procs[1]	SPEEDUPPROC$procs[2]	SPEEDUPPROC$procs[3]	SPEEDUPPROC$procs[4]\n";
				print datafile "$resI"."x"."$resJ";
				print datafile "	$speedup";
				
			}
			if(($resIfirstvar == $resI) && ($firstnamevar eq $name) && ($counter != 0)){
				$speedup = $CPU_1proc/$CPUtime;
				print datafile "	$speedup";
			}
			#do this after as we are setting resIfirst to be resI and we are starting from 1 proc again so we need to set CPU_1proc
			if($resIfirstvar != $resI){
				$CPU_1proc = $CPUtime;
				$speedup = $CPU_1proc/$CPUtime;
				print datafile "\n$resI"."x"."$resJ";
				print datafile "	$speedup";
				$resIfirstvar = $resI;
			}
		}
		#then efficiency = speedup / proc
		if($i == 2){
			print $time;
			
			if($counter == 0){
				#need to set CPU1_proc
				$CPU1_proc = $CPUtime;
				$efficiency = ($CPU1_proc/$CPUtime)/$procvar;
				print datafile "\n\n$name\n";
				print datafile "Res	EFFICIENCYProc$procs[0]	EFFICIENCYProc$procs[1]	EFFICIENCYProc$procs[2]	EFFICIENCYProc$procs[3]	EFFICIENCYProc$procs[4]\n";
				print datafile "$resI"."x"."$resJ";
				print datafile "	$efficiency";
				
			}
			if(($resIfirstvar == $resI) && ($firstnamevar eq $name) && ($counter != 0)){
				
				$efficiency = ($CPU1_proc/$CPUtime)/$procvar;
				print datafile "	$efficiency";
			}
			#do this after as we are setting resIfirst to be resI
			if($resIfirstvar != $resI){
				#need to set CPU1_proc as we are now going on to another resolution:
				$CPU1_proc = $CPUtime;
				$efficiency = ($CPU1_proc/$CPUtime)/$procvar;
				print datafile "\n$resI"."x"."$resJ";
				
				print datafile "	$efficiency";
				$resIfirstvar = $resI;
			}
		}
		$counter++;
		print "i is $i\n";
		if(i == 3){
			print datafile "\n";
		}
		close (datafile);
	}
	close(outpathfile);
}
close(storednamesfile);
#now put everything in the one file:
open(completespreadsheet,'>>',$spreadsheetfile);
open(storednamesfile,'storednames.txt');
foreach $mline (<storednamesfile>)
{
	$mline =~ s/\s+$//;
	print $mline;
	$modelfilename = $mline;
	open(datafile,$modelfilename);
	foreach $line (<datafile>)
	{
		#print "$line";
		print completespreadsheet "$line";
	}
	print completespreadsheet "\n";
	print completespreadsheet "\n";
	print completespreadsheet "\n";
	print completespreadsheet "\n";
	close(datafile);
}
close(spreadsheetfile);


print 'hello this script,';
print 'ran successfully';

