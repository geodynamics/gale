#! /usr/bin/perl
$outputfile = "outputPaths.txt";
$counter = 0;
open(outpathfile,$outputfile);
foreach $line (<outpathfile>) 
{
	
	($outputpath, $name, $procvar, $resI, $resJ) = split(' ',$line);
	#open the FrequentOutput.dat file from the outputpath:
	$outputvar = $outputpath . "/FrequentOutput.dat";
	print "$outputvar $name\n";
	if($counter == 0){
		$resIfirstvar = $resI;
		$firstnamevar = $name;
		print "first set of models are $firstnamevar and first resolution is $resIfirstvar\n";
	}
	
	#print "$firstnamevar $name\n";
	if($firstnamevar ne $name){
		$counter = 0;
		print "setting counter back to zero\n";
	}
	open(freqout, $outputvar);
	#read the first line
	($char, $itimestep, $itime, $ivrms, $iCPUtime) = split(' ',$line);
	print "$char, $itimestep, $itime, $inusselt, $ivrms, $iCPUtime\n";
		
	foreach $line (<freqout>)
	{
		#read each line of data and split it into variables
		($timestep, $time, $vrms, $CPUtime) = split(' ',$line);
		print "$timestep, $time, $nusselt, $vrms, $CPUtime\n";
		print "$resI $CPUtime\n";
		
	}
	print "$counter\n";
	$filename = 'datafrom' . $name . 'FrequentOutputs.dat';
	open (datafile, '>>', $filename);
	
	if($counter == 0){
		#store the name of the file in a separate file to access later to put into a spreadsheet
		$fileslist = 'alldatafilesforspreadsheet.dat';
		open (listoffiles, '>>', $fileslist);
		print listoffiles "$filename\n";
		print datafile "Res	CPUProc1	CPUProc4	CPUProc8	CPUProc16	CPUProc32\n";
		print datafile "$resI"."x"."$resJ";
		print datafile "	$CPUtime";
		print "$CPUtime\n";
	}
	if(($resIfirstvar == $resI) && ($firstnamevar eq $name) && ($counter != 0)){
		print datafile "	$CPUtime";
	}
	#do this after as we are resetting resIfirstvar
	if($resIfirstvar != $resI){
                print datafile "\n$resI"."x"."$resJ";
                print datafile "	$CPUtime";
                $resIfirstvar = $resI;
        }

	$counter++;
	close (datafile);
	
}

print 'hello this script,';
print 'ran successfully';

