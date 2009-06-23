#!/usr/bin/perl -w
#
use strict;

sub runTests;
sub testConvergence;
sub generateConvergence;
sub testNumbersAgainstExpected;
sub readOptionsFile;
our $progHelp = "To run checkpoint tests:

./checkpointTest.pl <xmlFile> [ OPTIONS ]

where OPTIONS:
	-optionsFile <fileName>   : where <fileName> is the options file. Command line agruments for the run are set in the options file
	-c                        : will \"create\" checkpointed data only. By default this flag in not set and the script only checks against previous checkpointed data.  
	-n                        : the timestep checkpoint writing (if -c is defined) or checkpoint testing will occur on. By default this is timestep 10.  
	-h                        : help

EXAMPLE:
  ./checkpointTest.pl testVelicSolS.xml -optionsFile OFile.dat\n\t(Runs with option file OFile.dat and checks against the expected file)\n\n";

######  Main program    ######

# 1) Run the xml 

# 2) Check against expected, checkpoint files 
exit &testConvergence( &runTests() );

######  End main program ######

sub runTests {
	my $res;
	my $command;
	my $createTest=0; # create an expected file

# read commandline args
	my $arg;
	my $ii = 0;
	my $xmlFile = " ";
	my $optFile = " ";
	my $testsToRun = 1;
	my $numberOfTimeSteps = 10; # testing Timestep is 10 by default
	my @procs = (1,1,1,1);
	my @commandLines = ("--elementResI=32 --elementResJ=32 " );
												
	# check if xml exists and options file is specified
	foreach $arg (@ARGV) {
		if( $arg =~ m/.*\.xml$/ ) { $xmlFile = $arg; }
		elsif( $arg =~ m/\-optionsFile/ ) { $optFile = $ARGV[$ii+1]; $ii++; }
		elsif( $arg =~ m/^\-h$/ ) { print $progHelp; exit }
		elsif( $arg =~ m/^\-\-help$/ ) { print $progHelp; exit }
		elsif( $arg =~ m/^\-c/ ) { $createTest=1; }
		elsif( $arg =~ m/^\-n/ ) { $numberOfTimeSteps=$ARGV[$ii+1]; $ii++; }
		$ii++;
	}
	if( $xmlFile eq " " ) { die "\nNo xml file specified, stopped" ; }
	if( !(-e $xmlFile) ) { die "\nCannot find input file: $xmlFile, stopped" ; }

	# check if options file is given, otherwise run default
	if( $optFile eq " " ) {
		warn "\nNo run options file specifed using commandline arg '-optionsFile xxx'".
				 "\nUsing default serial mode, with commandLines:\n"; foreach (@commandLines) { print "$_\n"; }
	}	else {
		if( !(-e $optFile) ) { die "\nCannot find run options file $optFile, stopped"; }
		# read in run options file
		$testsToRun = &readOptionsFile( $optFile, \@procs, \@commandLines );
	}


# do checks on executables and log files

	my $exec = "udw";
	my $stdout;
	my $stderr;

	my $xmlSegmentCreateTest = "<StGermainData xmlns=\"http://www.vpac.org/StGermain/XML_IO_Handler/Jun2003\">
		<param name=\"checkpointEvery\" mergeType=\"replace\">10</param>
		<param name=\"outputPath\" mergeType=\"replace\">./expected/$xmlFile</param>
		<param name=\"maxTimeSteps\" mergeType=\"replace\">$numberOfTimeSteps</param>
		<param name=\"dumpEvery\" mergeType=\"replace\">0</param>
</StGermainData>";
	my $xmlSegmentToTest = "<StGermainData xmlns=\"http://www.vpac.org/StGermain/XML_IO_Handler/Jun2003\">
	<struct name=\"components\" mergeType=\"merge\">
		<struct name=\"tester\">
			<param name=\"Type\">FieldTest</param>
		</struct>	
	</struct>

	<param name=\"outputPath\" mergeType=\"replace\">./output/$xmlFile</param>
	<param name=\"checkpointEvery\" mergeType=\"replace\">0</param>
	<param name=\"dumpEvery\" mergeType=\"replace\">0</param>
	<param name=\"maxTimeSteps\" mergeType=\"replace\">$numberOfTimeSteps</param>
	<struct name=\"pluginData\" mergeType=\"replace\">
		<list name=\"NumericFields\">
			<param>VelocityField</param> <param>0</param>
			<param>PressureField</param> <param>1</param>
		</list> 
		<param name=\"IntegrationSwarm\">gaussSwarm</param>
		<param name=\"ConstantMesh\">constantMesh</param>
		<param name=\"testTimestep\">$numberOfTimeSteps</param>
		<param name=\"ElementMesh\">linearMesh</param>
		<param name=\"normaliseByAnalyticSolution\">True</param>
		<param name=\"context\">context</param>
		<param name=\"appendToAnalysisFile\">True</param>
		<!-- reference soln stuff -->
		<param name=\"useReferenceSolutionFromFile\">true</param>
		<param name=\"referenceSolutionFilePath\">./expected/$xmlFile</param>
		<list name=\"ReferenceFields\">
			<param>VelocityField</param>
			<param>PressureField</param>
		</list> 
	</struct> 
</StGermainData>";

	# Need to check for an executable
	if( (-e "../../../build/bin/StGermain" eq 0 ) ) {
		die "\nCan't find ../../../build/bin/StGermain - the executable which runs the test, stopped";
	}
	print "\n--- Testing the $xmlFile ---\n";
	# is the symbolic link there, if not create it
	if( !(-e $exec) ) {
		$command = "ln -s ../../../build/bin/StGermain $exec";
		print "\n$command\n\n";
		`$command`;
	}	
	# check if there's a log dir
	if( !(-e "log/") ) {
		$command = "mkdir log";
		`$command`;
	}
	# declare stdout and stderr files, in log dir.
	$stdout = "log/$xmlFile"."_runs.stdout";
	$stderr = "log/$xmlFile"."_runs.stderr";

	# remove old log file, if it exists
	if( -e "$stdout" ) {
		$command = "rm $stdout";
		`$command`;
	}
	# remove old cvg file, if it exists
	$command = "ls *\.cvg 2>/dev/null";
	my $cvg = `$command`;
	chomp( $cvg );
	if( $cvg ne "" ) { $command = "rm $cvg"; `$command`;} 

	# create help.xml for setting up test
	if( $createTest ) {
		$command = "echo \'$xmlSegmentCreateTest\' > help.xml ";
	} else {
		$command = "echo \'$xmlSegmentToTest\' > help.xml ";
	}
	`$command`;

 # commence running 
	my $run_I = 0;

	for( $run_I = 0 ; $run_I < $testsToRun ; $run_I++ ) {
		# run test case
		$command = "mpirun -np $procs[$run_I] ./$exec $xmlFile help.xml $commandLines[$run_I] --pluginData.appendToAnalysisFile=True >$stdout";
		print "$command";
		$command .= " 2>$stderr";
		`$command`;

			# check error stream for error result
      open( ERROR, "<$stderr" );
			my $line;
			foreach $line (<ERROR>) {
				if( $line =~ m/[E|e]rror/ ) {
					close(ERROR); 
				 	die ("\n\tError: see $stderr or $stdout - stopped" ); 
				}
			}
			# if no error close file and delete it
			close(ERROR); 
			$command = "rm $stderr"; `$command`;

		print " ....finished\n\n";
	}

	# removing help.xml
	$command = "rm help.xml";
	print "$command\n";
  `$command`;

	# removing softlink
	$command = "rm $exec";
	print "$command\n";
	`$command`;
	print "--- Finished  ---\n\n";

	if( $createTest ) { exit(0); }

	$command = "ls *\.cvg 2>/dev/null";
	$cvg = `$command`;
	chomp( $cvg );
	return $cvg;
}

sub readOptionsFile {
	my ( $optFile, $procs, $commandLines ) = @_;
	my $line;
	# $line_I represents the number of tests to run
	my $line_I = 0;
	# open options file
	open OPTFILE, "<$optFile" || die "Can't open options file $optFile, stopped" ;
	foreach $line (<OPTFILE>) {
		chomp $line;
		# only process lines that start with np
		if( $line =~ m/^np\s+(\d+)\s+(.*)/ ) {
			$procs->[$line_I] = $1;
			$commandLines->[$line_I] = $2;
			$line_I++;
		} else { next; }
	}
	return $line_I;
}

sub testConvergence {
 	my $datFile = $_[0]; 
	my @keys;
	my $tolerance = 0.000001;
	my @errors;
	my $line;
	my $nKeys;
	my $nErrs;
	my $report;
	my $result;
	my $command;
	my $ii;
	# test convergence numbers
	open(INPUT, "<$datFile") || die "Can't open the expected file $datFile\n" ; 
	while ($line = <INPUT>) {
		chomp $line;
		if ( $line =~ m/^\#Res\s.*/ ) {
   		# parse for variable labels
			@keys = split (/\s+/, $line );
		}
		else {
			@errors = split(/\s+/, $line );
		}
	}

	# ensure the number of keys and error measures agree
	$nKeys = @keys;
	$nErrs = @errors;
	
	if( $nKeys != $nErrs ) { die "The number of keys against the number of errors in file $datFile don't agreed\n"; }

	$result = "Pass";
	$report = "";

	# go through all  errors and check if they're within tolerance
	for( $ii = 1 ; $ii < $nKeys ; $ii++ ) {
		if( abs($errors[$ii]) > 0.02 ) {
			$result = "Fail";
			$report .= "***BAD NEWS*** ... $keys[$ii] differs by more than " . $tolerance*100 . "\% tolerance from expected file\n";
		} else {
			$report .= "pass ... $keys[$ii] within a ". $tolerance*100 ."\% relative tolerance from expected file\n";
		}
	}

	close( INPUT );

	print "\n$report";
	print "\nResult = $result\n";

	# remove the used data file
	$command = "rm $datFile";
	`$command`;

	if( $result eq "Pass" ) {
		exit(0);
	} else {
		exit(1);
	}
}
