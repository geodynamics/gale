#!/usr/bin/perl -w
#
use strict;

sub runTests;
sub testConvergence;
sub generateConvergence;
sub testNumbersAgainstExpected;
sub readOptionsFile;
our $progHelp = "To run convergence tests:

./runAndTestConvergence.pl <xmlFile> [ OPTIONS ]

where OPTIONS:
  -optionsFile <fileName>   : where <fileName> is the option file. By default serial runs of resolution 16sq, 32sq, 64sq, 124sq are done.
  -againstExpected          : will check the generated numbers against a file in the 'expected' dir. By default check will not occur

EXAMPLE:
  ./runAndTestConvergence.pl testVelicSolS.xml -optionsFile OFile.dat
      Runs with option file OFile.dat and NO check against the expected file.

These scripts measure the convergence rates in numerical fields generated in Underworld against analytic solutions supplied by Mirko Velic.
The error of FEM fields should decrease with higher resolution - depending on the choice of element type that represents a given field.
These scripts run error measures on given field at several resolutions and then process how the error converges.

The optionsFile can be used to set how each run of the analytic problem differs from the next, see the OFile.dat for details.\n\n\n";

######  Main program    ######

# 1) Run the same xml at the different resolutions 
#    these will generate error data
&runTests();

# 2) Use a linear regression algorithm on the data generated
#    and check the results for good behavior and acceptable convergence
&testConvergence( &generateConvergence() );

######  End main program ######

sub runTests {
	my $res;
	my $command;

# read commandline args
	my $arg;
	my $ii = 0;
	my $xmlFile = " ";
	my $optFile = " ";
	my $testsToRun = 4;
	my @procs = (1,1,1,1);
	my @commandLines = ("--elementResI=16 --elementResJ=16 --elementResK=16",
											"--elementResI=32 --elementResJ=32 --elementResK=32",
											"--elementResI=64 --elementResJ=64 --elementResK=64",
											"--elementResI=128 --elementResJ=128 --elementResK=128");
												
	# check if xml exists and options file is specified
	foreach $arg (@ARGV) {
		if( $arg =~ m/.*\.xml$/ ) { $xmlFile = $arg; }
		elsif( $arg =~ m/\-optionsFile/ ) { $optFile = $ARGV[$ii+1]; $ii++; }
		elsif( $arg =~ m/^\-h$/ ) { print $progHelp; exit }
		elsif( $arg =~ m/^\-\-help$/ ) { print $progHelp; exit }
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
	# Need to check for an executable
	if( (-e "../../../build/bin/StGermain" eq 0 ) ) {
		die "\nCan't find ../../../build/bin/StGermain - the executable which runs the test, stopped";
	}
	print "\n--- Testing the convergence of $xmlFile ---\n";
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


# commence running 
	my $run_I = 0;

	for( $run_I = 0 ; $run_I < $testsToRun ; $run_I++ ) {
		# run test case
		$command = "mpirun -np $procs[$run_I] ./$exec $xmlFile $commandLines[$run_I] --pluginData.appendToAnalysisFile=True >$stdout";
		print "$command";
		$command .= " 2>$stderr";
		`$command`;

			# check error stream for error result
      open( ERROR, "<$stderr" );
			my $line;
			foreach $line (<ERROR>) {
				if( $line =~ m/[E|e]rror/ ) {
					close(ERROR); 
#	$command = "rm $stderr"; `$command`;
				 	die ("\n\tError: see $stderr or $stdout - stopped" ); 
				}
			}
			close(ERROR); 
			$command = "rm $stderr"; `$command`;

		print " ....finished\n\n";
	}
	# removing softlink
	$command = "rm $exec";
	print "$command\n";
	`$command`;
	print "--- Finished convergence runs ---\n\n";
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

sub testNumbersAgainstExpected {
	my $datFile = $_[0];
	my $expectedFile = $_[1];
	my $line;
	my $ii = 0;
	my $needLabels = 1;
	my $jj = 0;
	my @input;
	my @keys;
	my @expected;
	
	open( INPUT, "<$datFile");
	print "\n--- Testing the output against the expected file $expectedFile ---\n\n";

	# get results from inputfile
	foreach $line (<INPUT>) {
		if( $line =~ m/^\-\-\-\-\-\-\-\-\-/ ) {
    # exit input if ------ is detected
			last;
		} 
		if( ($needLabels == 1) && $line=~m/^#Res\s+(.*)/ ) {
			@keys = split(/\s+/, $1 );
			$needLabels = 0;
		}
		if( $line =~ m/^\s*\D/ ) {
		# disregard lines starting with non-Digits
			next;
		} 
		if( $line =~ m/^\s*\d/ ) {
		 	$input[$ii] = [ split (/\s/, $line) ];
			$ii++;
		}

	}
	close(INPUT);

	# get results from expected file
	$ii=0;
	open( EXPECTED, "<$expectedFile");
	foreach $line (<EXPECTED>) {
		if( $line =~ m/^\-\-\-\-\-\-\-\-\-/ ) {
    # exit input if ------ is detected
			last;
		} 
		elsif( $line =~ m/^\s*\D/ ) {
		# disregard lines starting with non-Digits
			next;
		} 
		elsif( $line =~ m/^\s*\d/ ) {
		 	$expected[$ii] = [ split (/\s/, $line) ];
			$ii++;
		}
	}
	close(EXPECTED);

	if( scalar(@input) != scalar(@expected) ) { 
		die "Error: The expected file \n\"$expectedFile\"\n".
			"has a different number of results than the test file \n\"$datFile\"\n".
			"so the results testing can't be done.\n".
			"Remove argument '-againstExpected' to not check the resulting errors agains the expected file.\n"; 
	}

	my $error;
	my $r_error;
	my $r_tol = 0.05;
	print "relative tolerance is set to $r_tol%\n";
	for( $ii = 0 ; $ii < scalar(@input) ; $ii++ ) {
		for( $jj = 1 ; $jj < scalar(@keys) ; $jj++ ) {
			$error = abs( $input[$ii][$jj] - $expected[$ii][$jj] );
			$r_error = 100 * ($error / $expected[$ii][$jj]);
			if( $r_error > $r_tol ) {
				die "Results of $keys[$jj-1] differs by $r_error%\n";
			} 
			# This is to verbose be could be useful to understand error
			#else {
				# print "Results of $keys[$jj-1] are $r_error\n";
			# }
		}
	}
	print "All values within tolerance\n";
	print "\n-- Finished testing the output against the expected file ---\n";
}

sub generateConvergence {
	# check if the cvg can be found
	my $command = "ls *\.cvg 2>/dev/null";
	my $datFile;
	my $cvgFile = `$command`;
	chomp $cvgFile;
	if( $cvgFile eq "" ) { die "There is no cvg file found analyse\n" ; }

	# run David A. May's convergence tester, cheers Mayhem
	$command = "./generate-convergence-rates.pl -errorfile $cvgFile";
	`$command`;

	$datFile = "$cvgFile\-convergence\-rates.dat";
	if( !(-e $datFile) ) { die "Can't find the convergence rate file: $datFile\n"; }

	# search ARGV for againstExpected flag
	my $arg;
	my $runAgainstExpected = 0;
	foreach $arg (@ARGV) {
		if( $arg =~ m/\-againstExpected/ ) {
			$runAgainstExpected = 1;
			last;
		}
	}


	if( $runAgainstExpected ) {
	if( -e "./expected/$datFile" ) { 
			&testNumbersAgainstExpected($datFile, "./expected/$datFile");
		}	else { 
			die "Can't find the expected file: ./expected/$datFile\n". 
			"Consider adding argument '-convergeOnly' to the executable to not check the resulting errors agains the expected file.\n"; 
		}
	}
	
	# remove the old file
	$command = "rm $cvgFile";
	`$command`;
	return $datFile;
}

sub testConvergence {
 	my $datFile = $_[0]; 
	my $command;
	my @keys;
	my @cvgRates;
	my @correlations;
	my $needLabels = 1;
	my $line;
	 #test convergence numbers
	open(INPUT, "<$datFile") || die "Can't open the expected file $datFile\n" ; 
	while ($line = <INPUT>) {
		chomp $line;
    # parse for convergence rates
		if ( $line =~ m/^cvg\.\srate\s+(.*)/ ) {
			@cvgRates = split (/\s+/, $1 );
		}
    # parse for correlation results
		elsif( $line =~ m/^corr\.\s+(.*)/ ) {
			@correlations = split (/\s+/, $1 );
		}
    # parse for variable labels
		elsif( ($needLabels == 1) && $line=~m/^#Res\s+(.*)/ ) {
			@keys = split(/\s+/, $1 );
			$needLabels = 0;
		}
	}

  # now check results and report
	my $ii;
	my $jj;
	my $length;
	my $status;
	my $result = "Pass";
	my $report;
	my $padding;
	print "--- Results of convergence ---\n\n";
	for( $ii = 0 ; $ii < scalar(@correlations) ; $ii++ ) {
		$status = "good";
		$report = "";
		$padding = "";

		$length = 25-length($keys[$ii]); 
		for( $jj = 0 ; $jj < $length ; $jj++ ) {
			$padding .= " ";
		}

		print "$keys[$ii] $padding- convergence rate = $cvgRates[$ii] - correlation $correlations[$ii] .... ";
		if( $correlations[$ii] < 0.99 ) {
			$status = "BAD!!!";
			$result = "Fail";
			$report = $report."\thas a correlation of $correlations[$ii], below the accepted threshold of 0.99\n";
		}
		if( $keys[$ii] =~ m/^Velocity/ ) {
			if( $cvgRates[$ii] < 1.6 ) {
				$status = "BAD!!!";
				$result = "Fail";
				$report = $report."\thas a convergence rate of $cvgRates[$ii], theoretically velocity fields should converge at a rate ~ 2\n";
			}
		} elsif ($keys[$ii] =~ m/Pressure/ ) {
			if( $cvgRates[$ii] < 0.9 ) {
				$status = "BAD!!!";
				$result = "Fail";
				$report = $report."\thas a convergence rate of $cvgRates[$ii], theoretically pressure fields should converge at a rate ~ 1\n";
			}
		} elsif ($keys[$ii] =~ m/StrainRate/ ) {
				if( $cvgRates[$ii] < 0.85 ) {
					$status = "BAD!!!";
					$result = "Fail";
					$report = $report."\thas a convergence rate of $cvgRates[$ii], theoretically pressure fields should converge at a rate ~ 1\n";
				}
		}
		print "$status\n$report";
	}

	close( INPUT );
	# move cvg dat to the log directory
	$command = "mv $datFile log/";
	`$command`;

	print "\nResult = $result\n";
}
