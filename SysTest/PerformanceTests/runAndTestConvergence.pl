#!/usr/bin/perl -w
#
use strict;

sub runTests;
sub testConvergence;
sub generateConvergence;
sub testNumbersAgainstExpected;

######  Main program    ######

# resolutions to run convergence tests at
my @resolutions = ( 16, 32, 64 ); #, 128, 256 );

# 1) Run the same xml at the resolutions specified in '@resolutions'
#    these will generate error data
&runTests();

# 2) Use a linear regression algorithm on the data generated
#    and check the results for good behavior and acceptable convergence
#&generateConvergence();
&testConvergence( &generateConvergence() );

######  End main program ######

sub runTests {
	my $arg;
	my $res;
	my $command;
	my $exec = "udw";
	my $stdout;
	my $stderr;
	my $xmlFile = "NULL";

# check if xml exists
	foreach $arg (@ARGV) {
		if( $arg =~ m/.*\.xml$/ ) { $xmlFile = $arg; }
	}

	if( $xmlFile eq "NULL" ) { die "\nNo xml file specified, stopped" ; }
	if( !(-e $xmlFile) ) { die "\nCannot find input file: $xmlFile, stopped" ; }

# declare stdout and stderr files, in log dir.
	$stdout = "log/$xmlFile"."_runs.stdout";
	$stderr = "log/$xmlFile"."_runs.stderr";

	print "\n### Testing the convergence of $xmlFile with resolutions:";
	foreach $res (@resolutions) {
		print " $res";
	}
	print " ###\n\n";

# Need to check for an executable
	if( (-e "../../../build/bin/StGermain" eq 0 ) ) {
		die "Can't find ../../../build/bin/StGermain - the executable which runs the test, stopped";
	}
	print "Making softlink to executable ../../../build/bin/StGermain called $exec\n\n";
	$command = "ln -s ../../../build/bin/StGermain $exec";
	`$command`;


# remove old log file if it exists
	if( (-e "$stdout" ne 0 ) ) {
		$command = "rm $stdout";
		`$command`;
	}

	foreach $res (@resolutions) {
		# run test case
		$command = "./$exec $xmlFile --elementResI=$res --elementResJ=$res --pluginData.appendToAnalysisFile=True >>$stdout";
		print "Runing $command";
		$command .= " 2>>$stderr";
		`$command`;

			# check error stream for error result
			open( ERROR, "<$stderr" );
			my $line;
			foreach $line (<ERROR>) {
				if( $line =~ m/[E|e]rror/ ) {
					close(ERROR); 
					$command = "rm $stderr"; `$command`;
				 	die ("\n\tError: see $stderr - stopped" ); 
				}
			}
			close(ERROR); 
			$command = "rm $stderr"; `$command`;

		print "\n\t\t\t\t....finished\n\n";
	}
# removing softlink
	print "Removing the softlink $exec\n";
	print "### Finished convergence runs ###\n\n";
	$command = "rm $exec";
	`$command`;
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
	print "\n### Testing the output against the expected file $expectedFile ###\n\n";

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
			"Add argument '-convergeOnly' to the executable to not check the resulting errors agains the expected file.\n"; 
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
	print "\n### Finished testing the output against the expected file ###\n";
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

	# search ARGV for convergeOnly flag
	my $arg;
	my $runAgainstExpected = 1;
	foreach $arg (@ARGV) {
		if( $arg =~ m/\-convergeOnly/ ) {
			$runAgainstExpected = 0;
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
	print "\nResults:\n";
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
			$report = $report . "\thas a correlation of $correlations[$ii], below the accepted threshold of 0.99\n";
		}
		if( $keys[$ii] =~ m/^Velocity/ ) {
			if( $cvgRates[$ii] < 1.6 ) {
				$status = "BAD!!!";
				$result = "Fail";
				$report = $report . "\thas a convergence rate of $cvgRates[$ii], theoretically velocity fields should converge at a rate ~ 2\n"
			}
		} elsif ($keys[$ii] =~ m/Pressure/ ) {
			if( $cvgRates[$ii] < 0.9 ) {
				$status = "BAD!!!";
				$result = "Fail";
				$report = $report . "\thas a convergence rate of $cvgRates[$ii], theoretically pressure fields should converge at a rate ~ 1\n"
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
