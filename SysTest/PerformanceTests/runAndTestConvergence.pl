#!/usr/bin/perl -w
#
use strict;

sub getInput;
sub runTests;
sub testConvergence;

######  Main program    ######

# resolutions to run convergence tests at
my @resolutions = ( 16, 32, 64, 128, 256 );

# 1) Run the same xml at the resolutions specified in '@resolutions'
#    these will generate error data
runTests();

# 2) Use a linear regression algorithm on the data generated
#    and check the results for good behavior and acceptable convergence
testConvergence();


######  End main program ######

sub runTests {
	my $arg;
	my $res;
	my $command;
	my $xmlFile = "NULL";

# check if xml exists
	foreach $arg (@ARGV) {
		if( $arg =~ m/.*\.xml$/ ) { $xmlFile = $arg; }
	}

	if( $xmlFile eq "NULL" ) { die "No xml file specified\n" ; }

######
# TODO: Need to check for an executable
#       or I can make a link to the executable
######

	print "Testing the convergence of $xmlFile with resolutions:";
	foreach $res (@resolutions) {
		print " $res";
	}
	print "\n\n";

	foreach $res (@resolutions) {
		$command = "./Underworld $xmlFile --elementResI=$res --elementResJ=$res --pluginData.appendToAnalysisFile=True";
		print "Runing $command" ;
		`$command >>runTests.log 2>>runTests.error`;
		print "\n\t\t\t\t....finished\n\n";
	}
}

sub testConvergence {
	# check if the cvg can be found
	my $command = "ls *\.cvg 2>/dev/null";
	my $cvgFile = `$command`;
	chomp $cvgFile;
	if( $cvgFile eq "" ) { die "There is no cvg file found analyse\n" ; }

	$command = "./generate-convergence-rates.pl -errorfile $cvgFile";
	`$command`;

	my $datFile = "$cvgFile-convergence-rates.dat"; 
	if( (-e $datFile) eq 0 ) { die "Can't find the expected file $cvgFile\n"; }
	
	# remove the old file
	$command = "rm $cvgFile";
  `$command`;

  #test convergence numbers
	open(INPUT, "<$datFile") || die "Can't open the expected file $datFile\n" ; 
	my @keys;
	my @cvgRates;
	my @correlations;
	my $needLabels = 1;
	my $line;
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

	my $ii;
	my $jj;
	my $length;
	my $status;
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
			$report = $report . "\thas a correlation of $correlations[$ii], below the accepted threshold of 0.99\n";
		}
		if( $keys[$ii] =~ m/^Velocity/ ) {
			if( $cvgRates[$ii] < 1.6 ) {
				$status = "BAD!!!";
				$report = $report . "\thas a convergence rate of $cvgRates[$ii], theoretically velocity fields should converge at a rate ~ 2\n"
			}
		} elsif ($keys[$ii] =~ m/Pressure/ ) {
			if( $cvgRates[$ii] < 0.9 ) {
				$status = "BAD!!!";
				$report = $report . "\thas a convergence rate of $cvgRates[$ii], theoretically pressure fields should converge at a rate ~ 1\n"
			}
		}
		print "$status\n$report";
	}
}
