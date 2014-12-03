#!/usr/bin/perl

# This script is to clean diff files from expected files in glucifer/DrawingObjects/tests that use the DummyOpenGL library
# This cleans out lines with the difference is only a small numerical value

sub min {
    if ($_[0]>$_[1]) {return $_[1]} else {return $_[0]};
}

$diffFile = $ARGV[0];

open( FILE1, "<$diffFile" ) || die "Couldn't open $diffFile";

$tolerance = 0.0001;
$hasFailed = 0;

while( <FILE1> ) {
	$line = $_;
	if ( $line eq substr( $line, 0, 1) ) {
		print $line;
		next;
	}
	# If this line in the diff file 
	if ( $line =~ /[0-9]+,[0-9]+c[0-9]+,[0-9]+/ || $line =~ /[0-9]+c[0-9]+/ ) {

		# Get number of lines to compare
		if ( $line =~ /([0-9]+),([0-9]+)c([0-9]+),([0-9]+)/ ) {
			$linesFromExpectedFileCount = $2 - $1 + 1;
			$linesFromOutputFileCount = $4 - $3 + 1;
			# Make sure that the comparison gave us the same number of lines changed
			if ( $linesFromExpectedFileCount != $linesFromOutputFileCount ) {
				print $line;
				next;
			}
		}
		else {
			$linesFromExpectedFileCount = 1;
			$linesFromOutputFileCount = 1;			
		}			
		
		# Get lines from expected file
		for ( $i = 0 ; $i < $linesFromExpectedFileCount ; $i++ ) {
			$linesFromExpectedFile[$i] = <FILE1>;
		}
		
		#Make sure next line is a '---'
		$dashedLine = <FILE1>;
		if ( $dashedLine ne "---\n" ) {
			die "Comparison is dodgy - got '$dashedLine' when expected '---'.";
		}

		# Get lines from output file
		for ( $i = 0 ; $i < $linesFromOutputFileCount ; $i++ ) {
			$linesFromOutputFile[$i] = <FILE1>;
		}		
		
		# Compare lines
		$hasFailed = 0;
		for ( $i = 0 ; $i < $linesFromOutputFileCount ; $i++ ) {		
			$linesFromExpectedFile[$i] =~ m/[><]\s*(.*)\((.*)\)/;				
			$expectedFileFunction 	   = $1;
			@expectedFileArgs          = split(/,/, $2);
			
			$linesFromOutputFile[$i] =~ m/[><]\s*(.*)\((.*)\)/;
			$outputFileFunction      = $1;
			@outputFileArgs          = split(/,/, $2);			
			
			# Check to make sure that the functions are the same
			if ( $expectedFileFunction ne $outputFileFunction ) {
				$hasFailed = 1;
			}
			
			# Check each argument and make sure that it is within a certain tolerance
			$numberOfArgs = @expectedFileArgs;
			for ( $j = 0 ; $j < $numberOfArgs ; $j++ ) {
				$error = abs($expectedFileArgs[ $j ] - $outputFileArgs[ $j ]);

				#scale by number
				$scaleFactor = min($expectedFileArgs[ $j ], $outputFileArgs[ $j ]);
				if ( $scaleFactor != 0 ) {
					$error /= $scaleFactor;
				}

				if ( $error > $tolerance ) {
					$hasFailed = 1;
				}
			}
			
			# if one of the tests failed then just print out all the text as is
			if ( $hasFailed == 1 ) {
				print $line;
			
				for ( $i = 0 ; $i < $linesFromOutputFileCount ; $i++ ) {
					print $linesFromExpectedFile[$i];
				}			
				print "---\n";
				for ( $i = 0 ; $i < $linesFromOutputFileCount ; $i++ ) {
					print $linesFromOutputFile[$i];
				}
			}
		}
	}
	else {
		print $line;
	}
}

close( FILE1 );
