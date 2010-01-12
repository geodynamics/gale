#!/usr/bin/perl

use POSIX qw(log10);
use strict;

sub CalculateLinearRegression_pearson;
sub test;
sub generate_convergence_rate;
sub trim($);

generate_convergence_rate();
#test();


## --------------------------------------------------------------------------------------
##      Options:      
##          -errorfile XXX, The file containing the errors you wish the determine convergence rates for.
##          -graphs, Flag to indicate you want to generate a seperate graph for each field. 
##                   Graph file format;        x_values[]    y_values[]     trendline_values[]
##                   If no labels were specified (see below) the file names will be 
##                     field_X-error.dat, where X is the column index for that field. 
##
##          A processed file will be generated with the filename XXX-convergence-rates.dat,
##          where XXX is the file specified with -errorfile XXX
##
##      Notes:
##          - File format. 
##              - Lines starting with a # will be ignored
##              - First column must contain x values.
##              - Subseqent columns contain y values.
##              - Data in each y-column can be given label using the following
##                # <Y_LABELS> vx vy vz p        
##                In this case there should be 4 columns, name vx,vy,vz and p respectively.
##                DO NOT LABEL THE X COLUMN
## --------------------------------------------------------------------------------------




sub generate_convergence_rate 
{
    my ( $arg, $produce_graphs );
    my $filename;
    my ( $I,$J, $NI,$NJ, @data_array );
    my $found_labels = 'false';
    my @data_labels;
    
    # Get file to parse 
    $I = 0;
    $filename = 'NULL';
    foreach $arg (@ARGV) {
        if( $arg =~ m/^\-errorfile/ ) {
            $filename = $ARGV[$I+1];
        }
        $I++;
    }
    if( $filename eq 'NULL' ) {
            print "Please specify a file containg errors via -errorfile xxx\n";
            exit;
    }

    print "Processing file $filename ... \n";



    # Check if want data which can be graphed 
    $I = 0;
    $produce_graphs = 'false';
    foreach $arg (@ARGV) {
        if( $arg =~ m/^\-graphs/ ) {
            $produce_graphs = 'true';
        }
        $I++;
    }
    if( $produce_graphs eq 'true' ) {
            print "Files containing graph data will be generated\n";
    }


    # ------------------------------------------------------------------------------ #

    open( FILE, "$filename" ) || die "Can't open '$filename': $!\n";
    my @data = <FILE>;

    my $out_filename = $filename . '-convergence-rates.dat';
    open( OUTFILE, ">$out_filename" ) || die "Can't open '$out_filename': $!\n";
    print "Output file $out_filename \n";

    
    $NI = 0;
    $NJ = 0;
    my $line;
    foreach $line (@data) {
        # print ALL original data into new file
        print OUTFILE "$line";
        
        chomp( $line );
        # collect numbers in 2d array. Any line not starting with a # 
        if( !($line =~ m/^\#/) ) {
            
            my @s_line = split( / /, $line) ; 
            
            $NJ = @s_line;
            for( $J=0; $J<$NJ; $J++ ) {
                $data_array[$NI][$J] = $s_line[ $J ];
            }
    
            $NI++;
        }
        else {
            # check for labels
            if( $line =~ m/\<Y\_LABELS\>/) {
                my $_line = $line;
                $_line =~ s/\,//g; # strip all commas
                $_line =~ s/^\#//; # strip starting hash
                $_line =~ s/\<Y\_LABELS\>//; # strip <Y_LABLELS>
                chomp( $_line ); # remove linebreak
                my $tline = trim( $_line ); # remove white space
                
                @data_labels = split( / /, $tline );
                
                $found_labels = 'true';
                print "Found y data tags: $tline \n";
            }
        }
    }
    close( FILE );

    print OUTFILE "----------------------------------------------------------------------------------------\n";
    print OUTFILE "                        Processed data    \n";
    print OUTFILE "----------------------------------------------------------------------------------------\n";

    if( $found_labels eq 'true' ) {
        print OUTFILE "# log10(h)   ";
        for( $J=1; $J<$NJ; $J++ ) {
           print OUTFILE "$data_labels[$J-1]     ";
        }
        print OUTFILE "\n";
    }

    for( $I=0; $I<$NI; $I++ ) {
        for( $J=0; $J<$NJ; $J++ ) {
            my $val = log10( $data_array[ $I ][ $J ] );
            $data_array[ $I ][ $J ] = $val;
            printf( OUTFILE "%1.5f ", $val);       
        }
        print OUTFILE "\n";
    }




    print OUTFILE "----------------------------------------------------------------------------------------\n";
    print OUTFILE "                        Convergence rates    \n";
    print OUTFILE "----------------------------------------------------------------------------------------\n";
    
    # get x values from first column 
    my ( @cvg_rate, @y_intercept, @corr );
    my @x;
    for( $I=0; $I<$NI; $I++ ) {
        $x[ $I ]  = $data_array[ $I ][ 0 ];
    }

    for( $J=1; $J<$NJ; $J++ ) {
        my (@y,$a,$b,$r);
        for( $I=0; $I<$NI; $I++ ) {
            $y[ $I ]  = $data_array[ $I ][ $J ];
        }
        
        CalculateLinearRegression_pearson( \@x,\@y, $a,$b,$r );
        $cvg_rate[ $J ] = $a;
        $corr[ $J ] = $r;
        $y_intercept[ $J ] = $b;
    }

    # write labels
    if( $found_labels eq 'true' ) {
        print OUTFILE "#             ";
        for( $J=1; $J<$NJ; $J++ ) {
           print OUTFILE "$data_labels[$J-1]     ";
        }
        print OUTFILE "\n";
    }

    print OUTFILE "cvg. rate   ";
    for( $J=1; $J<$NJ; $J++ ) {
        printf( OUTFILE "%1.2f     ", $cvg_rate[$J] );
    }
    print OUTFILE "\n";
    
    print OUTFILE "corr.       ";
    for( $J=1; $J<$NJ; $J++ ) {
#        printf( OUTFILE "%1.2e ", $corr[$J] );
        printf( OUTFILE "%f ", $corr[$J] );
    }
    print OUTFILE "\n";

    close( OUTFILE );



    # produce graphs if required
 if( $produce_graphs eq 'true' ) {
    my $graph_file;

    # PLOT RAW DATA 
    for( $J=1; $J<$NJ; $J++ ) {
        if( $found_labels eq 'true' ) {
            $graph_file = $data_labels[$J-1] . '-error.dat';
        }
        else {
            my $index = $J-1;
            $graph_file = 'field_' . $index . '-error.dat';
        }
        open( GRAPHFILE, ">$graph_file" ) || die "Can't open '$graph_file': $!\n";
        print "Genrating $graph_file \n";
    
        print GRAPHFILE "# Source data file: $filename \n";    
    
        if( $found_labels eq 'true' ) {
            print GRAPHFILE "# Field: $data_labels[ $J-1 ] \n";    
        }
        else {
            my $index = $J-1;
            print GRAPHFILE "# Field index: $index \n";    
        }
        print GRAPHFILE "# ----------------------------------------------------- \n";    
        print GRAPHFILE "#    log10(h)     log10(error)      trendline       \n";    
        print GRAPHFILE "# ----------------------------------------------------- \n";    
    
        for( $I=0; $I<$NI; $I++ ) {
            my ($XX, $YY);
            $XX = $data_array[ $I ][0];
            $YY = $cvg_rate[ $J ] * $XX + $y_intercept[ $J ];
            
            printf( GRAPHFILE "%1.5f          %1.5f          %1.5f\n", $data_array[ $I ][0], $data_array[ $I ][ $J ], $YY );
        }
    
        close( GRAPHFILE );
    }

 }


}




## Expected output from test
## Slope = 0.904273 
## Intercept = 3.46212 
## Regression coefficient = 0.808257 
## x = 0.00  y = 3.46212 
## x = 0.60  y = 4.00469 
## x = 1.20  y = 4.54725 
## x = 1.80  y = 5.08981 
## x = 2.40  y = 5.63238 
## x = 3.00  y = 6.17494 

sub test 
{
    my $i;
    my @x = ( 1.5, 2.4, 3.2, 4.8,  5.0, 7.0,  8.43 );
    my @y = ( 3.5, 5.3, 7.7, 6.2, 11.0, 9.5, 10.27 );
    my ($a,$b,$r);
    my ($XX, $YY); 
    
    my $N = @x;
    
    #print "@x \n";
    #print "@y \n";
    #print "$N \n";
    
    CalculateLinearRegression_pearson( \@x,\@y, $a,$b,$r );
    
    print "Slope = $a \n";
    print "Intercept = $b \n";
    print "Regression coefficient = $r \n";
        
    $XX = 0.0;
    for( $i=0; $i<6; $i++ ) {
        $YY = $a * $XX + $b;
        print "x = $XX  y = $YY \n";
        
        $XX = $XX + 0.6;
    }
}




# ==============================================================


## See 
##    http://en.wikipedia.org/wiki/Correlation
##
sub CalculateLinearRegression_pearson
{
    my $_x = @_[0];
    my $_y = @_[1];    
    
    my @x = @{$_x};
    my @y = @{$_y};
    
    my $N = @x; 

#    print "@x \n";
#    print "@y \n";
#    print "$N \n";



	my $i;
	my ($sum_sq_x, $sum_sq_y);
	my $sum_coproduct;
	my ($mean_x, $mean_y);
	my ($delta_x, $delta_y);
	my ($pop_sd_x, $pop_sd_y);
	my ($cov_x_y, $correlation);
	my  $sweep;
	
	if ($N<2) {
		print "Error: CalculateLinearRegression() \n";
		print "Cannot perform regression with less than 2 points \n";
		$_[2] = $_[3] = $_[4] = -6699.0;

		die;
	}
	
	
	$sum_sq_x = 0.0;
	$sum_sq_y = 0.0;
	$sum_coproduct = 0.0;
	$mean_x = $x[0];
	$mean_y = $y[0];
	for ($i=1; $i<$N; $i++) {
		$sweep = ( $i ) / ($i+1);
		$delta_x = $x[$i] - $mean_x;
		$delta_y = $y[$i] - $mean_y;
		$sum_sq_x += $delta_x * $delta_x * $sweep;
		$sum_sq_y += $delta_y * $delta_y * $sweep;
		$sum_coproduct += $delta_x * $delta_y * $sweep;
		$mean_x += $delta_x / ($i+1);
		$mean_y += $delta_y / ($i+1) ;
	}
	$pop_sd_x = sqrt( $sum_sq_x / $N );
	$pop_sd_y = sqrt( $sum_sq_y / $N );
	$cov_x_y = $sum_coproduct / $N;
	$correlation = $cov_x_y / ($pop_sd_x * $pop_sd_y);
	
	
	
	$_[2] = $sum_coproduct / $sum_sq_x;
	$_[3] = $mean_y - $_[2] * $mean_x;
	$_[4] = $correlation;
	
	return;
}


# Remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
