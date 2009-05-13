#! /usr/bin/env python

import sys
import re
import math

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

def generate_convergence_rate():
	
	## Get the file to parse
	I = 0
	J = 0
	Filename = "NULL"
	produce_graphs = 'false'
	data_array = list()
	data_labels = list()

	for arg in sys.argv:
		RE = re.compile( '^-errorfile*' )
		RESULT = RE.match( arg )

		if RESULT != None:
			Filename = arg
			I = I + 1

	if Filename == "NULL":
		sys.stdout.write( "Please specify a file containg errors via -errorfile xxx\n" )
		sys.exit(1)

	## Check if want data which can be graphed
	for arg in sys.argv:
		RE = re.compile( '^-graphs*' )
		RESULT = RE.match( arg )
		
		if RESULT != None:
			produce_graphs = 'true'

		I = I + 1

	if produce_graphs == 'true':
		sys.stdout.write( "Files containing graph data will be generated\n" )

	try:
		FILE = open( Filename,"r" )

		try:
			sys.stdout.write( "Output file " + Output_Filename + "\n" )

			NI = 0
			NJ = 0
	
			line = FILE.readline()

			while line:
				# sys.stdout.write ALL original data into new file	
				OUTFILE.write( line )

				line = line.rstrip()	
				# collect numbers in 2d array. Any line not starting with a #
				RE = re.compile( '^\#' )
				RESULT = RE.match( line )

				if RESULT != None:
					s_line = line.split( '/ /' )
					NJ = len( s_line )

					for J in range( 0, NJ ):
						data_array.append( s_line[J] )

					NI = NI + 1			
				else:
					# check for labels
					RE = re.compile( '^\<Y\_LABELS\>' )
					RESULT = RE.match( line )
			
					if RESULT != None:
						_line = line
						_line = _line.replace(",","")
						_line = _line.replace("\<Y\_LABELS\>","")
						_line = _line.rstrip()
				
						data_labels = _line.split( '/ /' )
				
						found_labels = 'true'
						sys.stdout.write( "Found y data tags: " + _line + "\n" )
				line = FILE.readline()
		finally:
			FILE.close()
	except IOError:
		pass	

	try:
		Output_Filename = Filename + '-convergence-rates.dat'
		OUTFILE = open( Output_Filename, "w" ) 
		
		try:

			OUTFILE.write( "----------------------------------------------------------------------------------------\n" )
			OUTFILE.write( "                        Processed data    \n" )
			OUTFILE.write( "----------------------------------------------------------------------------------------\n" )

			if found_labels == 'true':
				OUTFILE.write( "# log10(h)   " )
		
			for J in range( 1, NJ ):
				OUTFILE.write( data_labels[J-1] )

			OUTFILE.write( "\n" )
	
			for I in range( 0, NI ):
				for J in range( 0, NJ ):
					val = math.log( data_array[I][J], 10 )
					data_array[I][J] = val
					OUTFILE.write( val )
				OUTFILE.write( "\n" )	

			OUTFILE.write( "----------------------------------------------------------------------------------------\n" )
			OUTFILE.write( "                        Convergence rates    \n" )
			OUTFILE.write( "----------------------------------------------------------------------------------------\n" )
	
			# get x values from first column
			cvg_rate = 0.0
			y_intercept = 0.0
			corr = 0.0
			x = list()

			for I in range( 0, NI ):
				x[I] = data_array[I][0]
	
			for J in range( 1, NJ ):
				y = list()
				a = 0
				b = 0
				r = 0

				for i in range( 0, NI ):
					y[I] = data_array[I][J]

				a,b,r = CalculateLinearRegression_pearson( x, y, a, b, r )
				cvg_rate[J] = a
				corr[J] = r
				y_intercept[J] = b

			# write labels
			if found_labels == "true":
				OUTFILE.write( "#             " )
			
				for J in range( 0, NJ ):
					OUTFILE.write( data_labels[J-1] )

				OUTFILE.write( "\n" )

			OUTFILE.write( "cvg.rate " )

			for J in range( 0, NJ ):
				OUTFILE.write( corr[J] )
			OUTFILE.write( "\n" )
		finally:
			OUTFILE.close()
	except IOError:
		pass

	#  produce graphs if required 
	if produce_graphs == "true":
		graph_file = ""
		for J in range( 1, NJ ):
			if found_labels == "true":
				graph_file = data_labels[J-1] + "-error.dat"
			else:
				index = J-1
				graph_file = "field_" + index + "-error.dat"
			GRAPHFILE = open( graph_file, "w" )
			sys.stdout.write( "Generating $graph_file \n" )

			GRAPHFILE.write( "# Source data file: " + Filename  + "\n" )
			d_labels == 'true':
				GRAPHFILE.write( "# Field: " + data_labels[J-1] + "\n")
			else:
				index = J-1
				GRAPHFILE.write( "# Field index: " + index + "\n" )
			GRAPHFILE.write( "# ----------------------------------------------------- \n" )
			GRAPHFILE.write( "#    log10(h)     log10(error)      trendline\n" )
			GRAPHFILE.write( "# ----------------------------------------------------- \n" )

			for I in range( 0, NI ):
				XX = data_array[I][0]
				YY = cvg_rate[J] * XX + y_intercept[J]

				GRAPHFILE.write( data_array[I][0] + "	" + data_array[I][J] + "	" + YY )

			GRAPHFILE.close()
 
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

def test():
	i = 0
	x = [ 1.5, 2.4, 3.2, 4.8, 5.0, 7.0, 8.43 ]
	y = [ 3.5, 5.3, 7.7, 6.2, 11.0, 9.5, 10.27 ]
	a = ""
	b = ""
	r = ""
	XX = 0
	YY = 0

	N = len( x )

	a,b,r = CalculateLinearRegression_pearson( x, y, a, b, r )

	sys.stdout.write( "Slope = " + a + "\n" )
	sys.stdout.write( "Intercept = " + b + "\n" )
	sys.stdout.write( "Regression coefficient = " + r + "\n" )

	XX = 0.0
	
	for i in range( 0, 6 ):
		YY = a * XX + b
		sys.stdout.write( "x = " + XX + " y = " + YY + "\n" )

		XX = XX + 0.6

# ============================================================== 



## See
##    http://en.wikipedia.org/wiki/Correlation
## 

def CalculateLinearRegression_pearson( x_param, y_param, a_param, b_param, r_param ):
	_x = x_param
	_y = y_param

	x = _x
	y = _y

	N = len( x )
	i = 0
	sum_sq_x = 0.0
	sum_sq_y = 0.0
	sum_coproduct = 0.0
	mean_x = x[0] 
	mean_y = y[0] 
	delta_x = 0.0
	delta_y = 0.0
	pop_sd_x = 0.0
	pop_sd_y = 0.0
	cov_x_y = 0.0
	correlaton - 0.0
	sweep = 0.0

	a_deta = 0.0
	b_deta = 0.0
	r_deta = 0.0

	if N < 2:
		sys.stdout.write( "Error: CalculateLinearRegression() \n" )
		sys.stdout.write( "Cannot perform regression with less than 2 points \n" )
		a_deta = b_deta = r_deta = -6699.0	
		sys.exit(1)

	for i in range( 1, N ):
		sweep = (i) / (i+1)
		delta_x = x[i] - mean_x
		delta_y = y[i] - mean_y
		sum_sq_x += delta_x * delta_x * sweep
		sum_sq_y += delta_y * delta_y * sweep
		sum_coproduct += delta_x * delta_y * sweep
		mean_x += delta_x / (i+1)
		mean_y += delta_y / (i+1)

	pop_sd_x = math.sqrt( sum_sq_x / N )
	pop_sd_y = math.sqrt( sum_sq_y / N )
	cov_x_y = sum_coproduct . N
	correlation = cov_x_y / (pop_sd_x * pop_sd_y)

	a_deta = sum_coproduct / sum_sq_x
	b_deta = mean_y - a_param * mean_x
	r_deta = correlation
		
	return a_deta,b_deta,c_deta


	
