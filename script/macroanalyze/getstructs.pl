#!/usr/bin/perl

use Term::ANSIColor;
use Data::Dumper;

$dir=`pwd`;
$dir =~ /((\/\w+)+\/stgUnderworld\w*).*/;
$UWdir = $1;
$LIST="$UWdir"."/script/macroanalyze";

unshift(@INC, $LIST);

our $paren;
our $cparen;
our $getdefinesall;
our $reg;
require 'doomstr.pl';
require 'doominc.pl';
our $dir;
our $UWdir;
our $filename;
our $pc;
our $hugestr;
our $structsfile;

#################################################################################################################################
$dir=`pwd`; #print "$dir\n";
$dir =~ /((?:\/\w+)+\/stgUnderworld\w+)(.*)/;
$UWdir = $1;


$filename="$ARGV[0]";
open FILE, "$filename"  or die "can't find file: $!\n";
$hugestr='';

while (<FILE>){

    $hugestr .= "$_";

}

$structsfile = "$UWdir\/"."structs.txt";


# open FILE, "$structsfile";
# while (<FILE>){
#     /(\w+)\s(.*)/;
#     $getdefarg{ $1 } = $2;
# }
# close FILE;

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#### find all functions that are of the form ?X_?New ############################################################################
@functions=$hugestr=~/\w+(?:\**\s+|\s+\**|\s+\**\s+)(_\w+New\w*\s*\([^{}]+$cparen)/g; ## Only looking at _X_?New functions now i.e. starting with _
#########      is of form / (?: | | | )( ) / 2nd pair of brackets is the match  #################################################
#################################################################################################################################
# get fname fbody and fargs and fargslist from the function defn ################################################################
#######print "filename is $filename\n";  ########################################################################################
$i=0; $goforth = 1;# $numparents = 0;

$headerfile=$filename;   $headerfile =~ s/\.c/\.h/; # get the headerfile 
###########################################################################################################
if( -e  $headerfile){
    
    open FILE, "$headerfile";   $headerstr='';
    while (<FILE>){ $headerstr .= $_; }
    close FILE;
    $headerstr =~ s/\r//g;
    #$headerstr = &cleanstr($headerstr);
    $headerstr =~ s{/\* (?: (?!\*/). )* \*/}{}gxs;
    $headerstr =~ s{//.*}{}g; # get rid of any c++ style comments as well
    foreach (@functions){
	$fname[$i]=&getfname($_);
        ## Construct a macro name
	$Xstructname = $fname[$i];
	$Xstructname =~ s/(.*)_New/_$1/;
	
	print "Structname $Xstructname\n";
	$getmacrodefine=&createdefinestruct($Xstructname);
	
	while ($headerstr =~ /$getmacrodefine/g){
	    $hashdefine = $1;
	    $hashdefine = &cleanstr($hashdefine); # clean it up a bit and make it fit on one line
	    $hashdefine =~ s/\\//g;
	    $hashdefine = &cleanstr($hashdefine);
	    $hashdefine =~ s/^\s*//;
	    $hashdefine =~ s/(\W)struct(\W)/$1$2/g;
	    $hashdefine =~ s/\{//g;
	    $hashdefine =~ s/\}//g;
	    $hashdefine =~ s/\s*$//;
	    $hashdefine = &cleanstr($hashdefine);
	    $hashdefine =~ s/\#define//;
	    $hashdefine =~ s/^\s*//;
	    print "Define $hashdefine\n";
	    open FILE, ">>$structsfile";
	    ##print FILE "$Xstructname $hashdefine\n";
	    print FILE "$hashdefine\n";
	    close FILE;
	}
	$i++;
    }# foreach (@functions)

}# if headerfile exists
exit;

#################################################################################################################################
#################################################################################################################################
### Now check for macros
$i=0;
# cool it looks like ALL macros as args end in 'ARGS'
# Lets just worry about the DEFARGS part i.e. the args to X
# because we will construct PASSARGS from the DEFARGS anyway
foreach $args (@Xargs){
    @Xargset = split(/,/,$args);
    #print Dumper(@Xargset);
    $Xhasmacro[$i]=0;
    $Xisplugin[$i]=0;
    if( $#Xargset == 0 ){
	#print "Only has one argument\n";
	if( ($args !~ /\w+\W+\w+/) && ($args !~ /void/) ){
	    #print color 'blue'; print "Is a macro $args in $filename\n"; print color 'reset';
	    if($args !~ /ARGS/){
		open FILE, ">>log.txt";
		print FILE "Looks like macro but doesn't have ARGS in name $Xname[$i] $args in $filename\n";
		close FILE;
	    }
	    $Xhasmacro[$i] =1;
	}
	else{
	    if($args =~ /ARGS/){ 
		open FILE, ">>log.txt"; print FILE "has 2 parts to args but has ARGS in it $Xname[$i] $args $filename\n"; close FILE; }
	    #print "Is a plugin $Xname[$i] $args   $filename\n";
	    $Xisplugin[$i]=1;
	}
    }
    $i++;
}

#Now we know who has macros and who are probably plugins
## Process the functions that have macros now and build the defarg hash table
$i=0;
foreach $Xfuncname (@Xname){
    ###########################################################################################################
    if($Xhasmacro[$i]){
	$Xargs[$i] =~ /(\w+)/;
	$macname = $1;
	$propermacname = "$Xclass[$i]"."_DEFARGS";
	print ">>>>>>>>>$propermacname<<<<<<<<<<<< >>>>>>>>>>>>>>>$macname<<<<<<<<<<<<<<<<\n";
	$headerfile=$filename;   $headerfile =~ s/\.c/\.h/; # get the headerfile 
        ###########################################################################################################
   	if( -e  $headerfile){
	    open FILE, "$headerfile";   $headerstr='';
	    while (<FILE>){ $headerstr .= $_; }
	    close FILE;
	    $getmacrodefine=&createdefinestr($macname);
	    $headerstr =~ /$getmacrodefine/g;
	    $hashdefine = $1;
	    $hashdefine =~ s/\s+/ /g; # clean it up a bit and make it fit on one line
	    #print "$hashdefine\n";
            ###########################################################################################################
	    if($hashdefine eq ''){ 
		print "Could not find definition for $macname. Try again by trying proper macro name\n";
		# Try what the name should be
		$macname=$propermacname;
		$getmacrodefine=&createdefinestr($macname);
		$headerstr =~ /$getmacrodefine/g;
		$hashdefine = $1;
		$hashdefine =~ s/\s+/ /g; # clean it up a bit and make it fit on one line
		#print "$hashdefine\n";
		if($hashdefine eq ''){ 
		    #print "Could not find definition for $macname. Try again by remangling macro name to proper but with ARGS instead of DEFARGS\n";
		    # Try what it might be
		    $macname = "$Xclass[$i]"."_ARGS";
		    $getmacrodefine=&createdefinestr($macname);
		    $headerstr =~ /$getmacrodefine/g;
		    $hashdefine = $1;
		    $hashdefine =~ s/\s+/ /g; # clean it up a bit and make it fit on one line
		    #print "$hashdefine\n";
		    if($hashdefine eq ''){ 
			open FILE, ">>log.txt";	
			print FILE "Could not find definition for $macname used by $Xname[$i] in $filename. We can assume is equivalent to proper macro that doesn't exist yet!\n"; 
			close FILE;
			# OK...if a macroname was used but not defined in a function and was not found
			# then it must be defined elsewhere but can be taken to be the same as the proper macro name because it works
			if( !exists $getdefarg{$propermacname} ) {
			    $Xargs[$i] =~ /(\w+)/;  $macname = $1; # get original macro name from the Xargs again
			    open FILE, ">>$defargsfile"; print FILE "$propermacname $macname \#define $propermacname $macname\n"; close FILE;
			    # this one won't appear in the properdefargs file because we cannot define it in the proper way yet
			}
		    }
		}
	    }# if($hashdefine eq '')
	    ###########################################################################################################
	    if( !exists $getdefarg{$macname} && $hashdefine ne '' ){# then write out the current defn to defargs file
		open FILE, ">>$defargsfile";
		print FILE "$macname $propermacname $hashdefine\n";
		if($macname ne $propermacname){
		    print FILE "$propermacname $macname $hashdefine\n";
		}
		close FILE;
	    }# if macro name doesn't exist already and we found a defn for it
	    if( !exists $getdefarg{$propermacname} && $hashdefine ne '' ){# then write out the current defn to defargs file
		open FILE, ">>$defargsfile";
		print FILE "$propermacname $propermacname $hashdefine\n";
		close FILE;
	    }# if proper macro name doesn't exist already and we found a defn for it
	    ###########################################################################################################
	}# if headerfile exists
        ###########################################################################################################
	else{
	    open FILE, ">>log.txt"; print FILE "Where is the goddam header file for $filename??\n"; close FILE;
	}# if headerfile exists
	###########################################################################################################

    }## if function has a macro for args
    else{## function does not have a macro in its args but we are going to check headerfile anyway to see if a macro is defined there.
	## Some macros use macros that are defined but not used in their own class's args to its function
	if($Xisplugin[$i] == 0){## and not a plugin
	    $propermacname = "$Xclass[$i]"."_DEFARGS";
	    $macname = $propermacname;
	    $headerfile=$filename;   $headerfile =~ s/\.c/\.h/; # get the headerfile 
	    ###########################################################################################################
	    if( -e  $headerfile){
		open FILE, "$headerfile";   $headerstr='';
		while (<FILE>){ $headerstr .= $_; }
		close FILE;
		$getmacrodefine=&createdefinestr($macname);
		$headerstr =~ /$getmacrodefine/g;
		$hashdefine = $1;
		$hashdefine =~ s/\s+/ /g; # clean it up a bit and make it fit on one line
		#print "$hashdefine\n";
		###########################################################################################################
		if($hashdefine eq ''){ 
		    #print "Could not find definition for $macname. Try again by remangling macro name to proper but with ARGS instead of DEFARGS\n";
		    # Try what it might be
		    $macname = "$Xclass[$i]"."_ARGS";
		    $getmacrodefine=&createdefinestr($macname);
		    $headerstr =~ /$getmacrodefine/g;
		    $hashdefine = $1;
		    $hashdefine =~ s/\s+/ /g; # clean it up a bit and make it fit on one line
		    #print "$hashdefine\n";
		}# if($hashdefine eq '')
		## if we still don't find a macro here it doesn't matter..probably doesn't exist
		###########################################################################################################
		if( !exists $getdefarg{$macname} && $hashdefine ne '' ){# then write out the current defn to defargs file
		    open FILE, ">>$defargsfile";
		    print FILE "$macname $propermacname $hashdefine\n";
		    if($macname ne $propermacname){
			print FILE "$propermacname $macname $hashdefine\n";
		    }
		    close FILE;
		}# if macro name doesn't exist already and we found a defn for it
		if( !exists $getdefarg{$propermacname} && $hashdefine ne '' ){# then write out the current defn to defargs file
		    open FILE, ">>$defargsfile";
		    print FILE "$propermacname $propermacname $hashdefine\n";
		    close FILE;
		}# if proper macro name doesn't exist already and we found a defn for it
		###########################################################################################################
	    }# if headerfile exists
	    ###########################################################################################################
	    else{
		open FILE, ">>log.txt"; print FILE "Where is the goddam header file for $filename??\n"; close FILE;
	    }# if headerfile exists
	    ###########################################################################################################
	}
    }## does not have macro in functions args
    $i++;
    ###########################################################################################################
}
#################################################################################################################################
#################################################################################################################################
