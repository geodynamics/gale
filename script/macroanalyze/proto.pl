#!/usr/bin/perl
#use Tie::File;
use Term::ANSIColor;
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
our $defargsfile;
our $properdefargsfile;
our $protofile;
our %getparent;
our %getdefarg;
our %getproperdefarg;
our %getproto;
#################################################################################################################################
$dir=`pwd`; #print "$dir\n";
$dir =~ /((\/\w+)+\/stgUnderworld\w+).*/;
$UWdir = $1;

$filename="$ARGV[0]";

open FILE, "$filename"  or die "can't find file: $!\n";
$hugestr='';

while (<FILE>){

    $hugestr .= "$_";

}

$protofile = "$UWdir\/"."proto.txt";
#print "Protofile is $protofile\n";
open FILE, "$protofile";
while (<FILE>){
    /(\w+)\s(.*)/;
    $getproto{ $1 } = $2;
}
close FILE;
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#### find all functions that are of the form ?X_?New ############################################################################
@functions=$hugestr=~/\w+(?:\**\s+|\s+\**|\s+\**\s+)(\w+New\w*\s*\([^{}]+$cparen)/g;
#            is of form / (?: | | | )( ) / 2nd pair of brackets is the match
#################################################################################################################################
# get fname fbody and fargs and fargslist from the function defn
$i=0;
foreach (@functions){

    $fname[$i]=&getfname($_);
    $fbody[$i]=&getfbodyverbatim($_);
    $fargs[$i]=&getfargs($_);
    #print "$fname[$i] ( $fargs[$i] )\n";
    if($fname[$i] =~ /^_\w+/ && !exists $getproto{$fname[$i]} ){
	#print "Yarrrr!\n";
	open FILE, ">>$protofile";
	print FILE "$fname[$i] $fargs[$i]\n";
	close FILE;
    }
    $i++;

}
#################################################################################################################################
#################################################################################################################################
