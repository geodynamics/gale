#!/usr/bin/perl
use Term::ANSIColor;
$dir=`pwd`; #print "$dir\n";
$dir =~ /((\/\w+)+\/stgUnderworld\w*\/*).*/;
$UWdir = $1;
#print "UWDIR = $UWdir\n\n";
$properdefargsfile = "$UWdir\/"."properdefargs.txt";
#print "using $properdefargsfile\n\n";
open FILE, "$properdefargsfile";
while (<FILE>){
    /(\w+)\s+(.*)/;
    $getproperdefarg{ $1 } = $2;
}
close FILE;
$defargsfile = "$UWdir\/"."defargs.txt";
open FILE, "$defargsfile";
while (<FILE>){
    /(\w+)\s+(.*)/;
    $getdefarg{ $1 } = $2;
    #print "==============> $1 $2 <=============\n";
}
close FILE;

sub getall{
    my ($str) = @_;
    my $res;
    if($str=~/_(\w+)_New/){
	$str = uc $1;
	$str = "$str"."_DEFARGS";
    }
    #if($str eq "STG_CLASS_DEFARGS" || $str eq "STG_OBJECT_DEFARGS" || $str eq "STG_COMPONENT_DEFARGS"){
    if($str eq "STG_CLASS_DEFARGS"){
	print "$str\n";
	$str="SizeT _sizeOfSelf, Type type, Stg_Class_DeleteFunction* _delete, Stg_Class_PrintFunction* _print, Stg_Class_CopyFunction* _copy";
	return "$str";
    }
    else{
	my $nextpart = $getdefarg{ $str };
	print "$str\n";
	if($nextpart eq ''){
	    return "Failed!\n";
	}
	if($nextpart ne ''){
	    $nextpart =~ s/\n//g;
	    $nextpart =~ s/\s+/ /g;
	    $nextpart =~ s/\\//g;
	    $nextpart =~ /\#define\s*\\?\s*(\w+)\s+(\w+)\s*\,*\s*(.*)/;
	    my $rem = $3;
	    #print "$1 $2 $3 $rem\n";
	    my $tmp=$2;
	    $tmp=~s/\s*(\w+.*\w+)\s*/$1/;
	    my $str = &getall($tmp);
	    if($rem ne ''){
		#$res = "$str".", $rem";
		$res = "$str".";\n $rem";
	    }
	    else{
		$res = "$str";
	    }
	    return $res;
	}
	else{
	    return "str $str Failed";
	}
    }
    
}

$arg="$ARGV[0]";

$test = &getall($arg);
$test =~ s/\,\s+/\,\n/g;
$test =~ s/\;\s+/\;\n/g;
print "\n";
@ar = split(/;/,$test);
for $j (0 .. $#ar){
    if($j % 2 == 1){print color 'green'; print "$ar[$j]\n"; print color 'reset';}
    else { print color 'yellow'; print "$ar[$j]\n"; print color 'reset';}
}
#print color 'green';
#print "$test\n";
#print color 'reset';
