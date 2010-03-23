#!/usr/bin/perl

$pc = "parentchildhashtable.txt";

open FILE, "$pc";
while (<FILE>){
    /(\w+)\s(\w+)/;
    $getparent{ $1 } = $2;
}

$parent="$ARGV[0]";

print "  $parent";

while ($parent ne ''){
    if($getparent{ $parent } ne ''){
	print " => $getparent{ $parent }";
    }
    $parent = $getparent{ $parent };
}
print "\n";
