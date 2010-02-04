#################################################################################################################################
#################################################################################################################################
sub getargsfromdefarg{
    my ($str) = @_;
    my $res;
    if($str=~/_(\w+)_New/){
	$str = uc $1;
	$str = "$str"."_DEFARGS";
    }
    #if($str eq "STG_CLASS_DEFARGS" || $str eq "STG_OBJECT_DEFARGS" || $str eq "STG_COMPONENT_DEFARGS"){
    if($str eq "STG_CLASS_DEFARGS"){
	#print "str is $str\n";
	$str="SizeT _sizeOfSelf, Type type, Stg_Class_DeleteFunction* _delete, Stg_Class_PrintFunction* _print, Stg_Class_CopyFunction* _copy";
	return "$str";
    }
    else{
	my $nextpart = $getdefarg{ $str };
	#print "str is $str\n";
	if($nextpart ne ''){
	    #print "From defargs.txt \n $nextpart\n";
	    $nextpart =~ s/\n//g;
	    $nextpart =~ s/\\//g;
	    $nextpart =~ s/\s+/ /g;
	    $nextpart =~ /\#define\s+(\w+)\s+(\w+)\s*\,*\s*(.*)/;
	    if($nextpart eq ''){
		return "Failed!\n";
	    }
	    my $rem = $3;
	    #print "$1 $2 $3 \n";
	    my $tmp=$2;
	    $tmp=~s/\s*(\w+)\s*/$1/;
	    $tmp=~s/(\w+)_\w*ARGS/$1_DEFARGS/; ## some of the defs in this file have *_ARGS instead of *_DEFARGS
	    #print "XX$tmp -- $nextpart"."XX\n";
	    if($tmp =~ /SizeT/){
		return "$tmp"." $rem";
	    }
	    my $str = &getargsfromdefarg($tmp);
	    if($rem ne ''){
		$res = "$str".", $rem";
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
sub getall{
    my ($str) = @_;
    my $res;
    if($str=~/_(\w+)_New/){
	$str = uc $1;
	$str = "$str"."_DEFARGS";
    }
    #if($str eq "STG_CLASS_DEFARGS" || $str eq "STG_OBJECT_DEFARGS" || $str eq "STG_COMPONENT_DEFARGS"){
    if($str eq "STG_CLASS_DEFARGS"){
	#print "str is $str\n";
	$str="SizeT _sizeOfSelf, Type type, Stg_Class_DeleteFunction* _delete, Stg_Class_PrintFunction* _print, Stg_Class_CopyFunction* _copy";
	return "$str";
    }
    else{
	my $nextpart = $getproperdefarg{ $str };
	print "str is $str\n";
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
		$res = "$str".", $rem";
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
sub getlocationofArg{# returns first location of arg in list
    my($arg, @list) = @_;
    my $i;
    my $loc=-1;
    for $i (0 ..$#list){
	if($arg eq $list[$i]){
	    $loc = $i;
	    return $loc;
	}
    }
    return $loc;
}
sub getlocationArrayofArgs{
    my ($a, $l) = @_;
    my @args;
    my @list;
    my $i;
    my @loc;
    @args = @$a;
    @list = @$l;
    for $i (0 .. $#args){
	$loc[$i]=&getlocationofArg($args[$i],@list);
	if($loc[$i] != -1){
	    $list[$loc[$i]] = "$list[$loc[$i]]"."removed"; # remove arg from list so we don't find it twice if we have repeated args
	}
    }
    return @loc;
}
sub getlocationArrayofArgsNEWOLD{
    my ($ta, $a, $tb, $b) = @_;
    my @args;
    my @list;
    my $i;
    my @loc;
    my $tmp;
    @targsA=@$ta;
    @targsB=@$tb;
    @argsA = @$a;
    @argsB = @$b;
    for $i (0 .. $#targsA){
	$tmp = "$targsA[$i]"."$argsA[$i]";
	$tmp =~ s/_//g;
	$tmp =~ s/\*//g;
	$tmp = uc $tmp;
	$args[$i]=$tmp;
    }
    for $i (0 .. $#targsB){
	$tmp = "$targsB[$i]"."$argsB[$i]";
	$tmp =~ s/_//g;
	$tmp =~ s/\*//g;
	$tmp = uc $tmp;
	$list[$i]=$tmp;
    }  
    for $i (0 .. $#args){
	$loc[$i]=&getlocationofArg($args[$i],@list);
	if($loc[$i] != -1){
	    $list[$loc[$i]] = "$list[$loc[$i]]"."removed"; # remove arg from list so we don't find it twice if we have repeated args
	}
    }
    return @loc;
}
sub cleanargs{
    my (@args) = @_;
    my $i=0;
    foreach $term (@args){
	$term =~ s/^\s*(.+)/$1/;
	if($term =~ /^(.+)\s+$/){ $term = $1;}
	$args[$i] = $term;
	$i++;
    }
    return @args;
}
sub capandreduceargs{
    my (@args) = @_;
    my $i=0;
    foreach $term (@args){
	$term =~ s/^\s*(.+)/$1/;
	if($term =~ /^(.+)\s+$/){ $term = $1;}
	$term =~ s/_//g;
	$term = uc $term;
	$args[$i] = $term;
	$i++;
    }
    return @args;
}
sub gettypes{
    my ($key) = @_;
    my $args;
    my $argsarray;
    my $term;
    my $stars;
    my $word;
    my @type=();
    my $i;

    #print "$key in gettypes\n";
    #print Dumper(@type);
    #print " ..in gettypes\n";
    $args = $getproto{$key};
    $args =~ s/\r//g;
    #print "In getypes $key gives\n$args\n";
    @argsarray = split(/,/,$args);
    #print Dumper(@argsarray);
    $i=0;
    if($#argsarray > 0){

    foreach $term (@argsarray){
	#print "$term\n";
	$term =~ s/^\s*(.+)/$1/;
	if($term =~ /^(.+)\s+$/){ $term = $1;}
	if($term =~ /$reg/){
	    #print "$1 $2 $3\n";
	    $stars = $2;
	    $word = $1;

	    #$var[$i] = "$3";

            #print "$word >>$3<<\n";

	    $stars =~ s/\s+//g;
	    $word =~ s/\s+/ /g;
	    $type[$i] = "$word$stars";
	    #print "$type[$i] >$var[$i]<\n";
	}
	else{
	    print "Empty arg? >>$args<<\n";
	    
	    exit;
	}
	#print "$term\n";

	
	$i++;
    }#foreach
    }
    return @type;
}
sub getvars{
    my ($key) = @_;
    my $args;
    my $argsarray;
    my $term;
    my @var=();
    my $i;

    #print "$key in getvars\n";
    #print Dumper(@var);
    #print ".. in getvars\n";
    $args = $getproto{$key};
    $args =~ s/\r//g;
    #print "$args\n";
    @argsarray = split(/,/,$args);
    #print Dumper(@argsarray);
    $i=0;
    if($#argsarray > 0){

    foreach $term (@argsarray){
	#print "$term\n";
	$term =~ s/^\s*(.+)/$1/;
	if($term =~ /^(.+)\s+$/){ $term = $1;}
	if($term =~ /$reg/){
	    #print "$1 $2 $3\n";
	    #$stars = $2;
	    #$word = $1;

	    $var[$i] = "$3";

            #print "$word >>$3<<\n";

	    #$stars =~ s/\s+//g;
	    #$word =~ s/\s+/ /g;
	    #$type[$i] = "$word$stars";
	    #print "$type[$i] >$var[$i]<\n";
	}
	else{
	    print "Empty arg in getvars? >>$args<<\n";
	    print "filename $filename\n";
	    exit;
	}
	#print "$term\n";

	
	$i++;
    }#foreach
    }
    return @var;
}
sub gettypesfromstr{
    my ($args) = @_;
    my $argsarray;
    my $term;
    my $stars;
    my $word;
    my @type=();
    my $i;

    $args =~ s/\r//g;
    #print "$args\n";
    @argsarray = split(/,/,$args);
    #print Dumper(@argsarray);
    $i=0;
    if($#argsarray > 0){

    foreach $term (@argsarray){
	#print "$term\n";
	$term =~ s/^\s*(.+)/$1/;
	if($term =~ /^(.+)\s+$/){ $term = $1;}
	if($term =~ /$reg/){
	    #print "$1 $2 $3\n";
	    $stars = $2;
	    $word = $1;

	    #$var[$i] = "$3";

            #print "$word >>$3<<\n";

	    $stars =~ s/\s+//g;
	    $word =~ s/\s+/ /g;
	    $type[$i] = "$word$stars";
	    #print "$type[$i] >$var[$i]<\n";
	}
	else{
	    print "Empty arg? gettypesfromstr >>$args<<\n";
	    print "filename $filename\n";
	    exit;
	}
	#print "$term\n";

	
	$i++;
    }#foreach
    }
    return @type;
}
sub getvarsfromstr{
    my ($args) = @_;
    my $argsarray=();
    my $term;
    my @var=();
    my $i;

    $args =~ s/\r//g;
    if($args =~ /\,/){## need at least one comma for a list
	@argsarray = split(/,/,$args);
    }
    else{## handle case where we have a list of one element
	$argsarray[0] = $args;
	$term=$args;
	$term =~ s/^\s*(.+)/$1/;
	if($term =~ /^(.+)\s+$/){ $term = $1;}
	if($term =~ /$reg/){

	    $var[0] = "$3";

	}
	#print "HERE >>$args<< >>$var[0]<<\n\n";
	return @var;
    }
    #print Dumper(@argsarray);
    $i=0;
    if($#argsarray > 0){

    foreach $term (@argsarray){
	#print "$term\n";
	$term =~ s/^\s*(.+)/$1/;
	if($term =~ /^(.+)\s+$/){ $term = $1;}
	if($term =~ /$reg/){

	    $var[$i] = "$3";

	}
	else{
	    print "Empty arg in getvarsfromstr? >>$args<<\n";
	    print "filename $filename\n";
	    exit;
	}
	#print "$term\n";

	
	$i++;
    }#foreach
    }
    return @var;
}
sub getclass{
    my ($fn) = @_;
    my $class = $fn;

    if($class =~ /_*(\w+)_[A-Za-z]*New/){ 
	$class = $1;
	$class = uc $class;
    } else { $class = 'x'; }
    
    return $class;
}
sub getfargsverbatim{
    my ($fn) = @_;
    my $args = $fn;
    $args =~ s/\w+New[^()]*($paren)[^{}]*.*/$1/s;
    $args =~ s/\((.*)\)/$1/s;
    return $args;
}
sub cleanstr{
    my ($str) = @_;
    $str =~ s/^\((.*)\)$/$1/s;
    $str =~ s{/\* (?: (?!\*/). )* \*/}{}gxs;  #match /*any*/ as long as 'any' does not contain '*/'
    $str =~ s/assert\s*$paren\s*\;//; # get rid of pesky assert functions
    $str =~ s{//.*}{}g; # get rid of any c++ style comments as well
    $str =~ s{\n([ \t]*\n)*}{\n}g;
    $str =~ s{^([ \t]*\n)+}{}g;
    $str =~ s/\n//g;
    $str =~ s/[ \t]+/ /g;
    return $str;
}
sub getfargs{
    my ($fn) = @_;
    my $args = $fn;
    $args =~ s/\w+New[^()]*($paren)[^{}]*.*/$1/s;
    $args =~ s/^\((.*)\)$/$1/s;
    $args =~ s{/\* (?: (?!\*/). )* \*/}{}gxs;  #match /*any*/ as long as 'any' does not contain '*/'
    $args =~ s{//.*}{}g; # get rid of any c++ style comments as well
    $args =~ s{\n([ \t]*\n)*}{\n}g;
    $args =~ s{^([ \t]*\n)+}{}g;
    $args =~ s/\n//g;
    $args =~ s/[ \t]+/ /g;
    return $args;
}
sub getfargsany{
    my ($fn) = @_;
    my $args = $fn;
    $args =~ s/\w+[^()]*($paren)[^{}]*.*/$1/s;
    $args =~ s/^\((.*)\)$/$1/s;
    $args =~ s{/\* (?: (?!\*/). )* \*/}{}gxs;  #match /*any*/ as long as 'any' does not contain '*/'
    $args =~ s{//.*}{}g; # get rid of any c++ style comments as well
    $args =~ s{\n([ \t]*\n)*}{\n}g;
    $args =~ s{^([ \t]*\n)+}{}g;
    $args =~ s/\n//g;
    $args =~ s/[ \t]+/ /g;
    return $args;
}
sub getfname{
    my ($fn) = @_;
    my $name = $fn;
    $name =~ s/(^\w+New\w*)\s*\(.*$/$1/s;
    return $name;
}
sub getfnameany{
    my ($fn) = @_;
    my $name = $fn;
    $name =~ s/(^\w+)\s*\(.*$/$1/s;
    return $name;
}
sub getfbodyverbatim{
    my ($fn) = @_;
    my $body = $fn;
    $body =~ s/\w+New[^{}]*($cparen)[^{}]*/$1/s;
    return $body;
}
sub getfbodyverbatimany{
    my ($fn) = @_;
    my $body = $fn;
    $body =~ s/\w+[^{}]*($cparen)[^{}]*/$1/s;
    return $body;
}
sub getfbody{
    my ($fn) = @_;
    my $body = $fn;
    $body =~ s/\w+New[^{}]*($cparen)[^{}]*/$1/s; #match /*any*/ as long as 'any' does not contain '*/'
    $body =~ s/assert\s*$paren\s*\;//; # get rid of pesky assert functions
    $body =~ s{/\* (?: (?!\*/). )* \*/}{}gxs;
    $body =~ s{\#if 0(?:(?!\#endif).)*\#endif}{}gs;
    $body =~ s{\/\/.*}{}g;
    $body =~ s{^[ \t]*\n([ \t]*\n)*}{\n}g;
    #$body =~ s{^([ \t]*\n)+}{}g; #removes all blank lines
    return $body;
}
sub getfbodyany{
    my ($fn) = @_;
    my $body = $fn;
    $body =~ s/\w+[^{}]*($cparen)[^{}]*/$1/s; #match /*any*/ as long as 'any' does not contain '*/'
    $body =~ s/assert\s*$paren\s*\;//; # get rid of pesky assert functions
    $body =~ s{/\* (?: (?!\*/). )* \*/}{}gxs;
    $body =~ s{\#if 0(?:(?!\#endif).)*\#endif}{}gs;
    $body =~ s{\/\/.*}{}g;
    $body =~ s{^[ \t]*\n([ \t]*\n)*}{\n}g;
    #$body =~ s{^([ \t]*\n)+}{}g; #removes all blank lines
    return $body;
}
sub getYnotinX{#returns Y-(Y && X)
    my ( $Y, $X) = @_;
    my %original = ();
    my @subY  = ();
    map { $original{$_} = 'x' } @$X;
    @subY  = grep { !defined  $original{$_} } @$Y; ## i.e. $_ takes values from the @Y array and sees if they are in the original hash table
    return @subY ; # returns a subset of @Y
}
sub getYinX{#returns (Y && X)
    my ( $Y, $X) = @_;
    my %original = ();
    my @subY  = ();
    map { $original{$_} = 'x' } @$X;
    @subY  = grep { defined  $original{$_} } @$Y; ## i.e. $_ takes values from the @Y array and sees if they are in the original hash table
    return @subY ; # returns a subset of @Y
}
sub createdefinestr{
    my ($macro) = @_;
    my $def;
    $def = qr/ ([ \t]*\#define[ \t]+$macro.*\n ## the space between the [] is important; it is a literal space here
                                       (?:.*\\[ \t]*\n)* # anything that has a backslash as last character i.e. a continuation slash
                                             #^ needed a . here instead of \s because \s seems to match \n as well while . doesn't unless we use s on the
                                             # end of the regexp i.e. "xs" instead of "x"
                                       .*\n)  # then get very next line only
                                /x;
    return $def;
}
sub createdefinestruct{
    my ($macro) = @_;
    my $def;
    $def = qr/ ([ \t]*\#define[ \t]+$macro\s+[\s\\]* ## the space between the [] is important; it is a literal space here
                                       (?:.*\\[ \t]*\n)* # anything that has a backslash as last character i.e. a continuation slash
                                             #^ needed a . here instead of \s because \s seems to match \n as well while . doesn't unless we use s on the
                                             # end of the regexp i.e. "xs" instead of "x"
                                       .*\n)  # then get very next line only
                                /x;
    return $def;
}
sub repAwithBinCretC{#replace A with B in C and return new C
    my ($a,$b,$c)=@_;
    my $ameta=quotemeta($a);
    $c =~ s/$ameta/$b/;
    return $c;
}
1;
