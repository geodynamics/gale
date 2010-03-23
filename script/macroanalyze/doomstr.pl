#our ($paren);
$paren = qr/  ##THIS one seems to work OK.
      \(
        (?:
           [^()]+  # Not parens
         | 
           (??{ $paren })  # Another balanced group (not interpolated yet)
        )*
      \)
    /x;
$cparen = qr/  ##THIS one seems to work OK.
      \{
        (?:
           [^{}]+  # Not parens
         | 
           (??{ $cparen })  # Another balanced group (not interpolated yet)
        )*
      \}
    /x;
$getdefinesall = qr/ ([ \t]*\#define[ \t]+\w+ARGS.*\\[ \t]*\n
                   (?:.*\\[ \t]*\n)* # anything that has a backslash as last character i.e. a continuation slash
                    .*\n)  # then get very next line only
                   /x;
# $getdefinesrem = qr/ [ \t]*\#define[ \t]+\w+ARGS[ \t]*\\[ \t]*\n
#                    ((?:.*\\[ \t]*\n)* # anything that has a backslash as last character i.e. a continuation slash
#                     .*\n)  # then get very next line only
#                    /x;
$reg = qr/^\s*(\w+(?:[ \t]+\w+)*)
          ([ \t\*]+)
          ((?:\w|[\[\]]*)+)\s*$
          # (\w+)
         /x;
