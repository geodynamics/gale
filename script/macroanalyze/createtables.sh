#!/bin/sh
rm -f proto.txt defargs.txt parentchildhashtable.txt structs.txt
./script/macroanalyze/runproto.sh
./script/macroanalyze/rungetdefs.sh
./script/macroanalyze/rungetstructs.sh
