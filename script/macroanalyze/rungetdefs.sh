#!/bin/sh

find -name "*.c" -exec ./script/macroanalyze/getdefargs.pl \{\} \;
