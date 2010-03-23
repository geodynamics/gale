#!/bin/sh

find -name "*.c" -exec ./script/macroanalyze/getstructs.pl \{\} \;
