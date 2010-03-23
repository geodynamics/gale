#!/bin/sh

find -name "*.c" -exec ./script/macroanalyze/proto.pl \{\} \;
