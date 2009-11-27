#!/bin/bash
#Run mercurial identify on all repos
#Usage: ./identify-all.sh

#Process current dir and first level of subdirs excluding .hg
wd=`pwd`
for f in `find . -maxdepth 1 -type d \( ! -iname ".hg" \) | sort -f`
do
    #skip if no .hg folder
    ls $f/.hg &> /dev/null
    if [ "${?}" -ne "0" ]; then
        continue 
    fi

    cd $f
    echo "`hg identify` $f"
    cd $wd
done
