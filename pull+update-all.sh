#!/bin/bash
#Run mercurial pull and update in all directories
#Usage: ./pull+update-all.sh [username]
#Prompts for password

check-errs()
{
  # Function. Parameter 1 is the return code
  # Para. 2 is text to display on failure.
  if [ "${1}" -ne "0" ]; then
      echo "ERROR # ${1} : ${2}"
      exit ${1}
  fi
}

#check for username provided in hg paths url (look for @ symbol)
#if present don't rewrite repository urls with username/pass 
expr match "`hg paths`" '.*\(@\).*' &> /dev/null
if [ "${?}" -ne "0" ]; then
    #get username from command line argument, defaults to current user
    if [ $# -ne 1 ]
    then
        echo "Using currently logged in user name: `whoami`"
        user=`whoami`
    else
        user=$1
    fi

    #get password
    echo -n "Enter password for $user: "
    stty -echo
    read password
    stty echo
    echo ""

    #setup username and password
    login="https://$user:$password@"
else
    login="https://"
fi

#Process current dir and first level of subdirs excluding .hg
wd=`pwd`
for f in `find . -maxdepth 1 -type d \( ! -iname ".hg" \)`
do
    #skip if no .hg folder
    ls $f/.hg &> /dev/null
    if [ "${?}" -ne "0" ]; then
        continue 
    fi

    cd $f
    echo "-------- Processing [ $f ] ---------------------------------------"

    #strip start of hg paths output to get repository name
    paths=`hg paths | sed 's/default = https:\/\///'`

    hg pull $login$paths
    check-errs $? "hg pull failed"
    hg update 
    check-errs $? "hg update failed"
    cd $wd
done
