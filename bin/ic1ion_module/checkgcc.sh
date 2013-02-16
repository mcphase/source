#!/bin/sh

basedir=`dirname $(echo $0)`

if test ! -x $basedir/checkgcc.exe
then
   $2 -static-libgcc $basedir/checkgcc.cpp -o $basedir/checkgcc.exe 1>/dev/null 2>/dev/null
fi

$basedir/checkgcc.exe $1
if [ "$?" -ne "0" ]
then
  if [ $# -gt 2 ]
  then
    echo Sorry, you version of GCC is older than $1.
    echo ic1ion produces inaccurate numerical results with these compilers.
    echo So you must upgrade your compiler to a newer version.
    exit 1
  fi
fi
