#!/bin/bash

#Copyright (c) 2015 Brigham Young University

#See the file license.txt for copying permission.

# Recursively checks for convergence
# Case convergence found, creates a new .com file with the converged information
# Case convergence not found, the user can choose to run the job again


showhelp ()
{
   echo "USAGE:  $(basename $0) [options]"
   echo
   echo "Available options:"
   echo "   -d <directory>  set the directory you want to check for convergence"
   echo "                   still works if given current directory"

}

while getopts "hd:" opt; do
   case $opt in
      d) job_dir=$OPTARG;;
      h)
         showhelp
         exit 0
         ;;
      \?)
         echo "Invalid option: -$OPTARG"
         showhelp
         exit 1
         ;;
      :)
         echo "No argument given for -$OPTARG"
         showhelp
         exit 1
         ;;
   esac
done

_next="_next"
PROGNAME=$0

if test $job_dir ;
   then
   #echo "directory given"
   #dir=${PWD##*/}
if [ $job_dir != '.' ]
      then
      if test -d "$job_dir"
         then
            cd "$job_dir"
            $PROGNAME
            cd ".."
            exit 0
         else
            echo "$job_dir not found"
            exit 1
      fi
   else
      $PROGNAME
   fi

else
   jobname=${PWD##*/}          # working directory is the jobname
   if test -d "$jobname$_next"
      then
         cd "$jobname$_next"
         $PROGNAME
      else
         mecpstatus_helper
   fi


fi
