#!/bin/bash

showhelp ()
{
   echo "USAGE:  $(basename $0) [options]"
   echo
   echo "Available options:"
   echo "   -d <directory>  set the directory of mecp run outputs"

}

while getopts "hd:" opt; do
   case $opt in
      d) 
         job_dir=$OPTARG
         ;;
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

_next='_next'
#if -d flag left empty use current directory
if test -z $job_dir ; then
  job_dir=$(pwd)
fi

# if the job_dir variable given
if test -e $job_dir ; then
   #check if directory exists
    if test -d "$job_dir" ; then
        cd "$job_dir"
        subdircount=`find . -maxdepth 1 -type d | wc -l`
        if [ $subdircount -eq 1 ]; then
            mecpstatus_helper
        else
          for dir in `ls -d */` ; do 
          echo $dir
          cd "$dir"
          jobname=${PWD##*/}          # working directory is the jobname
          if test -d "$jobname$_next"; then
            mecpcheck -d "$jobname$_next"
          else
            mecpstatus_helper
          fi
          cd ..
          echo ""
          done
        fi
    else
      echo "Directory does not exist"
    fi
else
    echo 'No directory given'
    exit 1
fi
