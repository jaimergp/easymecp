#!/bin/bash
# Copyright (c) 2015 Brigham Young University
# See the file license.txt for copying permission.

# Created by Lily H. Carlson
# Parses through the .log file of a complete mecp job to extract data for a spreadsheet

while getopts "f:" opt; do
   case $opt in
      f) file=$OPTARG
           if test -f "$file"
           then
           echo "$file found"
          else
           echo "$file does not exist"
           exit 1
          fi
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


RESULTS="`basename $file .log`"
echo 'Step, Energy of First Spin State, Energy of Second Spin State, Difference in Energy, Max Gradient, RMS Gradient, Max Displacement, RMS Displacement' > $RESULTS.csv  


awk '{if ($1 == "Energy" && $3 =="first") print $7}' $file > energy1.data
awk '{if ($1 == "Energy"&& $3=="second") print $7}' $file > energy2.data 
awk '{if ($1 == "Max" && $2=="Gradient") print $4 }' $file > maxGrad.data 
awk '{if ($1 == "RMS" && $2=="Gradient") print $4 }' $file > rmsGrad.data 
awk '{if ($1 == "Max" && $2=="Displacement:") print $3 }' $file > maxDisp.data 
awk '{if ($1 == "RMS" && $2=="Displacement:") print $3 }' $file > rmsDisp.data 


paste energy1.data energy2.data maxGrad.data rmsGrad.data maxDisp.data rmsDisp.data|awk -v step=0 '{print step++,",",$1,",",$2,",",$1-$2,",",$3,",",$4,",",$5,",",$6}' >> $RESULTS.csv
rm *.data
