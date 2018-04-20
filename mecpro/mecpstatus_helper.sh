#!/bin/bash

name=${PWD##*/} 
if test -f "$name.mecp"; then
  if ls *.out > /dev/null 2>&1; then
    if test -f "$name.log"; then
      if grep -q 'CONVERGENCE FOUND!' "$name.log"; then
        echo "Converged"|column -t
      else
        awk -v N=1 -v pattern="Run terminated" '{i=(1+(i%N));if (buffer[i]&& $0 ~ pattern) print buffer[i]; buffer[i]=$0;}' "$name.log"
        #awk '/STEP/{a=$0}END{print a}' OCH.log
      fi
    else
      echo "No logs. Either running or had a format error."
    fi
  else
    echo 'Never submitted or pending on queue'
  fi
else 
  echo "No matching .mecp file"
fi
