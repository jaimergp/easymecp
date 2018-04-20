#!/bin/bash

#Copyright (c) 2017 Brigham Young University

#See the file license.txt for copying permission.

#echo "$1"
result="$1"
dotnext=".next"
mecp=".mecp"
old="OldChk"
A_chk="_A.chk"
B_chk="_B.chk"
origin=".origin"

echo $result
echo "To receive an email when the job finishes, enter your email. Otherwise, press enter."
	read email

echo "We will continue searching for crossing points"

# copy original .mecp into an origin file if it does not exist yet
if ! test -f $result$mecp$origin
	then
		touch $result$mecp$origin
		cp -v $result$mecp $result$mecp$origin
fi
# copy the .next.mecp into the .mecp file
cp $result$dotnext$mecp $result$mecp
# remove the .next .mecp
rm $result$dotnext$mecp

# case preopt has already been run, remove preopt
if grep -q 'pre_opt' "./$result$mecp";
	then
		echo "Optimization previously completed. Removing preopt keyword."
		sed -i '/pre_opt*/c\' ./$result$mecp
fi

# TODO: should we use previous walltime/memory/nproc?

# handling user email
if [ -z "$email" ]
then
    mecpstart -d $result -l -n $result
else
    mecpstart -e $email -l -n $result 
fi


