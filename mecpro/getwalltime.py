#!/usr/bin/python

#Copyright (c) 2015 Brigham Young University

#See the file license.txt for copying permission.

# This exists because doing this in python is 10x easier than doing it in bash.
# You can use it from a bash script and assign the result to variables using `backquotes`

import sys

def hours2time(time):
    time = float(time)
    hours = int(time)

    time = (time - hours) * 60
    minutes = int(time)

    seconds = int ((time - minutes) * 60)
    return hours, minutes, seconds

if __name__ == '__main__':
    print "%02d:%02d:%02d" % hours2time(sys.argv[1])
