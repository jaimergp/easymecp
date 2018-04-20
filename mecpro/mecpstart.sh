#!/bin/bash

# Copyright (c) 2017 Brigham Young University

# See the file license.txt for copying permission.

# Given a .mecp file for submission mecpstart will handle the creation of the .sh submission script

showhelp ()
{
	echo "USAGE:  $(basename $0) [options] [job]"
	echo
	echo "Available options:"
	echo "	-H <hours>	Sets the walltime, and doesn't ask for it anymore."
	echo "	-c <cores>	Sets how many cores you wish to request."
	echo "	-m <mem>	Sets how much memory (in GB) is needed per processor"
	echo "	-d <dir>	Sets the directory that this will create."
	echo "	-e <email>	Emails the given email when the job completes or is stopped."
	echo "	-l		Sets the job submission directory to the current directory."
	echo "	-n		Do not clean up the submission directory."
	echo "	-i		Interactive mode. Will prompt for nodes, cores and"
	echo "			memory and confirm sending the script. (default behavior)"
	echo "	-I		Non-interactive mode: Use the default options instead of"
	echo "			prompting for them."
	echo "	-s		Generate script only. Do not submit it to SLURM."
	echo "	-t		Designates this as a test job"
	echo "	-h		Show this help screen and exit."
	echo "If you do not provide a job, you will be prompted for one instead."
}

interactive=1
script_only=0
test_job=0
clean=1

while getopts "H:c:m:d:e:hIistln" opt; do
	case $opt in
		i) interactive=1;;
		I) interactive=0;;
		s) script_only=1;;
		t) test_job=1;;
		H) hours=$OPTARG;;
		c) cores=$OPTARG;;
		m) origmem=$OPTARG;;
		d) job_dir=$OPTARG;;
		l) job_dir='.'; clean=0;;
		n) clean=0;;
		e) email=$OPTARG;;
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

shift $(($OPTIND - 1))
job=$1

if test -z $job ; then
	echo "What is the name of your job?"
	read job
fi

# check for correct mecp file format
python ~/scripts/mecpinput.py -f $job.mecp
ret=$?
echo $ret
# TODO: python numpy import error when running python < 2.7
#if [ $ret -ne 0 ]; then
# 	echo 'Correct the MECPFormat Exception before submitting again.'
# 	exit 1
#fi

# check for cores and memory in job.mecp
if test -z $cores ; then
	cores=$(awk -F '=' '/^%nproc/ {gsub(/[ \t\r]+/, "", $2); print $2}' $job.mecp)
fi

if test -z $origmem ; then
	origmem=$(awk -F '=|[^0-9]+$' '/^%mem/ {gsub(/[ \t\r]+/, "", $2); print $2}' $job.mecp)
fi

# if interactive, and cores and memory were not in job.mecp
if test $interactive -eq 1 ; then
	if test -z $hours ; then
		echo "How many hours do you expect this job to take? (Default is 24)"
		read hours
	fi
	
	if test -z $cores ; then
		echo "How many cores do you need? (Default is 12)"
		read cores
	fi

	if test -z $origmem ; then
		echo "How much memory do you need, in GB? (Default is 12)"
		read origmem
	fi
fi

# if not interactive, set cores and memory to DEFAULTS
if test -z $hours ; then
	hours=24
fi

if test -z $cores ; then
	cores=12
fi

if test -z $origmem ; then
	origmem=12
fi

mem=$(( (origmem * 1024 + 1536) / cores ))M  # Give an extra 1.5 GB

walltime=`hours2time $hours`

echo "Using the following settings: "
  echo "Walltime:	$walltime"
  echo "Memory:		$origmem GB"
  echo "Cores:		$cores"

if [ $interactive == 1 ]; then
	echo "Is this correct? yes, no or quit. (y/n/q)"
	read submit
	case $submit in
		y) ;;
		n) 
			echo "How many hours do you expect this job to take? (Default is 24)"
			read hours
			echo "How many cores do you need? (Default is 12)"
			read cores
			echo "How much memory do you need, in GB? (Default is 12)"
			read origmem
			mem=$(( (origmem * 1024 + 1536) / cores ))M  # Give an extra 1.5 GB
			walltime=`hours2time $hours`;;
		*) exit;; 
	esac
fi
  
if test -z $job_dir ; then
	job_dir=$job
fi

if test ! -e $job_dir; then
	echo "Creating folder for your job..."
	mkdir $job_dir
fi

if [ $job_dir != '.' ]; then
	if [ $clean == 1 ]; then
		echo "Cleaning up from previous submission of the same name..."
		rm -rfv $job_dir/*
	fi

	cp "$job.mecp" "$job_dir/"
fi

echo "Generating job script"
sed "s/~WALLTIME~/$walltime/g;s/~JOB~/$job/g;s/~MEM~/$mem/g;s/~CORES~/$cores/g;"\
 ~/.templates/mecp.template > "$job_dir/$job.sh"

if [ $test_job == 1 ]; then
	sed 's/~ISTEST~/#SBATCH --qos=test/g' "$job_dir/$job.sh" -i
else
	sed 's/~ISTEST~//g' "$job_dir/$job.sh" -i
fi

if test -z $email ; then
	sed 's/~EMAIL~//g' "$job_dir/$job.sh" -i
else
	sed "s/~EMAIL~/\
#SBATCH --mail-user=$email   # email address\n\
#SBATCH --mail-type=END\n\
#SBATCH --mail-type=FAIL\
/g" "$job_dir/$job.sh" -i
fi

if [ $interactive == 1 ]; then
	echo "File generated:"
	head -n 20 "$job_dir/$job.sh"
	echo '...'
	echo "Would you like to sbatch the file as well? (y/n)"
	read submit
	case $submit in
		y) script_only=0;;
		*) script_only=1;; 
	esac
fi

if [ $script_only == 0 ]; then
	echo "Submitting job..."
	pushd $job_dir
	sbatch "$job.sh"
	popd
fi
