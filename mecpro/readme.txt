==---
Foreword
========---

If you have any questions, comments, bug reports, or feature requests, you can
email me at lily.carlson@byu.edu
    - Lily Carlson

======---
Installation
============---

Extract the files from the archive (you have probably already done that, since
you are reading this file) and put them on the supercomputer somewhere. 

If you have a previous version of mecpro, please put the new version in a new
directory. The installation script will put files in all the right
places for you. Just note that the source code is contained in the new
directory. You can remove the old source code by deleting the old source 
directory.

Once installed, the MECPro commands will work from any directory. We 
recommend keeping your source directory separate from your job directory.

To install, navigate to the directory containing the extracted source code 
and run the following command:
 > source install.sh

You can also run 'bash install.sh', but in doing so, you will have to run:
 > source ~/.bash_profile
Afterward, or log out and back in again.

If you ever change the mecpro source code, scripts or templates you can 
simply navigate the the source code directory and run the following command 
to update all the scripts and templates:
 > source install.sh
Running this won't "Double-install" or anything silly like that. There are
guards which check for that, so you don't need to delete folders or anything.


===============---
Tools in this package
=====================---

The following commands will be installed on your account. You can just type them
all by themselves. There is no need to precede them with ~/scripts/.

mecpstart : Creates submission scripts and instantly submits jobs to the
supercomputer. Eliminates the need to manually create job scripts, thus saving
you a minute or two per job. Very flexible options and easy to interface with
other scripts. To run mecpstart, navigate the to the directory that contains
the .mecp file you want to run, and type the command where the jobname does
not include the .mecp suffix:
  > mecpstart <jobname>

More information on mecpstart usage to follow.

mecpstatus : Returns the status of a job. You can choose to check the status
of a single directory using the -d flag. Without the flag, the script will 
check the status of all the subdirectories in the current directory.
To run mecpstatus navigate to the directory containing the job you ran:
  > mecpstatus
  > mecpstatus -d [path/to/job/directory]
If you have a previous version of MECPro, we mecpstatus will still work.

mecpnext : If a convergence was not found, or the job did not finish, mecpnext
will use mecpstart to continue the job from the last step in the same directory, 
and concatenate to the same <jobname>.log file. You must run mecpnext within the 
directory of the job you want to continue.

mecpdata : This script will generate a <jobname>.csv file from the <jobname>.log 
file. It will extract Step, Energy of First Spin State, Energy of Second Spin State, 
Difference in Energy, Max Gradient, RMS Gradient, Max Displacement, and RMS 
Displacement.
  > mecpdata -f [jobname.log]

mecpopt is also installed to a directory called bin in your home directory, but
chances are you won't ever be running it directly.

All other files in this package are libraries/modules and other dependencies.

===========---
Input file format
=================---

Input files for this program are somewhat free-format, but can easily resemble a
gaussian input file. At the same time, they were designed with the idea of being
easy to generate- in fact, unfinished jobs will create an input file you can run
to resume the calculation. They will look ugly and will be in a strange order.
That is normal behavior that happens because of the way inputs were implemented.

An input file is divided into "sections," which are most easily designated with
square brackets. (e.g. [geometry]) Section headers are not case sensitive (e.g.
GEOMETRY is the same as geometry and GeoMeTrY). A section ends where another
begins. Some sections can be "inlined" into other sections or can be started
with alternate notations so that inputs can closely resemble typical inputs for
Gaussian.

Some sections have no special format and are just "parameter" sections. They
simply contain a list of parameters and their values. Many of these have
defaults. Parameters can be set with :, =, or whitespace. For example, all of
these formats are valid ways to set a parameter:
basis = um06/gen
maxsteps: 20
maxstepsize 0.2

We highly recommend using the examples as a reference when creating your job.

SECTIONS:

[General] (Parameter Section)
This allows you to set general parameters for the run.
Available parameters:
basis[A|B]:     The basis set used in calculations, one for each of the spin
                states. If you don't give a or b, it will use the same basis set
                for both spin states. The program expects a "u" DFT method to be
                used. (defaults to um06/gen)
charge:         The charge of the molecule (defaults to 0)
spinstates:     The spin states to use in calculation. This value should be two
                numbers delimited by a forward slash (/). (defaults to 1/3)
read_later:     When set to true/on, causes all steps after the first to have
                'guess=read' set for the route. This will NOT override your initial
                guess setting and will reflect itself in the resumable inputs
                (the '.next.mecp' file). A singlet for example will have 
                guess=(mix, read) in the route line. (defaults to false/off)
max_steps:      The maximum number of steps to find a convergence (default: 40)
max_stepsize:   The maximum geometry displacement per step (defaults to 0.1)
pre_opt:        Do a standard geometry optimization before doing the rest of the
                calculation. Give it a spin state and a method/basis set to work
                with, space separated. You can provide 'a' or 'b' for the spin
                state, indicating to use whichever one is used by the first or
                second states, respectively. If you have included the guess 
                keyword in your route, it will delete that keyword to perform
                the optimization.(defaults to "none m06/gen")
show_hessian:   Shows information about the hessian matrix at each step. Valid
                settings are 'none', 'full', and 'sign'. The 'sign' setting will
                show only the sign of the hessian, rather than the values.
                (defaults to none; this is primarily a debug feature.)

[Link] (Parameter Section)
This allows you to specify the LINK0 (% lines) section in the gaussian inputs.
You can also "inline" lines from this section by starting a line with a %.

You do not need to specify the chk. It will be automatically named, and in fact,
will be overridden even if you try to set it.

[Route]
You can now include your basis here, so inputs are no longer much different at
all from Gaussian input files. The program will check for and eliminate
redundancy, however, this operates under the assumption that anything with a /
in it is the method/basis, so this may screw up some existing inputs.

This allows you to customize the route card (# line) in the generated gaussian
input files. It defaults to include "force" and "integral=ultrafinegrid". In 
addition, the default will add "guess=mix" to any singlet. If you specify a 
different guess option it will be combined with mix for the singlet. For example,
if you use guess=read for both spin states, the singlet will force 
guess=(mix, read) to be used on the route line. 

All other options should be specified manually.

Each line includes a single option for the route line. If an option should
only be used for one of the spin states, you can prefix a line with A: or 
B: to signify the first and the second spin states, respectively.

You can also "inline" this section into anywhere in the file by using a line
that starts with # for both states a# for the first state only, and b# for the
second state only. Such lines may include multiple options.

[Cutoffs] (Parameter Section)
This section allows you to customize your convergence criteria.
The available parameters are:
max_grad: the largest component of the gradient (defaults to 0.0007)
rms_grad: the root mean square of the gradient  (defaults to 0.0005)
max_chg:  the largest component of displacement (defaults to 0.004)
rms_chg:  the root mean square of displacement  (defaults to 0.0025)
energy_diff: the difference in energy of the two states (defaults to 0.00005)

[Geometry]
The most important part and the only essential part of the file. Include your
atoms and their coordinates in their usual way. It will work with both atomic
numbers and atomic symbols. (You should even be able to mix and match if you're
crazy like that)

An alternate way to begin this section is with a line matching the pattern
# #/#
Where the #s are numbers. The first one sets the charge, the second two set the
spin states.

Unlike other sections, this one may be ended with a blank line, at which point,
it starts an [extra] section.

[Extra]
This section allows you to put extra lines after your geometry.

=========---
mecpstart Usage
===============---

To submit a Gaussian job (.com) to the supercomputer, run the following command:
 > mecpstart <jobname>
(Where jobname is the name of your .com file, without the .com)

As soon as you run it, it will prompt you for how many hours you need to run
your job, how many cores you need, and how much memory you want.
When this is run, it will make a new folder for your job, copy your input file
into it, prepare a .sh file, and then submit the .sh to the supercomputer.

Here are some "Advanced" command-line options for those of you who are
comfortable with running complex commands...

Here are the most notable options:
 > mecpstart -e <email> <jobname>  # Emails you when the job is done/terminated
 > mecpstart -l <jobname>  # If you don't want it to make a folder

Other options are mostly ideal for interfacing with scripts, as many of them
allow you to bypass the prompting of parameters.

This is the "detailed" usage statement given when -h is given as an option:
USAGE:  mecpstart [options] [job]

Available options:
        -H <hours>      Sets the walltime, and doesn't ask for it anymore.
        -c <cores>      Sets how many cores you wish to request.
        -m <mem>        Sets how much memory (in GB) is needed per processor
        -d <dir>        Sets the directory that this will create.
        -e <email>      Emails the given email when the job completes or is stopped.
        -l              Sets the job submission directory to the current directory.
        -n              Do not clean up the submission directory.
        -i              Interactive mode. Will prompt for nodes, cores and
                        memory and confirm sending the script. (default behavior)
        -I              Non-interactive mode: Use the default options instead of
                        prompting for them.
        -s              Generate script only. Do not submit it to SLURM.
        -t              Designates this as a test job
        -h              Show this help screen and exit.
If you do not provide a job, you will be prompted for one instead.

=========---
Examples
===============---

Two example input files are included with this package. Both examples take less
than an hour to run. (You can really expect it to take around 10 minutes)

example1.mecp
This lists all the available parameters with their titles. Run it with the command:
 > mecpstart example1
Follow the prompts you are given. The job will finish in about 10 minutes and will
converge. Check for convergence with the command
 > mecpstatus -d example1

example2.mecp
This is the same geometry as example1, we've just removed some of the parameters. The
program will use the same defaults as in example1. By changing the max_steps to 7, you
can see how to restart a job. To start example2 run the command:
 > mecpstart -e <your email> example2 
Follow the prompts. When you get the email that your job finishes, run the command:
 > mecpstatus -d example2
You will see that example2 has not finished.
Navigate to ./example2 and run the command:
 > mecpnext
Follow the prompts and the job will continue and converge.

example1-error.mecp and example2-error.mecp are also examples of files that will
fail. You may test them if you want to see what an error looks like.

=========---
Output Files
===============---

Suppose you have a run named foo.mecp. 
You may see any of the following output files:

  foo.mecp          Submission file
  foo.mecp.origin   If the job does not converge and you choose to continue, this 
                    file contains your original submission, while foo.mecp will 
                    contain the point the job most recently continued from.
  foo.log           Information about the calculations at each step and 
                    whether convergence criteria has been met is listed here, and
                    the script mecpdata will pull its data from here.
  foo_A0.com        There will be a .com and .log file for each molecule at each step 
                    (foo_A0.log, foo_B0.com, foo_B0.log, foo_A1.com, ...) If the job
                    converges, you can expect to look at the com file of the last step
                    of either spinstate to find the final geometry.
  foo-opt.com       If pre_opt is chosen, there will be a -opt.com and -opt.log file
                    from the resulting optimization calculation.
  foo_A.chk         For both molecule A and B a checkpoint file is used for 
                    calculations.
  foo.csv           Generated by mecpdata
    
