                    BEAST v1.10.4 2002-2018
        Bayesian Evolutionary Analysis Sampling Trees
                              by
      Alexei J. Drummond, Andrew Rambaut & Marc Suchard

                Department of Computer Science
                     University of Auckland
                    alexei@cs.auckland.ac.nz

              Institute of Evolutionary Biology
                    University of Edinburgh
                      a.rambaut@ed.ac.uk

               David Geffen School of Medicine
            University of California, Los Angeles
                      msuchard@ucla.edu


Last updated: a.rambaut@ed.ac.uk - 11th November 2018

Contents:
1) INTRODUCTION
2) INSTALLING BEAST
3) CONVERTING SEQUENCES
4) RUNNING BEAST
5) ANALYZING RESULTS
6) NATIVE LIBRARIES
7) SUPPORT & LINKS
8) ACKNOWLEDGMENTS

___________________________________________________________________________
1) INTRODUCTION

BEAST (Bayesian evolutionary analysis sampling trees) is package for
evolutionary inference from molecular sequences.

BEAST uses a complex and powerful input format (specified in XML) to
describe the evolutionary model. This has advantages in terms of
flexibility in that the developers of BEAST do not have to try and predict
every analysis that researchers may wish to perform and explicitly provide
an option for doing it. However, this flexibility means it is possible to
construct models that don't perform well under the Markov chain Monte Carlo
(MCMC) inference framework used. We cannot test every possible model that
can be used in BEAST. There are two solutions to this: Firstly, we  supply
a range of recipes for commonly performed analyses that we know should work
in BEAST and provide example input files for these (although, the actual
data can also produce unexpected behavour). Secondly, we provide advice and
tools for the diagnosis of problems and suggestions on how to fix them:

<http://beast.community/>

BEAST is not a black-box into which you can put your data and expect an
easily interpretable answer. It requires careful inspection of the output
to check that it has performed correctly and usually will need tweaking,
adjustment and a number of runs to get a valid answer. 

___________________________________________________________________________
2) INSTALLING BEAST

BEAST requires a Java Virtual Machine to run. Many systems will already
have this installed. It requires at least version 1.6 of Java to run. The
latest versions of Java can be downloaded from:

<http://java.sun.com/>

If in doubt type "java -version" to see what version of java is installed
(or if it is installed at all).

Mac OS X will already have a suitable version of Java installed.

Within the BEAST.v1.10.x package will be the following directories:
Directory       Contents
doc/            documentation of BEAST
examples/       some example NEXUS and XML files
lib/            Java & native libraries used by BEAST
native/         some C code to compile into native libraries
bin/            Scripts for the corresponding OS

___________________________________________________________________________
3) CONVERTING SEQUENCES

A program called "BEAUti" will import data in NEXUS format, allow you to
select various models and options and generate an XML file ready for use in
BEAST.

To run BEAUti simply double-click the "BEAUti v1.10.x" application in the BEAST
folder. If this doesn't work then you may not have Java installed correctly.
Try opening an MS-DOS window and typing:

	java -jar lib/beauti.jar

See also the separate BEAUti README.txt document.

__________________________________________________________________________
4) RUNNING BEAST

To run BEAST simply double-click the "BEAST v1.10.x" application in the BEAST
folder. You will be asked to select a BEAST XML input file.

Alternatively open a Command window and type:

	java -jar lib/beast.jar input.xml

Where "input.xml" is the name of a BEAST XML format file. This file can
either be created from scratch using a text editor or be created by the
BEAUti program from a NEXUS format file.

For documentation on creating and tuning the input files look at the
documentation and tutorials on-line at:

Help -      <http://beast.bio.ed.ac.uk/>
FAQ -       <http://beast.bio.ed.ac.uk/FAQ/>
Tutorials - <http://beast.bio.ed.ac.uk/tutorials/>

The latest manual can be downloaded from here:

<http://code.google.com/p/beast-mcmc/>

BEAST arguments:
     -verbose        "Give verbose XML parsing messages"
     -warnings       "Show warning messages about BEAST XML file"
     -strict         "Fail on non-conforming BEAST XML file"
     -window         "Provide a console window"
     -options        "Display an options dialog"
     -working        "Change working directory to input file's directory"
     -seed           "Specify a random number generator seed"
     -prefix         "PREFIX", "Specify a prefix for all output log filenames"
     -overwrite      "Allow overwriting of log files"
     -errors         "Specify maximum number of numerical errors before stopping"
     -threads        "The number of computational threads to use (default auto)"
     -java           "Use Java only, no native implementations"
     -threshold      "Full evaluation test threshold (default 0.1)"
     -beagle         "Use BEAGLE library if available (default on)"
     -beagle_off     "Don't use BEAGLE library"
     -beagle_info          "BEAGLE: show information on available resources"
     -beagle_auto          "BEAGLE: automatically select fastest resource for analysis"
     -beagle_order         "BEAGLE: set order of resource use"
     -beagle_instances     "BEAGLE: divide site patterns amongst instances"
     -beagle_CPU           "BEAGLE: use CPU instance"
     -beagle_GPU           "BEAGLE: use GPU instance if available"
     -beagle_SSE           "BEAGLE: use SSE extensions if available"
     -beagle_SSE_off       "BEAGLE: turn off use of SSE extensions"
     -beagle_cuda          "BEAGLE: use CUDA parallization if available"
     -beagle_opencl        "BEAGLE: use OpenCL parallization if available"
     -beagle_single        "BEAGLE: use single precision if available"
     -beagle_double        "BEAGLE: use double precision if available"
     -beagle_async         "BEAGLE: use asynchronous kernels if available"
     -beagle_scaling       "BEAGLE: specify scaling scheme to use"
     -beagle_rescale       "BEAGLE: frequency of rescaling (dynamic scaling only)"
     -mc3_chains    "number of chains"
     -mc3_delta     "temperature increment parameter"
     -mc3_temperatures     "a comma-separated list of the hot chain temperatures"
     -mc3_swap      "frequency at which chains temperatures will be swapped"
     -load_dump     "Specify a filename to load a dumped state from"
     -dump_state    "Specify a state at which to write a dump file"
     -dump_every    "Specify a frequency to write a dump file"
     -version       "Print the version and credits and stop"
     -help"         "Print this information and stop"

For example:

     java -jar lib/beast.jar -seed 12345 -overwrite input.xml


Note: MCMC is not a black-box, and poor prior choices and starting values
      can affect inference, publication probability and career outcomes.

___________________________________________________________________________
5) ANALYZING RESULTS

We have produced a powerful graphical program for analysing MCMC log files
(it can also analyse output from MrBayes and other MCMCs). This is called
'Tracer' and is available from the Tracer web site:

<http://tree.bio.ed.ac.uk/software/tracer>

We have now included the "loganalyser" program again in order to analyse
log and tree files without the need for tracer.

Additionally, two new programs are distributed as part of the BEAST
package: LogCombiner & TreeAnnotator. LogCombiner can combine log or tree
files from multiple runs of BEAST into a single combined results file
(after removing appropriate burn-ins). TreeAnnotator can summarize a sample
of trees from BEAST using a single target tree, annotating it with
posterior probabilities, HPD node heights and rates. This tree can then be
viewed in a new program called 'FigTree' which is available from:

<http://tree.bio.ed.ac.uk/software/figtree>

___________________________________________________________________________
6) NATIVE LIBRARIES

BEAST uses the BEAGLE library to provide fast, native, likelihood 
calculations. This library is distributed independently and should
be installed prior to using BEAST.

BEAGLE can be downloaded from here:
<https://github.com/beagle-dev/beagle-lib>

Further information about using BEAGLE with BEAST can be found here:
<http://beast.community/beagle>

___________________________________________________________________________
7) SUPPORT & LINKS

BEAST is an extremely complex program and as such will inevitably have
bugs. Please email us to discuss any problems:

<alexei@cs.auckland.ac.nz>
<a.rambaut@ed.ac.uk>
<msuchard@ucla.edu>

We would encourage you to join the BEAST users' mailing-list to get
notifications of updates and bugs. At a later date this may be expanded to
be a discussion-list so that users can exchange ideas and help. You can
join the mailing list here:

<http://groups.google.com/group/beast-users>

The website for beast is here:

<http://beast.community/>

Source code distributed under the GNU Lesser General Public License:

<https://github.com/beast-dev/beast-mcmc/>

___________________________________________________________________________
8) ACKNOWLEDGMENTS

Thanks to the following for supplying code or assisting with the creation
or testing of BEAST (this is a woefully incomplete list - apologies to the
inadvertently absent names):

	Alex Alekseyenko, Daniel Ayres, Guy Baele, Trevor Bedford, Peter Beerli,
	Filip Bielejec, Erik Bloomquist, Aaron Darling, Gytis Dudas, Roald Forsberg,
	Matthew Hall, Jess Hedge, Joseph Heled, Simon Ho, Sebastian Hoehna, Mark Holder,
	Denise Kuehnert, Philippe Lemey, Paul Lewis, Wai Lok Sibon Li, Gerton Lunter,
	Sam Lycett, Sidney Markowitz, Vladimir Minin, Michael Defoin Platel,
	Oliver Pybus, Tulio de Oliveira, Jayna Ragwani, Beth Shapiro,
	Korbinian Strimmer, David Swofford, Chieh-Hsi Wu, Walter Xie, Ziheng Yang
