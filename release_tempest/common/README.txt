                      Path-O-Gen v1.2 2009
                Temporal Signal Investigation Tool
                              by
                       Andrew Rambaut

              Institute of Evolutionary Biology
                    University of Edinburgh
                      a.rambaut@ed.ac.uk

UNIX / Mac OS X / Linux / Windows README 
a.rambaut@ed.ac.uk - 27 November 2009

Contents:
1) INTRODUCTION
2) INSTALLING AND RUNNING PATH-O-GEN
3) ANALYSING TREES
4) VERSION HISTORY
5) SUPPORT & LINKS
6) ACKNOWLEDGMENTS

___________________________________________________________________________
1) INTRODUCTION

Path-O-Gen is a tool for investigating the temporal signal and 'clocklikeness' of molecular phylogenies. It can read and analyse contemporaneous trees (where all sequences have been collected at the same time) and dated-tip trees (where sequences have been collected at different dates). It is designed for analysing trees that have not been inferred under a molecular-clock assumption to see how valid this assumption may be. It can also root the tree at the position that is likely to be the most compatible with the assumption of the molecular clock.

___________________________________________________________________________
2) INSTALLING AND RUNNING PATH-O-GEN

Mac OS X: To install Path-O-Gen, simply drag the program file to where you normally put applications. Then double click to run.

Windows: To install Path-O-Gen, simply drag the program file to where you normally put applications. Then double click to run.

Linux / UNIX: Copy or move the folder to where you normally put applications and then double click the "pathogen.jar" file (in the lib/ directory) to run or type "./pathogen" at the command-line. 

___________________________________________________________________________
3) ANALYSING TREES

Once Path-O-Gen is running it will ask for a tree file to load. This should be in NEXUS format and should have been constructed using a phylogenetic method that does not assume a molecular clock (such as Neighbor-Joining or Maximum Likelihood or Bayesian methods with the molecular clock option off. It is also important that the trees contain branch lengths as genetic distance (substitutions per site). 

When the tree is loaded you will see a table containg all the taxa (sequence labels). If the sequences are contemporaneous (i.e., not sampled through time) then you can leave this as it is. If the sequences have dates associated with them you can enter them into this table. If the taxon labels have the dates encoded in them, you can use the "Guess Dates" button to try and extract them. The final thing you need to set here is whether the dates are "Since some time in the past" - which they will be for calendar dates or days since some time etc. or "Before the present" - most likely the case for carbon dated samples.

You can now select the "Trees" tab at the top of the window and you will see your tree, along with a table of statistics on the left. The nature of the statistics will depend on whether the tree has contemporaneous tips or dated tips. If it is a contemporaneous tree then the statistics will include the mean and variance of the root-to-tip distances for all tips. If it has dated tips then the table will contain various details of a regression of root-to-tip distances against dates of sampling. For example the slope of this regression is an estimate of the rate of evolution, the X-Intercept an estimate of the time of the most recent common ancestor of the sample. The variance for the contemporaneous trees and the correlation coefficient for the dated tips will give you some idea about how 'clocklike' the data are (i.e., how much variation in rate there is).

Selecting the "Best-fitting root" button at the top left of the window will attempt to find the root of the tree that gives the best fit to the hypothesis that the data have a roughly constant rate of evolution. For contemporaneous trees this will find the root which minimizes the variance of root-to-tip variance. For dated tips this will be the root which maximizes the correlation of root-to-tip distance to sampling date. You can also select the "Root-to-tip" tab which will show you a chart of the distribution of root-to-tip distances (a regression against sampling date for dated tips).

Finally, you can export the tree (rooted as displayed) using the "Export Tree..." option in the file menu and export the raw root to tip data using the "Export Data..." option. To obtain a graphic of the displayed tree or chart, you can use the Print option and then "Save as PDF..." or similar option depending on the operating system being used.

___________________________________________________________________________
4) VERSION HISTORY

---Version 1.2 27 November 2009---

* Added the ability to select points in the plots and the equivalent taxa will be highlighted in the tree (and vice-versa).

* Added a residual plot for time-sampled trees. This shows the distribution of residual from the regression line to look for outliers.

---Version 1.1 23 February 2009---

* Added a more flexible tree viewing component (based on FigTree)

* Tips of a dated tip tree are now shown coloured by their residual from the root to tip regression line (blue: above, red: below, black on the regression).

---Version 1.0 12 February 2009---

* First released verson

___________________________________________________________________________
5) SUPPORT & LINKS

Please email me to discuss any problems:

a.rambaut@ed.ac.uk

___________________________________________________________________________
6) ACKNOWLEDGMENTS

Thanks to the following for supplying code or assisting with the creation or testing of BEAST and its associated software:

	Alexander Alekseyenko
	Erik Bloomquist
	Roald Forsberg
	Joseph Heled
	Simon Ho
	Philippe Lemey
	Gerton Lunter
	Sidney Markowitz
	Tulio de Oliveira
	Oliver Pybus
	Beth Shapiro
	Korbinian Strimmer
	Marc Suchard
	+ numerous other users who have kindly helped make BEAST better.
