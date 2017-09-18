# Version History and notes

Scripts written by: Tim Thurman. timothy.j.thurman@gmail.com
Created Feb. 2016

##Scripts put onto GitHub Sept. 18, 2017

I haven't worked on this program in a bit and it isn't under active development, but wanted to put it in another place where people can get to it.

##V1.1 uploaded to the Barrett google drive on June 3, 2016-

Fixes the problem of overestimating the # of bp sequenced per individual and % of genome coverage when doing long-read, paired-end sequencing on short fragments.

##v1.0 uploaded to the Barrett google drive on Feb. 5, 2016-

The two R scripts in the folder allow you to simulate the # of ddRAD fragments generated with various enzyme pairs and size selection ranges, and to understand how that relates to total amount of sequencing and % of genome sequenced. These functions require the SimRAD and dplyr packages. NB- As the script is now, when estimating genome coverage for paired-end sequencing with small fragments, the calculations do not account for possible overlap of sequencing reads in the middle of the fragment. Thus, the % genome coverage will be overestimated.

Full notes and instructions are in the ddRAD_enzyme_optimization R script. Please email me to report any bugs/issues.

Planned updates:
function to allow for calculating # of sequencing lanes needed and cost of sequencing based on the number of individuals in a library and the desired level of coverage.
Fixing overestimation of genome coverage for paired-end sequencing of small fragments.
