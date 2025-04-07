# flicker-neuralcodes

This repository contains custom code and for the analyses performed in "40 Hz sensory stimulation enhances CA3-CA1 coordination and prospective coding during navigation in a mouse model of Alzheimer’s disease", Abigail L. Paulson, Lu Zhang, Ashley M. Prichard, Annabelle C. Singer, 2025. Data collected from this experiment can be found here: https://doi.org/10.6084/m9.figshare.28736681

### hardware and software

This code was developed and tested on a Windows 10 machine using Matlab 2019b and RStudio 2022.12.0+353 with R version 4.2.2. 

### code

The main scripts used for analysis and figure generation are `cf_RunAnalyses.m ` and `cf_AllFigures.m`. <br>
Additional analysis scripts can be found in `/scripts/`. <br>
Analysis for PPC, WPLI, and PSD can be found in `/scripts/PPC_WPLI_PSD/`. <br>
R code for statistical analysis can be found in `/scripts/R_code/` <br>

Most functions are located in `/functions/`, but some may be found in the following repositories:<br>
  https://github.com/gmbcrazy/GenMatCode <br>
  https://github.com/singerlabgt/singer-lab-common-functions

