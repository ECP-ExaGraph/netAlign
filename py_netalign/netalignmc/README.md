---
title: "Multicore codes for network alignment"
layout: project
---

Multicore codes for network alignment
==================================

### [Arif Khan](http://www.cs.purdue.edu/homes/khan58/)
### [David F. Gleich](http://www.cs.purdue.edu/homes/dgleich)
### [Mahantesh Halappanavar](http://www.pnl.gov/science/staff/staff_info.asp?staff_num=7452)
### [Alex Pothen](http://www.cs.purdue.edu/homes/pothen)

_These codes are research prototypes and may not work for you. No promises. But do email if you run into problems._

Download
--------

* [netalignmc.tar.gz](netalignmc.tar.gz) 


Prereqs
-------

* (Matlab) working mex compiler
* (Multi-core codes) Intel C++ compiler


Overview
--------

The package is organized by directory

`data`  
: datasets for the experiments

`data-mtx`  
: datasets in matrix market format for netalign c++ codes

`experiments`
: script files that use and evaluate the algorithms

`matlab`
: code for all the main algorithms

`netalign`
: the multicore C++ network alignment codes 

`private_data`
: a historically named directory for our data files



Usage
-----

With a recent version of Matlab (2011a used for development), 
the best way to use the package is to open Matlab, navigate 
to the `matlab` directory, and then execute the 
following commands to solve the network alignment problem for Figure 1.

    % load A, B, and L for figure 1
    >> load('../data/example_overlap.mat');

    % create S, w, and a variant of L from graph A, B, and L
    >> [S,w,li,lj] = netalign_setup(A,B,L);

    % use S, w, L, and alpha=0, beta=1 and the bp algorithm
    >> x = netalignbp(S,w,0,1,li,lj);

    % x is just a heuristic, so we need to round it
    >> [ma mb mi overlap weight] = mwmround(x,S,w,li,lj);
    % [ma mb] give pairs of matched vertices 
    % mi is the binary matching indicator for L (as li, lj)
    % overlap and weight are the overlap and weight objectives
    
If you want to get more information about how to use the tools,
see the `matlab/demo.m` script.
    

Algorithms
----------

All of the algorithms and codes are contained in the `matlab` directory.

|Algorithm            |Code                  |Source           |
|:--------------------|:---------------------|:----------------|
|IsoRank              |`full_isorank.m`      |Singh et al. 2007|
|SpaIsoRank           |`isorank.m`           | |
|NetAlignBP           |`netalignbp.m`        | |
|Exact enumeration    |`netalign_exact.m`    | |
|NetAlignMR           |`netalignmr.m`        |Klau 2009|
|LP                   |`netalign_lp_prob.m`  |Klau 2009|
|Lagrangean LP        |`netalign_llp.m`      |Klau 2009|

Data
----

|Dataset            |File                                 |Source             | 
|:------------------|:------------------------------------|:------------------|
|`lcsh-rameau`  | `data/lcsh2rameau`   | |
|`lcsh2wiki-small`  |`data-mtx/lcsh-small`   | |
|`lcsh2wiki-full`   |`private_data/lcsh2wiki_full.mat`    | |
|`dmela-scere`      | `data/dmele-scere.mat` | Singh et al. 2007 |
|`musm-homo`        | `data/natalie_graphs.mat`        | Klau 2009 |


#### `lcsh2wiki`

These data sets only come with S, w, li, and lj.  The original datasets 
are considerably larger (hundreds of megabytes).  Please contact 
us if you want them.  None of the experiments require them.

#### `dmela-scere`

These datasets are distributed with IsoRank in the 
file `http://groups.csail.mit.edu/cb/mna/packages/multiway_kpartite.tgz` .

Use the python program `experiments/bioinfo/convert_isorank_data.py` 
to generate dmela-scere.smat dmela.smat and scere.smat which are
L, A, and B, respectively.  The program `readSMAT.m` will load these 
files into Matlab.

See http://groups.csail.mit.edu/cb/mna/ http://groups.csail.mit.edu/cb/mna/ for more about IsoRank.

After converting the data to smat with the above script, execute

    A = readSMAT('dmela.smat');
    B = readSMAT('scere.smat');
    L = readSMAT('dmela-scere.smat');
    [S,w,li,lj] = netalign_setup(A,B,L);
    save '../../data/dmela-scere.mat' A B L S w li lj
    
to save the data for use with the experiments.

#### `musm-homo`

These datasets are distributed with 
[Natalie](https://www.mi.fu-berlin.de/w/LiSA/Natalie) in the file
https://www.mi.fu-berlin.de/wiki/pub/LiSA/Natalie/natalie-0.9.tgz .

Use the perl script `experiments/bioinfo/parse_natalie.pl` in the 
to extract the files and then `load_natalie.m` to read the data in 
Matlab.

In the experiments/bioinfo directory, execute

    load_natalie
    [S,w,li,lj] = netalign_setup(A,B,L);
    save '../../data/natalie_graphs.mat' A B L S w li lj
    
to save the data for use with the experiments.    

Experiments
------------

* The synthetic experiments require the `gaimc` package : 
[gaimc](http://www.mathworks.com/matlabcentral/fileexchange/24134), but
it is distributed in the `experiments/misc/gaimc` directory
* For linear programming, we used Clp and a matlab interface, see 
[Clp](http://www.stanford.edu/~dgleich/notebook/2009/03/coinor_clop_for_matlab.html)
* Before running the power law experiments, please compile RandomPowerLaw.cpp
to a.out (`g++ RandomPowerLaw.cpp` ought to work)


|Experiment|Description|Figure|
|:------------------|:------------------------------------|:------------------|
|`misc/figure_2.m` | Generate figure 1 | Figure 1 |
|`powerlaw/powerlaw.m` | Synthetic power law experiments, plot with `plot_powerlaw.m` | Figure 2 |
|`evaluation/evaluate_all_problems.m` | Generate results for figure 3 plot with `evaluation_approx_plots.m` | Figure 3 |
|`numaperf/allrun.sh` | Generate speedup data, plot with `make_speedup_plots.m` and `make_speedup_plots_rameau.m` | Figures 4, 5, 6 |
|`numaperf/make_step_scaling_plots.sh` | Plot step scaling data | Figure 7 and 8 |


References
----------

* IsoRank http://isorank.csail.mit.edu
* Natalie https://www.mi.fu-berlin.de/w/LiSA/Natalie


