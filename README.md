DpBayesianNet
=============

C liblary of the Bayesian network including learning algorithm which was optimized by dynamic programming of Silander and Myllymaki(2006).   

## Features
* Both learning and inference was implemented on exact algorithms
* Cross platform
* Scoring:BIC
* Maximum number of node:16
* Source code simply written.
* License: New BSD  


## Install
DpBayesianNet uses the following dependencies:
cmake

**[Case of Linux]**  
1. Input this following on shell.  

    $ cmake .     #(Be careful dot.)
    $ make && sudo make install

**[Case of MinGW on Windows]**  
[1] Input this following on command prompt.

    $ PATH=C:\MinGW\bin;%PATH%
    $ cmake -G "MinGW Makefiles"
    (input one more)
    $ cmake -G "MinGW Makefiles"
    $  mingw32-make

[2] Add "libdpBayesianNet.a" and "dpBayesianNetwork.h" , "laa.h" to your project.

## Usage
Read dpBayesianNetwork.h and example.c
