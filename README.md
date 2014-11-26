DpBayesianNet
=============

C liblary of the Bayesian network including learning algorithm which was optimized by dynamic programming of Silander and Myllymaki(2006).   

## Features
* Improve parformance


## Install
You require cmake and GSL

**[Case of Linux]**  
1. Input this following on shell.  

    $ cmake .     #(Be careful dot.)
    $ make && make install

**[Case of MinGW on Windows]**  
1. Input this following on command prompt.

    $ PATH=C:\MinGW\bin;%PATH%
    $ cmake -G "MinGW Makefiles"
    (input one more)
    $ cmake -G "MinGW Makefiles"
    $  mingw32-make

2. Add "libdpgmm.a" and "dpgmm.h" , "laa.h" to your project.

## Usage
Read example.c
