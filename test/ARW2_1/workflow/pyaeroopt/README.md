# README #

### What is this repository for? ###

* pyaeroopt is a Python package for scripting and integrating FRG codes, especially for optimization.
* Initially developed by Matt Zahr
* version 1.1

### How do I get set up? ###

* Clone this repository onto your machine
* [Install an anaconda3 distribution](https://www.continuum.io/downloads) onto your machine
* Copy the pyaeroopt34 conda environment from my Independence directory.  This has a bunch of 
  python packages already compiled for Independence, and some proprietary optimizers that I
  can't put on this repository.
* Before using pyaeroopt, be sure to run `source activate pyaeroopt34` to switch into this
  conda environment.  I called the environment pyaeroopt34 specifically to remind you that
  if you launch python and see a version other than 3.4, you are not in the right environment.
* In your .bashrc file, define environment variables that simply specify the name of the FRG executables
  as they appear in your search path.  For instance, to run Sdesign on my Independence setup, I type
  `sdesign.Linux.opt`, so I'd add a line to my .bashrc that says `export SDESIGN=sdesign.Linux.opt`.
  You should do this for the environment variables SOWER, PARTMESH, XP2EXO, MATCHER, AEROF, AEROS,
  SDESIGN, and BLENDER.  If you don't set these variables, reasonable defaults are coded into PyAeroOpt.
* pyaerooopt/test/naca0012 is a good simple example of how to use pyaeroopt to script FRG codes.
  It runs a constant-reconstruction flow simulation followed by a linear-reconstruction simulation
  automatically, and allows you to change some flow angles.
  pyaeroopt/test/ARW2 is a very involved aeroleastic optimization problem that uses most of the 
  functionality of pyaeroopt to minimize the L/D of a full wing subject to flutter constraints.

### Contribution guidelines ###

* Follow PEP 8 and the rest of the Python style guide! At the very least document classes,
  function, and modules with docstrings consistent with the rest of the package.  This makes 
  it much easier to use pyaeroopt from a text editor with intellisense or similar functionality.  
* Remember to test your code from the pyaeroopt34 conda environment (if you type `python` at the 
  command line and see a python version other than 3.4, you're not in the right environment)
* at the top level directory of pyaeroopt, run `./clean` immediately before `hg add` and `hg commit`
  to remove \*.pyc files and \_\_pycache\_\_ directories. 

### Who do I talk to? ###
* Spenser Anderson is the main point of contact (aspenser at stanford dot edu)
