# README #

### What is this repository for? ###

* pyaeroopt is a Python package for scripting and integrating FRG codes, especially for optimization.
* Initially developed by Matt Zahr
* version 1.1

### How do I get set up? ###

* In your .bashrc file, define environment variables that specify the name of the FRG executables
  as they appear in your search path.  For instance, to run Sdesign on my Independence setup, I type
  `sdesign.Linux.opt`, so I'd add a line to my .bashrc that says `export SDESIGN=sdesign.Linux.opt`.
  You should do this for the environment variables SOWER, PARTMESH, XP2EXO, MATCHER, AEROF, AEROS,
  SDESIGN, and BLENDER.  If you don't set these variables, reasonable defaults are coded into PyAeroOpt.
* pyaerooopt/test/naca0012 is a good simple example of how to use pyaeroopt to script FRG codes.
  It solves a simple optimization problem using parameterized AERO-F simulations to
  evaluate the objective function.
* pyaeroopt/test/ARW2 is a very involved aeroleastic optimization problem that uses most of the 
  functionality of pyaeroopt to minimize the L/D of a full wing subject to flutter constraints.
* some stuff

### Contribution guidelines ###

* Follow PEP 8 and the rest of the Python style guide. At the very least document classes,
  function, and modules with docstrings consistent with the rest of the package.  This makes 
  it much easier to use pyaeroopt from a text editor with intellisense or similar functionality.  

