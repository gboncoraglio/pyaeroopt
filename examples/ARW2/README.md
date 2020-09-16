# ARW2 Problem README #

The purpose of this folder (pyaeroopt/test/ARW2) is to provide a well-documented example of how to use 
PyAeroOpt to solve an optimization problem using FRG codes.   

The layout of the code in this example is as follows. 

* pyao              : defines classes for running codes
	* aerof_blk.py  : defines input blocks for the AERO-F AeroelasticSensitivity problem
	* aerof_cls.py  : defines a class to run the AERO-F AeroelasticSensitivity problem
	* aeros_blk.py  : defines input blocks for the AERO-S AeroelasticSensitivity problem
	* aeros_cls.py  : defines a class to run the AERO-S AeroelasticSensitivity problem
	* sdesign_blk.py: defines a input blocks do the shape deformation and material changes
	* sdesign_cls.py: defines a class to do shape deformation and material changes
	* ARW2opt.py    : defines a class to scripts the codes for evaluating objectives and constraints
* workflow          : code to actually run the problem
 	* __init__.py   : define compute objects and other variables for the problem
 	* clean.py		: cleans files produced by running the FRG codes
 	* makebin.py    : create FRG binaries from meshes
 	* workflow.py   : run the optimization problem
 	* binaries/     : hold FRG binaries produced by makebin.py
 	* Sresults/     : holds AERO-S outputs
 	* AeroF-Files/  : holds mesh topology
 	* AeroS-Files/  : holds files defining AERO-S model
 	* Fdata/        : holds intermediate AERO-F data
 	* Fresults/     : holds intermediate AERO-F data and restart data
 	* Sdesign/      : holds Sdesign files
 	* optExamples/  : a few examples of using pyoptsparse for more general problems
 	* out/          : holds the output of the latest objective/constraint evaluation


The basic organizational philosophy in this example is that the workflow/ folder holds the scripts
that a person simply running the optimization problem would actually run, like the makebin.py script
to create the FRG binaries, the clean.py file to clean up the folder if they want to return to the 
folder's original state, and the workflow.py script to run the problem.  The results of running the
codes are all placed in this folder as well.

The classes that do all of the lower-level work are in the pyao/ folder.  There are classes defined
here for each type of problem run with AERO-F, AERO-S, or Sdesign.  