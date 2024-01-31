# NLvib
NLvib is a Matlab tool for nonlinear vibration analysis.

An overview of its capabilities, included examles, the monograph on Harmonic Balance, and some further presentation material can be found on https://www.ila.uni-stuttgart.de/nlvib/.

Please also see the manual in the SRC folder.

The tool should work well with a wide range of Matlab releases. It mainly relies on the optimization toolbox.
To get it to run under OCTAVE, you need to change the line(s) that thes the solver options to something like:
   Solopt = optimset(optimset ("fsolve"),'Display','off',... 'Jacobian','on','MaxIter',50);
Also, some 'legend' calls (for figures) might require modification.

# NLvib - Basic and NLvib - PEACE 
The main branch of NLvib is now split in two options:
1. NLvib - Basic: traditional, compact implementation
2. NLvin - PEACE: adapted implementation for model refinement capabilities

All EXAMPLES included in `NLvib - Basic` can be run in `NLvib - PEACE`, but not vice versa.
