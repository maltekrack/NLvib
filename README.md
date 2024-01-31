# NLvib
 NLvib is a Matlab tool for nonlinear vibration analysis.

An overview of its capabilities, included examles, the monograph on Harmonic Balance, and some further presentation material can be found on https://www.ila.uni-stuttgart.de/nlvib/.

Please also see the manual in the SRC folder.

The tool should work well with a wide range of Matlab releases. It mainly relies on the optimization toolbox.
To get it to run under OCTAVE, you need to change the line(s) that thes the solver options to something like:
   Solopt = optimset(optimset ("fsolve"),'Display','off',... 'Jacobian','on','MaxIter',50);
Also, some 'legend' calls (for figures) might require modification.

# Changes w.r.t. NLvib - BASIC
- `solve_and_continue.m` now additionaly stores 
   - local tangents on converged solution points
   - solution point dependent DSCALE 
   - time spent per solution point
- computation of the local tangents (`compute_tangent.m`) and the extended residuum (`extended_residual.m`) function are now excluded from `solve_and_continue.m`
- `peace.m` performs *H-refinement* and *M-refinement* according to **paper-reference-to-be-added**
- EXAMPLE folder includes test cases according to the publication **paper-reference-to-be-added**
