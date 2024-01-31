# NLvib
For general information, please see README file of NLvib-Basic branch.

# Changes w.r.t. NLvib - Basic
- `solve_and_continue.m` now additionaly stores 
   - local tangents on converged solution points
   - solution point dependent DSCALE 
   - time spent per solution point
- computation of the local tangents (`compute_tangent.m`) and the extended residuum (`extended_residual.m`) are now separate from `solve_and_continue.m`
- `peace.m` performs *H-refinement* and *M-refinement* according to **paper-reference-to-be-added**
- EXAMPLE folder includes test cases according to **paper-reference-to-be-added**
