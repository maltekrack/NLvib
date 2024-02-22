# NLvib
For general information, please see README file of `NLvib - Basic` branch.

# Branch `NLvib - PEACE`: Changes w.r.t. `NLvib - Basic`
- `solve_and_continue.m` now additionaly stores 
   - local tangents on converged solution points
   - solution point dependent DSCALE 
   - time spent per solution point
- computation of the local tangents (`compute_tangent.m`) and the extended residuum (`extended_residual.m`) are now separate from `solve_and_continue.m`
- `peace.m` performs *H-refinement* and *M-refinement* according to [1]
- EXAMPLE folder includes test cases according to [1]

[1] Gross, Johann; Gupta, Vasudev; Berthold, Christian; Krack, Malte (2024): A new paradigm for multi-fidelity continuation using parallel model refinement. Computer Methods in Applied Mechanics and Engineering. (https://doi.org/10.1016/j.cma.2024.116860)https://doi.org/10.1016/j.cma.2024.116860
