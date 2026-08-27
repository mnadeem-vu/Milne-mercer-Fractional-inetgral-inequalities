MATLAB CODE FOR "MILNE-MERCER INEQUALITIES FOR H-CONVEX FUNCTIONS BUILT ON
THE NON-CONFORMABLE FRACTIONAL INTEGRAL"
============================================================================

This folder contains MATLAB scripts that reproduce every numerical table
and figure in the paper. All scripts were translated from, and checked
against, independent high-precision (mpmath, 20-30 significant digit)
Python verification of the same closed-form propositions; the printed
numbers each script produces match the values quoted in the manuscript.

HOW TO RUN
----------
Keep all files in the same folder and add that folder to the MATLAB path
(or just cd into it). Each fig_*.m script is run independently:

    >> fig_exampleA

It will print a small table to the Command Window and open a figure
window with the corresponding plot. Requires only base MATLAB (no
toolboxes) except fig_bessel_vs_p.m, which uses besseli (Statistics-free,
part of base MATLAB) and fig_qdigamma_vs_p.m / fig_matrix_vs_lambda.m,
which use the Symbolic-free "integral" adaptive-quadrature function
(base MATLAB, no toolbox required).

FILE LIST
---------
Core (shared) functions -- used by the Example A-D and Theorem-comparison
figures, implementing the general h-convex framework of Lemma 3.1 and
Theorems 3.5 / 3.14 / 3.22 directly:
    LHS_exact.m   - exact left-hand side of the Milne-Mercer identity
    Kfun.m        - kernel K(t) = 4/3 - (1-t)^(1-beta)
    RHS_thm1.m    - Theorem 3.5 bound (triangle inequality)
    RHS_thm2.m    - Theorem 3.14 bound (Hoelder-type)
    RHS_thm3.m    - Theorem 3.22 bound (improved power-mean)

Figure-producing scripts (Section 4, Numerical Verification):
    fig_exampleA.m                 -> exampleA.png
    fig_exampleB.m                 -> exampleB.png
    fig_exampleC.m                 -> exampleC.png
    fig_exampleD.m                 -> exampleD.png
    fig_comparison.m               -> comparison.png
    fig_thm1_vs_thm2.m             -> thm1_vs_thm2.png
    fig_thm1_vs_thm2_vs_thm3.m     -> thm1_vs_thm2_vs_thm3.png

Figure-producing scripts (Section 5, Applications) -- these use the
closed-form classical (beta -> 0+) propositions directly rather than the
generic engine, since each special-function application has its own
explicit derivative formula:
    fig_bessel_vs_p.m              -> bessel_vs_p.png    (Section 5.1)
    fig_means_vs_p.m               -> means_vs_p.png     (Section 5.2)
    fig_qdigamma_vs_p.m            -> qdigamma_vs_p.png  (Section 5.3)
    fig_quadrature_convergence.m   -> quadrature_convergence.png (Section 5.4)
    fig_matrix_vs_lambda.m         -> matrix_vs_lambda.png       (Section 5.5)

NOTES
-----
- The x-axis in the Example A-D and comparison/theorem-comparison figures
  is labelled "alpha" in the plotted image but is mathematically the
  fractional order beta in (0,1) from the paper; the variable name in the
  scripts (beta_grid) reflects the correct mathematical meaning.
- Every script prints its computed values at the exact parameter points
  quoted in the paper's tables/remarks (e.g. beta = 0.2, 0.4, 0.6, 0.8, or
  p = 2, 3, 5, 10) so the printed numbers can be checked directly against
  the manuscript text.
- fig_bessel_vs_p.m, fig_means_vs_p.m, and fig_qdigamma_vs_p.m implement
  the classical (beta -> 0+) limit of the general framework directly
  through the closed-form coefficients derived in the paper (e.g.
  C(s) = 4(1+2^(1-s))/(3(s+1)) - 1/((s+1)(s+2)) - 2^(1-s)/(s+2)), rather
  than calling RHS_thm1/2/3 with beta forced to 0, because those
  closed forms are what the corresponding propositions state explicitly.
