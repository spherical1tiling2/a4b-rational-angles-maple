# A4BR recursive root-of-unity elimination

## Main files

- `canonical_phase_engine.mpl`: converts the affine angle model to integral-exponent canonical phase variables. It handles half-integer printed exponents by using `X=exp(i*pi*u1/D1)` and `Y=exp(i*pi*u2/D2)`.
- `recursive_torsion_engine.mpl`: factor-by-factor sign/square transformations, resultants, cyclotomic coefficient norm, gcd diagnostics, and recursive 2D/3D reduction.
- `phase_backsubstitution.mpl`: restores phase lifts and the affine angle variables after roots of unity have been found.
- `run_recursive_231.mpl`: runs the recursive tree over the 230 primary cases in the corrected manifest, skipping its 27 OR branches. Its historical filename is retained. It writes `recursive_231_tree.log` and `recursive_231_phase_maps.log`.
- `run_recursive_3variable_demo.mpl`: demonstrates the 3-variable route: 15 branches from 3D to 2D, then 7 branches per 2D factor. After a candidate `x` (hence a candidate `f`) is found, set `_candidate_x` and run the remaining 2-variable tree in `(u,w)`; these remain two independent angle variables.

## Mathematical order

For a 2-variable factor, the engine first factors it, then computes the seven nontrivial sign/square resultants. Each nonzero resultant is factored again and sent to the one-variable cyclotomic gcd test.

For a 3-variable factor, it computes fifteen 3-to-2 resultants. Every nonzero 2-variable factor is then sent through its own seven 2-to-1 resultants. The nominal count is therefore (15\times7=105), before duplicate-factor removal.

The coefficient-field norm is applied after each resultant rather than to the original high-degree polynomial. This preserves original roots while avoiding unnecessary degree inflation. All candidates must still be back-substituted into the original factor, because norms and resultants can add conjugate or projection-only solutions.

When a root determines (f), (q=1/f) is substituted only at the back-substitution stage. In the three-variable example, fixing (f) leaves the two phase variables (y,z) (or (u,w)) as two angle variables; they are not collapsed to one variable merely because (f) is fixed.
