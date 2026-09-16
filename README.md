# Spherical Tilings: Maple Computations

Maple programs, algebraic data, and current computation reports for edge-to-edge spherical tilings by almost equilateral pentagons with edge lengths $a^4b$.

The accompanying constructions and figures are in [spherical-tilings-geogebra](https://github.com/m1i2p3u4/spherical-tilings-geogebra).

## Contents

| Directory | Contents |
|---|---|
| [cases_231](cases_231) | Affine angle models, Laurent polynomials, recursive root-of-unity elimination, and appendix comparisons for 231 primary cases |
| [three_variable_demo](three_variable_demo) | A three-variable recursive-elimination example |
| [four_variable_experiment](four_variable_experiment) | Four-variable extensions, intermediate results, and an experimental status report |
| [rational_tetrahedra](rational_tetrahedra) | A separate Python prototype for rational dihedral angles of tetrahedra |

The `.mpl` files contain Maple source. The `.m` files are Maple saved data and should be read with Maple. CSV files contain comparison results. Intermediate algebraic data are included to retain the archived computation state.

## Running the Maple programs

Use a local Maple installation. Set the working directory to the directory containing the selected entry point; internal file references are relative to that directory. For example, after changing into `cases_231`, run:

```maple
read "run_recursive_231_solution_audit.mpl":
```

The full audit can be expensive. For a smaller run, create `cases_231/solution_audit_range.mpl` with:

```maple
first_case := 1:
last_case := 1:
```

Remove that range file before a full run. Generated reports may overwrite reports in the working directory; use a separate checkout for new experiments. The PowerShell extraction scripts take the manuscript source as an explicit argument; the manuscript itself is not included.

The demo entry point is `three_variable_demo/run_recursive_3variable_demo.mpl`. The four-variable experiments are described in [EXPERIMENT_STATUS.md](four_variable_experiment/EXPERIMENT_STATUS.md).

## Python prototype

From the repository root, with Python 3.10 or later:

```sh
python -m rational_tetrahedra audit
python -m rational_tetrahedra search 12
python -m unittest discover -s rational_tetrahedra/tests
```

The prototype uses the Python standard library.

## Scope of the results

The current three-variable summary is [current_formula_solution_audit_v2.md](cases_231/current_formula_solution_audit_v2.md). Superseded comparison reports and obsolete result tables are omitted. Intermediate Maple data required by the experiments remain available.

The four-variable and tetrahedron programs are experimental. Their status reports distinguish verified candidates and families from unresolved completeness questions.
