# Rational-angle spherical pentagons: Maple computations

Programs and results for the 257 three-vertex combinations and the separate two-vertex case in the accompanying paper on spherical tilings by congruent pentagons with edge combination `a^4 b`.

The appendix contains 230 primary three-vertex headings and 27 alternative (`OR`) combinations: **257 + 1 = 258 distinct cases**. These are case records processed by shared Maple programs; the number of `.mpl` files is not the case count. Case identifiers are determined by the vertex combinations. APP-133 covers the repeated heading APP-158, and APP-177-OR1 is included.

## Contents

- `cases/`: Maple programs, manuscript inputs, and extraction scripts.
- `results/`: outputs, logs, input hashes, and execution metadata from the run described in [RUN_REPORT.md](RUN_REPORT.md).
- `run_all.ps1`: execution and verification from empty working directories.
- `computation-record.json` and `FILES.sha256`: reproducibility metadata and file checksums.

All angles in the data files are measured in units of pi. Vertex vectors use the order `(alpha, beta, gamma, delta, epsilon)`. Phase vectors additionally use `q=1/f`.

## Run

The included driver requires PowerShell 7 on Windows and a licensed Maple command-line installation. This release was run with Maple 2024.2.

```powershell
pwsh -File ./run_all.ps1 -Maple 'C:/Program Files/Maple 2024/bin.X86_64_WINDOWS/cmaple.exe' -OutputDirectory ./results-new
```

Choose a directory that does not exist. The driver creates a separate working directory, runs all 257 three-vertex records, checks the 230 printed primary polynomials, verifies the parameter families and finite candidates, and calculates the two-vertex resultants and candidate checks. It stops if a required comparison fails. Four independent Maple processes are used by default, with one Maple CPU thread per search process; use `-Workers` to change the process count.

To check the inputs against a manuscript source as well:

```powershell
pwsh -File ./run_all.ps1 -OutputDirectory ./results-new -ManuscriptPath ./a4br.tex
```

The manuscript source is not included. Its SHA-256 is recorded in `results/manuscript.sha256`. With `-ManuscriptPath`, the driver checks extraction, compares the three input datasets, and refuses to run if they differ.

## Euler bound

The search applies the identity

`f = 12 + 2 * sum((k-3) * v_k, k >= 4)`.

Each required vertex type occurs at least once. For required types of degrees `d_1,d_2,d_3`, the program therefore requires

`f >= 12 + 2 * ((d_1-3) + (d_2-3) + (d_3-3))`.

In particular, any required vertex of degree at least four gives `f >= 14`. Cases with only degree-three vertices may still have `f=12`. Algebraic candidates rejected by this bound are recorded separately in `results/degree_bound_exclusions.csv`.

## Verification scope

The computations check the manuscript's listed angle candidates and symbolic equations. Finite candidates undergo exact vertex, angle-sum, and cyclotomic checks of equation (2.6), after the search's numerical screening. The two printed parameter families are checked symbolically.

Factor logs retain `rank1_family` and `one_variable_factor` entries. The general classification of all positive-dimensional components and the geometric and combinatorial arguments for tiling existence remain part of the paper. The two-vertex program recalculates the 15 resultants and checks the four listed candidates; it does not independently enumerate all three-variable torsion points. See the run report for the precise results.
