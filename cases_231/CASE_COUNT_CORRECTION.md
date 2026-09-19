# Case-count correction

The manuscript contains **257 three-vertex combinations and one two-vertex combination**, for **258 cases**. Its 230 three-vertex headings include 27 additional OR branches. This counts the listed combinations without further quotienting by geometric symmetry.

The previous manifest had 257 rows but only 256 distinct combinations: one duplicate and one omitted OR branch.

| Change | Historical identifier | Vertex combination |
|---|---|---|
| Retire duplicate; use APP-133 | APP-158 | beta delta epsilon, alpha^2 gamma, alpha beta^2 |
| Add missing OR branch | APP-177-OR1 | beta delta epsilon, alpha beta gamma, alpha delta epsilon^3 |
| Already present separately | B2V-1 | alpha delta epsilon, beta^2 gamma |

APP-158 and APP-133 had identical vertices, phase bases, and reference polynomials. APP-158's expected solution and coefficient-field norm handling are transferred to APP-133. Other APP identifiers remain unchanged; numbering is intentionally not contiguous. Reference lookup uses identifiers rather than heading positions.

For APP-177, write its three vertex rows as

```
v1 = (0,1,0,1,1)
v2 = (1,1,1,0,0)
v3 = (2,0,1,0,2)
w3 = (1,0,0,1,3)
```

The added branch has `w3 = v3 + v1 - v2`. All vertex angle sums equal 2, so this replacement gives the same linear system in both directions. Its phase basis equals APP-177's. Thus the branch shares the primary case's affine model and phase polynomial, while remaining a distinct listed vertex combination.

## Reproduction and validation

From this directory, run in Maple:

```maple
read "validate_case_manifest.mpl":
```

This checks 257 distinct records, the 230/27 split, primary reference IDs, all affine vertex and total-angle equations, the separate two-vertex record, the transferred candidate, and exact equality of the APP-177 models and phase polynomials. It does not rerun the full recursive root-of-unity search.

Both PowerShell extractors accept `Case` and `Case.`. The shared `appendix_case_ids.ps1` resolves stable IDs from the reference manifest. Unknown, missing, or duplicated identities stop extraction instead of shifting IDs or overwriting partial output. The current manuscript yields 257 three-vertex records, 230 reference-polynomial records, and one separately skipped two-vertex heading.

Two archived phase bases differ from the current manuscript (APP-048, APP-139), and six printed reference-polynomial expressions have changed (APP-002, APP-048, APP-075, APP-076, APP-119, APP-230). This count correction preserves the archived algebraic data for retained IDs. Regeneration from a newer manuscript must be followed by polynomial comparison; matching case identities does not certify matching formulas.

## Archived results

The version-2 CSVs and computation report remain historical results with their original 231-primary count and APP-158 identifiers. They are not a new run. Old APP-158 results map to APP-133. Start fresh runs in a clean output directory because raw row positions and checkpoint compatibility have changed.

The `cases_231` directory and `*_231*` filenames remain valid for compatibility. The root `FILES.sha256` is updated for changed and added files.
