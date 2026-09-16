# Generic exact elimination engine for Maple polynomial/rational systems.
# Projection and certification are kept separate on purpose.

if not assigned(a4b_solver) then
    a4b_solver := table():
end if:

a4b_solver[NormalizeSystem] := proc(eqs::list, vars::list)
    local p, q, out, denoms;
    out := []:
    denoms := []:
    for p in eqs do
        q := p:
        if type(q, equation) then
            q := lhs(q)-rhs(q):
        end if:
        q := normal(q):
        denoms := [op(denoms), denom(q)]:
        out := [op(out), expand(numer(q))]:
    end do:
    return table([equations=out, variables=vars, denominators=denoms]):
end proc:

a4b_solver[VariablesExcept] := proc(vars::list, elim)
    local out, v;
    out := []:
    for v in vars do
        if v <> elim then out := [op(out), v] end if:
    end do:
    return out
end proc:

a4b_solver[ProjectSystem] := proc(eqs::list, vars::list, keep::list)
    local nrm, polys, elim, order, basis, projected, p, v;
    nrm := a4b_solver[NormalizeSystem](eqs, vars):
    polys := nrm[equations]:
    elim := []:
    for v in vars do
        if not member(v, keep) then elim := [op(elim), v] end if
    end do:
    order := [op(elim), op(keep)]:
    if nops(polys)=0 then
        basis := []
    else
        basis := Groebner:-Basis(polys, plex(op(order))):
    end if:
    projected := []:
    for p in basis do
        if (indets(p) intersect {op(elim)}) = {} then
            projected := [op(projected), expand(p)]
        end if
    end do:
    return table([
        input_equations=polys, input_variables=vars, eliminated=elim,
        retained=keep, order=order, groebner_basis=basis,
        equations=projected, denominators=nrm[denominators]
    ]):
end proc:

a4b_solver[EliminateOne] := proc(eqs::list, vars::list, elim)
    local keep;
    if not member(elim, vars) then
        error "elimination variable is not in the variable list"
    end if:
    keep := a4b_solver[VariablesExcept](vars, elim):
    return a4b_solver[ProjectSystem](eqs, vars, keep)
end proc:

a4b_solver[EliminateSequence] := proc(eqs::list, vars::list, order::list)
    local current_eqs, current_vars, stages, elim, stage, k;
    current_eqs := eqs:
    current_vars := vars:
    stages := []:
    for k to nops(order) do
        elim := order[k]:
        stage := a4b_solver[EliminateOne](current_eqs, current_vars, elim):
        stages := [op(stages), stage]:
        current_eqs := stage[equations]:
        current_vars := stage[retained]:
    end do:
    return table([
        initial_equations=eqs, initial_variables=vars,
        elimination_order=order, stages=stages,
        final_equations=current_eqs, final_variables=current_vars
    ]):
end proc:

a4b_solver[CheckCandidate] := proc(eqs::list, vars::list, candidate::list)
    local vals, p, r, out;
    if nops(candidate)<>nops(vars) then
        error "candidate length must equal variable list length"
    end if:
    vals := zip((v,w)->v=w, vars, candidate):
    out := []:
    for p in subs(vals, eqs) do
        r := simplify(normal(p)):
        out := [op(out), r]
    end do:
    return table([
        substitutions=vals, residuals=out,
        certified=andmap(r->evalb(r=0), out)
    ]):
end proc:

a4b_solver[CheckConditions] := proc(conditions::list, vars::list, candidate::list)
    local vals, c, r, out;
    if nops(candidate)<>nops(vars) then
        error "candidate length must equal variable list length"
    end if:
    vals := zip((v,w)->v=w, vars, candidate):
    out := []:
    for c in conditions do
        r := simplify(normal(subs(vals,c))):
        out := [op(out), r]
    end do:
    return table([
        substitutions=vals, residuals=out,
        satisfied=andmap(r->evalb(r<>0), out)
    ]):
end proc:

# Convenience wrapper for the old trigonometric workflow.  The caller first
# supplies a Laurent expression/system and a monomial-coordinate change; after
# denominators are cleared, the same exact elimination engine is used.
a4b_solver[LaurentProject] := proc(eqs::list, vars::list, keep::list)
    return a4b_solver[ProjectSystem](eqs, vars, keep)
end proc:

a4b_solver[LaurentEliminateSequence] := proc(eqs::list, vars::list, order::list)
    return a4b_solver[EliminateSequence](eqs, vars, order)
end proc:

# Apply a monomial coordinate change.  Row j of M gives the exponents of the
# old variable oldvars[j] in the new variables newvars.  Rational entries are
# allowed, which covers the square-root change used in the old worksheets.
a4b_solver[MonomialSubstitute] := proc(expr::anything,
        oldvars::list, newvars::list, M::Matrix)
    local rules, j, i, mon;
    if LinearAlgebra:-RowDimension(M) <> nops(oldvars) or
       LinearAlgebra:-ColumnDimension(M) <> nops(newvars) then
        error "monomial matrix dimensions do not match variable lists"
    end if:
    rules := []:
    for j to nops(oldvars) do
        mon := 1:
        for i to nops(newvars) do
            mon := mon*newvars[i]^M[j,i]
        end do:
        rules := [op(rules), oldvars[j]=mon]
    end do:
    return simplify(subs(rules,expr))
end proc:

a4b_solver[PrintElimination] := proc(result::table)
    local k, st;
    printf("variables: %a\n", result[initial_variables]):
    printf("order: %a\n", result[elimination_order]):
    for k to nops(result[stages]) do
        st := result[stages][k]:
        printf("step %d: eliminate %a; retain %a\n", k, st[eliminated], st[retained]):
        printf("  projected equations: %a\n", st[equations]):
    end do:
    printf("final equations: %a\n", result[final_equations]):
    printf("final variables: %a\n", result[final_variables]):
    return NULL
end proc:

NULL:
