# Exact phase back-substitution.
# If X=exp(i*pi*u1/D1)=zeta_N^r and
# Y=exp(i*pi*u2/D2)=zeta_M^s, then
#   u1 = D1*(2*r/N + 2*k1),
#   u2 = D2*(2*s/M + 2*k2).
# The integer branch variables are retained until angle inequalities select
# the admissible lifts.  This is also the mechanism used when f is fixed:
# the remaining phase rows still represent two independent angles.

`a4b_solver/CanonicalRootBackSubstitute` := proc(T,Q,N1,r1,N2,r2)
    local u1,u2,k1,k2,aff,sol,vals,j;
    u1 := '_u1': u2 := '_u2': k1 := '_k1': k2 := '_k2':
    aff := T["affine_expressions"]:
    vals := [u1=Q["D1"]*(2*r1/N1+2*k1),
             u2=Q["D2"]*(2*r2/N2+2*k2)]:
    sol := [seq(normal(subs(vals,aff[j])),j=1..6)]:
    return table(["affine_angles"=sol,"branch_variables"=[k1,k2],
                  "root_data"=[[N1,r1],[N2,r2]]]):
end proc:

`a4b_solver/PhaseCongruenceSystem` := proc(affine,phase_rows,orders,residues)
    local k,j,r,rhs,lhs,eqs,branchvars,sol;
    if nops(phase_rows)<>nops(orders) or nops(orders)<>nops(residues) then
        error "phase rows, orders, and residues must have equal length"
    end if:
    eqs := []: branchvars := []:
    for j from 1 to nops(phase_rows) do
        k := cat('_phase_branch_',j):
        branchvars := [op(branchvars),k]:
        lhs := add(phase_rows[j][r]*affine[r],r=1..6):
        rhs := 2*residues[j]/orders[j]+2*k:
        eqs := [op(eqs),lhs=rhs]:
    end do:
    sol := solve(eqs,['_u1','_u2']):
    return table(["equations"=eqs,"branch_variables"=branchvars,
                  "solution_parameters"=sol]):
end proc:

NULL:
