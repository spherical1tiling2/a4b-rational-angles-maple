if not assigned(a4b_solver) then a4b_solver := table() end if:

`a4b_solver/CanonicalScalarPolynomial` := proc(P,vars::list)
    local G,terms,term,best,bestexp,exps,j,better,c,mon;
    G := `a4b_solver/PrimitivePolynomial`(P,vars):
    if G=0 then return 0 end if:
    terms := [op(expand(G))]:
    best := terms[1]:
    bestexp := [seq(degree(best,vars[j]),j=1..nops(vars))]:
    if nops(terms)>1 then
        for term in terms[2..-1] do
            exps := [seq(degree(term,vars[j]),j=1..nops(vars))]:
            better := false:
            for j from 1 to nops(vars) do
                if exps[j]>bestexp[j] then better:=true: break
                elif exps[j]<bestexp[j] then break
                end if:
            end do:
            if better then best:=term: bestexp:=exps end if:
        end do:
    end if:
    mon := mul(vars[j]^bestexp[j],j=1..nops(vars)):
    c := normal(best/mon):
    return expand(G/c):
end proc:

`a4b_solver/FindProjectivePolynomial` := proc(P,L::list,vars::list)
    local C,j;
    C := `a4b_solver/CanonicalScalarPolynomial`(P,vars):
    for j from 1 to nops(L) do
        if evalb(expand(C-L[j])=0) then return j end if:
    end do:
    return 0:
end proc:

NULL:
