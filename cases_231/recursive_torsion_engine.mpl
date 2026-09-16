# Recursive root-of-unity elimination engine.
#
# The engine follows the reduction tree used in the paper:
#   n variables -> 2^(n+1)-1 resultants in n-1 variables,
#   then recursively repeats.  Resultants give necessary conditions;
#   every candidate is subsequently back-substituted into the preceding
#   polynomial/factor.

with(numtheory):
with(PolynomialTools):

`a4b_solver/VariablesExcept` := proc(vars::list,elim)
    local out,v;
    out := []:
    for v in vars do if v<>elim then out := [op(out),v] end if end do:
    return out:
end proc:

`a4b_solver/SignSquareTransforms` := proc(F, vars::list)
    local n,mask,i,sgn,rules,out,id;
    n := nops(vars): out := []:
    # Sign changes, excluding the identity.
    for mask from 1 to 2^n-1 do
        rules := []: id := "sign":
        for i from 1 to n do
            sgn := if irem(iquo(mask,2^(i-1)),2)=1 then -1 else 1 end if:
            rules := [op(rules),vars[i]=sgn*vars[i]]:
            id := cat(id,if sgn=-1 then "-" else "+" end if):
        end do:
        out := [op(out),table(["id"=id,"expression"=expand(subs(rules,F))])]:
    end do:
    # Squared variables, including all sign choices.
    for mask from 0 to 2^n-1 do
        rules := []: id := "square":
        for i from 1 to n do
            sgn := if irem(iquo(mask,2^(i-1)),2)=1 then -1 else 1 end if:
            rules := [op(rules),vars[i]=sgn*vars[i]^2]:
            id := cat(id,if sgn=-1 then "-" else "+" end if):
        end do:
        out := [op(out),table(["id"=id,"expression"=expand(subs(rules,F))])]:
    end do:
    return out:
end proc:

# The seven auxiliary curves in Beukers--Smyth depend on the coefficient
# field.  For coefficients in Q(zeta_N), the square branch is not always
# F(+-x^2,+-y^2) with the same coefficients.  If N is odd, the required
# branch is F^sigma(+-x^2,+-y^2), sigma(zeta_N)=zeta_N^2.  If 4|N, the
# required branch is F^tau(+-x,+-y), tau(zeta_N)=-zeta_N.  The three
# nonidentity sign branches with the original coefficients are retained in
# both cases.  For N=1 or 2 the coefficient field is Q and the usual square
# branches apply.
`a4b_solver/CyclotomicSignPowerTransforms` := proc(F, vars::list, N, z)
    local n,mask,i,sgn,rules,out,id,H,coeffmap;
    n:=nops(vars): out:=[]:
    for mask from 1 to 2^n-1 do
        rules:=[]: id:="sign":
        for i from 1 to n do
            sgn:=if irem(iquo(mask,2^(i-1)),2)=1 then -1 else 1 end if:
            rules:=[op(rules),vars[i]=sgn*vars[i]]:
            id:=cat(id,if sgn=-1 then "-" else "+" end if):
        end do:
        out:=[op(out),table(["id"=id,"expression"=expand(subs(rules,F))])]:
    end do:
    if N<>1 and irem(N,4)=0 then
        # tau(zeta_N)=zeta_N^(1+N/2)=-zeta_N.
        coeffmap:=subs(z=-z,F):
        for mask from 0 to 2^n-1 do
            rules:=[]: id:="tau":
            for i from 1 to n do
                sgn:=if irem(iquo(mask,2^(i-1)),2)=1 then -1 else 1 end if:
                rules:=[op(rules),vars[i]=sgn*vars[i]]:
                id:=cat(id,if sgn=-1 then "-" else "+" end if):
            end do:
            out:=[op(out),table(["id"=id,"expression"=expand(subs(rules,coeffmap))])]:
        end do:
    else
        for mask from 0 to 2^n-1 do
            rules:=[]: id:="square":
            for i from 1 to n do
                sgn:=if irem(iquo(mask,2^(i-1)),2)=1 then -1 else 1 end if:
                rules:=[op(rules),vars[i]=sgn*vars[i]^2]:
                id:=cat(id,if sgn=-1 then "-" else "+" end if):
            end do:
            if N<>1 and irem(N,2)=1 then
                H:=subs(z=z^2,F):
            else
                H:=F:
            end if:
            out:=[op(out),table(["id"=id,"expression"=expand(subs(rules,H))])]:
        end do:
    end if:
    return out:
end proc:

# Compare two polynomials up to a nonzero scalar in Q(zeta_N).  The scalar
# freedom is essential because an irreducible factor is only defined up to a
# coefficient-field unit.
`a4b_solver/CyclotomicProjectiveEqual` := proc(F,G,vars::list,N,z)
    local terms,t,c,cg,v,lam,rootvar,exps,j;
    if N=1 then return evalb(expand(G-F)=0) end if:
    terms:=[op(expand(F))]:
    if nops(terms)=0 then return evalb(G=0) end if:
    t:=terms[1]: exps:=[seq(degree(t,vars[j]),j=1..nops(vars))]:
    c:=F: cg:=G:
    for j from 1 to nops(vars) do
        c:=coeff(c,vars[j],exps[j]): cg:=coeff(cg,vars[j],exps[j])
    end do:
    # The monomial support is fixed, so the ratio of corresponding first
    # coefficients determines the projective scalar.
    if c=0 then return false end if:
    lam:=normal(cg/c):
    rootvar:='_a4b_projective_root':
    return evalb(evala(normal(subs(z=RootOf(cyclotomic(N,rootvar)),
                                  expand(G-lam*F))))=0)
end proc:

# Find the smallest cyclotomic conductor M|N for which the projective
# coefficient field of F is contained in Q(zeta_M).  Testing the subgroup
# a=1 (mod M) is exact and also handles factors whose displayed zeta_N
# notation uses a nonminimal ambient field.
`a4b_solver/CyclotomicProjectiveConductor` := proc(F,vars::list,N,z)
    local M,a,ok,Fa;
    if N=1 then return 1 end if:
    for M from 1 to N do
        if irem(N,M)=0 then
            ok:=true:
            for a from 1 to N do
                if igcd(a,N)=1 and (M=1 or irem(a,M)=1) then
                    Fa:=subs(z=z^a,F):
                    if not `a4b_solver/CyclotomicProjectiveEqual`(F,Fa,vars,N,z) then
                        ok:=false:
                    end if:
                end if:
            end do:
            if ok then return M end if:
        end if:
    end do:
    return N:
end proc:

`a4b_solver/LiftCyclotomicAutomorphism` := proc(N,M,target)
    local a;
    if M=1 then return 1 end if:
    for a from 1 to N do
        if igcd(a,N)=1 and irem(a-target,M)=0 then return a end if
    end do:
    error "cannot lift the requested cyclotomic automorphism":
end proc:

# Field-aware version of the seven auxiliary curves.  It first detects the
# actual conductor of the factor rather than blindly using the ambient N.
`a4b_solver/CyclotomicMinimalTransforms` := proc(F,vars::list,N,z)
    local M,n,mask,i,sgn,rules,out,id,H,a;
    M:=`a4b_solver/CyclotomicProjectiveConductor`(F,vars,N,z):
    n:=nops(vars): out:=[]:
    for mask from 1 to 2^n-1 do
        rules:=[]: id:="sign":
        for i from 1 to n do
            sgn:=if irem(iquo(mask,2^(i-1)),2)=1 then -1 else 1 end if:
            rules:=[op(rules),vars[i]=sgn*vars[i]]:
            id:=cat(id,if sgn=-1 then "-" else "+" end if):
        end do:
        out:=[op(out),table(["id"=id,"expression"=expand(subs(rules,F))])]:
    end do:
    if irem(M,2)=1 then
        if M=1 then a:=1 else a:=`a4b_solver/LiftCyclotomicAutomorphism`(N,M,2) end if:
        H:=subs(z=z^a,F):
        for mask from 0 to 2^n-1 do
            rules:=[]: id:="sigma2":
            for i from 1 to n do
                sgn:=if irem(iquo(mask,2^(i-1)),2)=1 then -1 else 1 end if:
                rules:=[op(rules),vars[i]=sgn*vars[i]^2]:
                id:=cat(id,if sgn=-1 then "-" else "+" end if):
            end do:
            out:=[op(out),table(["id"=id,"expression"=expand(subs(rules,H))])]:
        end do:
    elif irem(M,4)=0 then
        a:=`a4b_solver/LiftCyclotomicAutomorphism`(N,M,1+iquo(M,2)):
        H:=subs(z=z^a,F):
        for mask from 0 to 2^n-1 do
            rules:=[]: id:="tau":
            for i from 1 to n do
                sgn:=if irem(iquo(mask,2^(i-1)),2)=1 then -1 else 1 end if:
                rules:=[op(rules),vars[i]=sgn*vars[i]]:
                id:=cat(id,if sgn=-1 then "-" else "+" end if):
            end do:
            out:=[op(out),table(["id"=id,"expression"=expand(subs(rules,H))])]:
        end do:
    else
        H:=F:
        for mask from 0 to 2^n-1 do
            rules:=[]: id:="square":
            for i from 1 to n do
                sgn:=if irem(iquo(mask,2^(i-1)),2)=1 then -1 else 1 end if:
                rules:=[op(rules),vars[i]=sgn*vars[i]^2]:
                id:=cat(id,if sgn=-1 then "-" else "+" end if):
            end do:
            out:=[op(out),table(["id"=id,"expression"=expand(subs(rules,H))])]:
        end do:
    end if:
    return out:
end proc:

`a4b_solver/PrimitivePolynomial` := proc(F,vars::list)
    local G,den,i,mon,terms,t,dx,mins,content;
    G := normal(F):
    den := denom(G):
    G := expand(numer(G)):
    if G=0 then return 0 end if:
    # Remove monomial factors; all torus variables are nonzero.
    terms := [op(expand(G))]:
    mins := []:
    for i from 1 to nops(vars) do
        mins := [op(mins),min(seq(degree(t,vars[i]),t=terms))]:
    end do:
    mon := 1:
    for i from 1 to nops(vars) do mon := mon*vars[i]^mins[i] end do:
    G := expand(G/mon):
    content := 1:
    return G:
end proc:

`a4b_solver/ReduceZetaPowers` := proc(F,N,z)
    local terms,t,minz,G;
    if N=1 then return expand(F) end if:
    terms := [op(expand(F))]:
    minz := min(seq(degree(t,z),t=terms)):
    G := expand(z^(-minz)*F):
    return expand(rem(G,z^N-1,z)):
end proc:

`a4b_solver/CyclotomicNorm` := proc(F,N,z)
    local terms,t,minz,H,Phi,phiN,unit;
    if N=1 then return expand(F) end if:
    Phi := cyclotomic(N,z):
    phiN := degree(Phi,z):
    # Clear a possible Laurent power in z.  This multiplies the norm by
    # the explicit unit Norm(z)^(-minz); restore that harmless sign below.
    terms := [op(expand(F))]:
    minz := min(seq(degree(t,z),t=terms)):
    H := expand(z^(-minz)*F):
    H := expand(rem(H,z^N-1,z)):
    H := expand(rem(H,Phi,z)):
    # For monic Phi_N, Norm(z)=(-1)^phi(N).  The returned polynomial is
    # therefore the exact norm up to the restored unit, with no root loss.
    unit := (-1)^(phiN*minz):
    return expand(unit*resultant(Phi,H,z)):
end proc:

`a4b_solver/FindCyclotomicOrders` := proc(s,var)
    local G,fl,entry,H,d,n,phi,out;
    G := expand(numer(normal(s))):
    if G=0 then return [] end if:
    # Factor over Q first.  Every irreducible cyclotomic factor H has
    # degree phi(n), so invphi(degree(H)) gives the complete finite list of
    # possible orders; no arbitrary scan bound is used.
    fl := factors(G)[2]: out := []:
    for entry in fl do
        H := entry[1]: d := degree(H,var):
        if d>0 then
            for n in invphi(d) do
                phi := cyclotomic(n,var):
                if rem(H,phi,var)=0 and not member(n,out) then
                    out := [op(out),n]
                end if
            end do:
        end if:
    end do:
    return out:
end proc:

`a4b_solver/SquareLiftCyclotomicOrders` := proc(baseOrders,var)
    local m,P,fl,entry,H,d,n,phi,out;
    out := []:
    for m in baseOrders do
        P := expand(cyclotomic(m,var^2)):
        fl := factors(P)[2]:
        for entry in fl do
            H := entry[1]: d := degree(H,var):
            if d>0 then
                for n in invphi(d) do
                    phi := cyclotomic(n,var):
                    if rem(H,phi,var)=0 and not member(n,out) then
                        out := [op(out),n]
                    end if
                end do:
            end if:
        end do:
    end do:
    return out:
end proc:

`a4b_solver/UnitRootDiagnostics` := proc(s,var)
    local G,Gm,G2,terms,t,allEven,n,orders,U,E,d,baseOrders;
    G := `a4b_solver/PrimitivePolynomial`(s,[var]):
    if G=0 then return table(["status"="zero"]) end if:
    Gm := gcd(G,expand(subs(var=-var,G))):
    G2 := gcd(G,expand(subs(var=var^2,G)*subs(var=-var^2,G))):
    terms := [op(expand(G))]: allEven := true:
    for t in terms do
        if irem(degree(t,var),2)<>0 then allEven := false end if:
    end do:
    if degree(Gm,var)=0 and degree(G2,var)=0 then
        orders := []:
    elif allEven then
        # G(var)=E(var^2).  Work with E, whose degree is half as large,
        # then lift each cyclotomic order through Phi_m(var^2).
        U := '_a4b_even_variable':
        d := degree(G,var):
        E := expand(add(coeff(G,var,2*n)*U^n,n=0..iquo(d,2))):
        baseOrders := `a4b_solver/FindCyclotomicOrders`(E,U):
        orders := `a4b_solver/SquareLiftCyclotomicOrders`(baseOrders,var):
    else
        orders := `a4b_solver/FindCyclotomicOrders`(G,var):
    end if:
    return table(["status"="ok","polynomial"=G,
                  "gcd_sign"=Gm,"gcd_square"=G2,
                  "all_exponents_even"=allEven,
                  "cyclotomic_orders"=orders,
                  "no_root_test"=evalb(Gm=1 and G2=1)]):
end proc:

`a4b_solver/ResultantStep` := proc(F,vars::list,elim)
    local G,others,keep,trans,h,R,ret,v;
    G := `a4b_solver/PrimitivePolynomial`(F,vars):
    others := `a4b_solver/SignSquareTransforms`(G,vars):
    keep := [seq(v,v in vars)]:
    ret := []:
    for h in others do
        R := resultant(G,h["expression"],elim):
        R := `a4b_solver/PrimitivePolynomial`(R,`a4b_solver/VariablesExcept`(vars,elim)):
        ret := [op(ret),table(["id"=h["id"],"transform"=h["expression"],"resultant"=R])]:
    end do:
    return table(["input"=G,"variables"=vars,"eliminated"=elim,
                  "retained"=`a4b_solver/VariablesExcept`(vars,elim),
                  "branches"=ret,"branch_count"=nops(ret)]):
end proc:

# Exact resultant over the cyclotomic coefficient field.
# The formal symbol z is replaced by a primitive N-th root satisfying
# Phi_N(z)=0 before elimination.  The resultant is then mapped back to the
# paper's formal zeta notation.  This is an exact coefficient-field
# homomorphism, not a floating-point approximation; multiplying by a
# nonzero algebraic scalar does not change the zero set.
`a4b_solver/CyclotomicResultant` := proc(A,B,elim,N,z)
    local t,zr,Aalg,Balg,Ralg,R;
    if N=1 then return resultant(A,B,elim) end if:
    t := '_a4b_cyclotomic_root':
    zr := RootOf(cyclotomic(N,t)):
    Aalg := subs(z=zr,A):
    Balg := subs(z=zr,B):
    Ralg := Algebraic:-Resultant(Aalg,Balg,elim,makeindependent=false):
    R := subs(zr=z,Ralg):
    return expand(R):
end proc:

`a4b_solver/EvenPartInVariable` := proc(F,var,u)
    local terms,t,d,j,allEven,H;
    H := expand(F):
    terms := [op(expand(F))]:
    if nops(terms)=0 then return table(["all_even"=true,"polynomial"=0]) end if:
    allEven := true:
    for t in terms do
        if irem(degree(t,var),2)<>0 then allEven:=false end if
    end do:
    if not allEven then return table(["all_even"=false,"polynomial"=F]) end if:
    d := degree(F,var): H := expand(add(coeff(F,var,2*j)*u^j,j=0..iquo(d,2))):
    return table(["all_even"=true,"polynomial"=H])
end proc:

`a4b_solver/CyclotomicFieldNoUnitRoot` := proc(F,var,N,z)
    local t,zr,G,w,u,e,step,Gm,G2;
    if N=1 then
        G:=`a4b_solver/PrimitivePolynomial`(F,[var]): w:=var:
    else
        t:='_a4b_raw_diag_root': zr:=RootOf(cyclotomic(N,t)):
        G:=`a4b_solver/PrimitivePolynomial`(subs(z=zr,F),[var]): w:=var:
    end if:
    # Repeatedly compress all-even exponents before applying the two gcd
    # tests.  A root y of G(y^2) corresponds exactly to a root u=y^2.
    for step from 1 to 8 do
        u:=cat('_a4b_raw_even_',step): e:=`a4b_solver/EvenPartInVariable`(G,w,u):
        if not e["all_even"] then break end if:
        G:=e["polynomial"]: w:=u:
    end do:
    Gm:=gcd(G,expand(subs(w=-w,G))):
    G2:=gcd(G,expand(subs(w=w^2,G)*subs(w=-w^2,G))):
    return table(["no_unit_root"=evalb(degree(Gm,w)=0 and degree(G2,w)=0),
                  "gcd_sign"=Gm,"gcd_square"=G2,"compressed_polynomial"=G,
                  "compressed_variable"=w])
end proc:

# The previous field-only shortcut is not a complete torsion test when the
# coefficient field is Q(zeta_N): y -> y^2 need not preserve a minimal
# factor over that field.  Take the exact coefficient norm to Q first, then
# apply the gcd certificate over Q.  A root of unity of F is a root of its
# norm, so a Q-certificate of no root is conclusive for F.
`a4b_solver/CyclotomicNormNoUnitRoot` := proc(F,var,N,z)
    local RN,G,w,u,e,step,Gm,G2;
    RN := `a4b_solver/CyclotomicNorm`(F,N,z):
    G := `a4b_solver/PrimitivePolynomial`(RN,[var]):
    w := var:
    for step from 1 to 8 do
        u := cat('_a4b_norm_even_',step):
        e := `a4b_solver/EvenPartInVariable`(G,w,u):
        if not e["all_even"] then break end if:
        G := e["polynomial"]: w := u:
    end do:
    Gm := gcd(G,expand(subs(w=-w,G))):
    G2 := gcd(G,expand(subs(w=w^2,G)*subs(w=-w^2,G))):
    return table(["no_unit_root"=evalb(degree(Gm,w)=0 and degree(G2,w)=0),
                  "gcd_sign"=Gm,"gcd_square"=G2,
                  "normed_polynomial"=G,"compressed_variable"=w])
end proc:

# If both inputs are polynomials in elim^2, replace elim^2 by a new
# variable.  The original resultant is the square of the reduced one;
# the reduced polynomial has exactly the same zero set in the retained
# variables and is sufficient for torsion-candidate recovery.
`a4b_solver/CyclotomicResultantEvenAware` := proc(A,B,vars::list,elim,N,z)
    local u,v,keep,AA,BB,ea,eb,eka,ekb,ka,kb,R;
    u := '_a4b_even_elimination_variable':
    v := '_a4b_even_retained_variable':
    ea := `a4b_solver/EvenPartInVariable`(A,elim,u):
    eb := `a4b_solver/EvenPartInVariable`(B,elim,u):
    if ea["all_even"] and eb["all_even"] then
        AA:=ea["polynomial"]: BB:=eb["polynomial"]:
        keep := `a4b_solver/VariablesExcept`(vars,elim)[1]:
        ka := `a4b_solver/EvenPartInVariable`(AA,keep,v):
        kb := `a4b_solver/EvenPartInVariable`(BB,keep,v):
        if ka["all_even"] and kb["all_even"] then
            R:=`a4b_solver/CyclotomicResultant`(ka["polynomial"],kb["polynomial"],u,N,z):
            return expand(subs(v=keep^2,R))
        else
            return `a4b_solver/CyclotomicResultant`(AA,BB,u,N,z)
        end if:
    elif nops(vars)=2 then
        keep := `a4b_solver/VariablesExcept`(vars,elim)[1]:
        ka := `a4b_solver/EvenPartInVariable`(A,keep,v):
        kb := `a4b_solver/EvenPartInVariable`(B,keep,v):
        if ka["all_even"] and kb["all_even"] then
            R:=`a4b_solver/CyclotomicResultant`(ka["polynomial"],kb["polynomial"],elim,N,z):
            return expand(subs(v=keep^2,R))
        end if:
        return `a4b_solver/CyclotomicResultant`(A,B,elim,N,z)
    else
        return `a4b_solver/CyclotomicResultant`(A,B,elim,N,z)
    end if:
end proc:

`a4b_solver/CyclotomicFactorList` := proc(F,vars::list,N,z)
    local t,zr,Galg,ff,out,j,H;
    if N=1 then return factors(F)[2] end if:
    t := '_a4b_factor_cyclotomic_root':
    zr := RootOf(cyclotomic(N,t)):
    Galg := subs(z=zr,F):
    ff := factors(Galg)[2]: out := []:
    for j from 1 to nops(ff) do
        H := `a4b_solver/PrimitivePolynomial`(subs(zr=z,ff[j][1]),vars):
        out := [op(out),[H,ff[j][2]]]
    end do:
    return out:
end proc:

`a4b_solver/FactorUnivariateEvenAware` := proc(F,var,N,z)
    local u,e,ff,out,j,H,fac;
    u := '_a4b_even_factor_variable':
    e := `a4b_solver/EvenPartInVariable`(F,var,u):
    if e["all_even"] then
        fac := `a4b_solver/CyclotomicFactorList`(e["polynomial"],[u],N,z):
        out := []:
        for j from 1 to nops(fac) do
            H := expand(subs(u=var^2,fac[j][1])):
            out := [op(out),[H,fac[j][2]]]
        end do:
        return out
    else
        return `a4b_solver/CyclotomicFactorList`(F,[var],N,z)
    end if:
end proc:

`a4b_solver/ResultantStepCyclotomic` := proc(F,vars::list,elim,N,z)
    local G,others,h,R,ret,kept,v;
    G := `a4b_solver/PrimitivePolynomial`(F,vars):
    others := `a4b_solver/CyclotomicSignPowerTransforms`(G,vars,N,z):
    kept := `a4b_solver/VariablesExcept`(vars,elim):
    ret := []:
    for h in others do
        R := `a4b_solver/CyclotomicResultantEvenAware`(G,h["expression"],vars,elim,N,z):
        R := `a4b_solver/PrimitivePolynomial`(R,kept):
        ret := [op(ret),table(["id"=h["id"],"transform"=h["expression"],"resultant"=R])]:
    end do:
    return table(["input"=G,"variables"=vars,"eliminated"=elim,
                  "retained"=kept,"branches"=ret,
                  "branch_count"=nops(ret),"coefficient_field"=cat("Q(zeta_",N,")")]):
end proc:

`a4b_solver/ResultantStepGCDSplit` := proc(F,vars::list,elim)
    local G,others,h,R,H,Hraw,A,B,ret,kept,v,degH,unitH;
    G := `a4b_solver/PrimitivePolynomial`(F,vars):
    others := `a4b_solver/SignSquareTransforms`(G,vars):
    kept := `a4b_solver/VariablesExcept`(vars,elim):
    ret := []:
    for h in others do
        Hraw := gcd(G,h["expression"]):
        if Hraw=0 then Hraw:=G end if:
        H := `a4b_solver/PrimitivePolynomial`(Hraw,vars):
        degH := add(degree(H,v),v in vars):
        unitH := evalb(degH=0):
        if unitH then
            A := G: B := h["expression"]:
        else
            A := `a4b_solver/PrimitivePolynomial`(normal(G/Hraw),vars):
            B := `a4b_solver/PrimitivePolynomial`(normal(h["expression"]/Hraw),vars):
        end if:
        R := resultant(A,B,elim):
        R := `a4b_solver/PrimitivePolynomial`(R,kept):
        ret := [op(ret),table(["id"=h["id"],
            "transform"=h["expression"],"resultant"=R,
            "common_factor"=if unitH then 1 else H end if,
            "reduced_input"=A,"reduced_transform"=B,
            "has_common_factor"=not unitH])]:
    end do:
    return table(["input"=G,"variables"=vars,"eliminated"=elim,
        "retained"=kept,"branches"=ret,"branch_count"=nops(ret)]):
end proc:

`a4b_solver/ResultantStepGCDSplitCyclotomic` := proc(F,vars::list,elim,N,z)
    local G,others,h,R,H,Hraw,A,B,ret,kept,v,degG,degH,unitH,
          t,zr,Galg,Balg,Halg,Aalg,BalgRed,Ralg;
    G := `a4b_solver/PrimitivePolynomial`(F,vars):
    degG := add(degree(G,v),v in vars):
    others := `a4b_solver/CyclotomicMinimalTransforms`(G,vars,N,z):
    kept := `a4b_solver/VariablesExcept`(vars,elim):
    ret := []:
    if N<>1 then
        t := '_a4b_gcd_cyclotomic_root':
        zr := RootOf(cyclotomic(N,t)):
        Galg := subs(z=zr,G):
    end if:
    for h in others do
        if N=1 then
            Hraw := gcd(G,h["expression"]):
            if Hraw=0 then Hraw:=G end if:
            H := `a4b_solver/PrimitivePolynomial`(Hraw,vars):
        else
            Balg := subs(z=zr,h["expression"]):
            Halg := gcd(Galg,Balg):
            if Halg=0 then Halg:=Galg end if:
            H := `a4b_solver/PrimitivePolynomial`(subs(zr=z,Halg),vars):
        end if:
        degH := add(degree(H,v),v in vars):
        unitH := evalb(degH=0):
        if degH=degG then
            # The transformed equation contains the whole input factor.
            # Its residual quotient is (1,1), so this branch contributes no
            # residual resultant; the original factor remains represented
            # by the other branches and by the input tree itself.
            A := 1: B := 1: R := 1:
        elif unitH then
            A := G: B := h["expression"]:
            R := `a4b_solver/CyclotomicResultantEvenAware`(A,B,vars,elim,N,z):
        else
            if N=1 then
                A := `a4b_solver/PrimitivePolynomial`(normal(G/Hraw),vars):
                B := `a4b_solver/PrimitivePolynomial`(normal(h["expression"]/Hraw),vars):
                if degree(A,elim)=0 and degree(B,elim)=0 then R:=1 else R := `a4b_solver/CyclotomicResultantEvenAware`(A,B,vars,elim,N,z) end if:
            else
                Aalg := evala(normal(Galg/Halg)):
                BalgRed := evala(normal(Balg/Halg)):
                A := `a4b_solver/PrimitivePolynomial`(subs(zr=z,Aalg),vars):
                B := `a4b_solver/PrimitivePolynomial`(subs(zr=z,BalgRed),vars):
                if degree(Aalg,elim)=0 and degree(BalgRed,elim)=0 then
                    R := 1
                else
                    R := `a4b_solver/CyclotomicResultantEvenAware`(A,B,vars,elim,N,z):
                end if:
            end if:
        end if:
        R := `a4b_solver/PrimitivePolynomial`(R,kept):
        ret := [op(ret),table(["id"=h["id"],
            "transform"=h["expression"],"resultant"=R,
            "common_factor"=if unitH then 1 else H end if,
            "reduced_input"=A,"reduced_transform"=B,
            "has_common_factor"=not unitH])]:
    end do:
    return table(["input"=G,"variables"=vars,"eliminated"=elim,
        "retained"=kept,"branches"=ret,"branch_count"=nops(ret),
        "coefficient_field"=cat("Q(zeta_",N,")")]):
end proc:

`a4b_solver/ResultantStepNormed` := proc(F,vars::list,elim,N,z)
    local G,others,h,R,RN,ret,kept;
    G := `a4b_solver/PrimitivePolynomial`(F,vars):
    others := `a4b_solver/CyclotomicSignPowerTransforms`(G,vars,N,z):
    kept := `a4b_solver/VariablesExcept`(vars,elim):
    ret := []:
    for h in others do
        # Keep the coefficient field during elimination.  Take the norm
        # only after the dimension has dropped; this prevents degree blow-up.
        R := resultant(G,h["expression"],elim):
        RN := `a4b_solver/CyclotomicNorm`(R,N,z):
        RN := `a4b_solver/PrimitivePolynomial`(RN,kept):
        ret := [op(ret),table(["id"=h["id"],"transform"=h["expression"],"resultant"=RN])]:
    end do:
    return table(["input"=G,"variables"=vars,"eliminated"=elim,
                  "retained"=kept,"branches"=ret,
                  "branch_count"=nops(ret),"norm_order"=N]):
end proc:

`a4b_solver/RecursiveTorsion2D` := proc(F,vars::list)
    local step,b,diag,ret,elim,keep,fl,idxf;
    if degree(F,vars[1])<=degree(F,vars[2]) then
        elim := vars[1]: keep := [vars[2]]:
    else
        elim := vars[2]: keep := [vars[1]]:
    end if:
    step := `a4b_solver/ResultantStep`(F,vars,elim):
    keep := step["retained"]:
    ret := []:
    for b in step["branches"] do
        if b["resultant"]<>0 then
            fl := factors(b["resultant"])[2]:
            for idxf from 1 to nops(fl) do
                diag := `a4b_solver/UnitRootDiagnostics`(fl[idxf][1],keep[1]):
                ret := [op(ret),table(["id"=cat(b["id"],"/F",idxf),
                                       "source_id"=b["id"],
                                       "transform"=b["transform"],
                                       "resultant"=fl[idxf][1],"diagnostics"=diag])]:
            end do:
        end if:
    end do:
    return table(["dimension"=2,"step"=step,"univariate_branches"=ret]):
end proc:

`a4b_solver/RecursiveTorsion2DNormed` := proc(F,vars::list,N,z)
    local step,b,diag,ret,elim,keep,fl,idxf;
    if degree(F,vars[1])<=degree(F,vars[2]) then
        elim := vars[1]:
    else
        elim := vars[2]:
    end if:
    step := `a4b_solver/ResultantStepNormed`(F,vars,elim,N,z):
    keep := step["retained"]:
    ret := []:
    for b in step["branches"] do
        if b["resultant"]<>0 then
            fl := `a4b_solver/FactorUnivariateEvenAware`(b["resultant"],keep[1],N,z):
            for idxf from 1 to nops(fl) do
                diag := `a4b_solver/UnitRootDiagnostics`(fl[idxf][1],keep[1]):
                ret := [op(ret),table(["id"=cat(b["id"],"/F",idxf),
                                       "source_id"=b["id"],
                                       "transform"=b["transform"],
                                       "resultant"=fl[idxf][1],"diagnostics"=diag])]:
            end do:
        end if:
    end do:
    return table(["dimension"=2,"step"=step,"univariate_branches"=ret]):
end proc:

# Delayed norm variant.  Factor the resultant over the cyclotomic
# coefficient field first, then norm each univariate factor separately.
# This is algebraically equivalent for root containment, but avoids taking
# the norm of a large product and then factoring the inflated polynomial.
`a4b_solver/RecursiveTorsion2DDelayedNorm` := proc(F,vars::list,N,z)
    local step,b,diag,ret,elim,keep,fl,idxf,RN;
    if degree(F,vars[1])<=degree(F,vars[2]) then
        elim := vars[1]:
    else
        elim := vars[2]:
    end if:
    step := `a4b_solver/ResultantStepCyclotomic`(F,vars,elim,N,z):
    keep := step["retained"]:
    ret := []:
    for b in step["branches"] do
        if b["resultant"]<>0 then
            fl := `a4b_solver/FactorUnivariateEvenAware`(b["resultant"],keep[1],N,z):
            for idxf from 1 to nops(fl) do
                RN := `a4b_solver/CyclotomicNorm`(fl[idxf][1],N,z):
                RN := `a4b_solver/PrimitivePolynomial`(RN,keep):
                diag := `a4b_solver/UnitRootDiagnostics`(RN,keep[1]):
                ret := [op(ret),table(["id"=cat(b["id"],"/F",idxf),
                                       "source_id"=b["id"],
                                       "transform"=b["transform"],
                                       "resultant"=fl[idxf][1],
                                       "normed_resultant"=RN,
                                       "diagnostics"=diag])]:
            end do:
        end if:
    end do:
    return table(["dimension"=2,"step"=step,
                  "univariate_branches"=ret,
                  "norm_order"=N,"norm_strategy"="delayed_factorwise"]):
end proc:

# Adaptive 2D tree: use the formal-field gcd no-root test first.
# Only branches with a nontrivial torsion diagnostic are normed.
# Strict 2D tree: split each transformed pair into a common-factor
# subtree and a residual-resultant subtree.
`a4b_solver/RecursiveTorsion2DGCDSplitNorm` := proc(F,vars::list,N,z)
    local step,b,diag,ret,subtrees,elim,keep,fl,idxf,RN,H,degG,degH,v,
          process_resultant,rawdiag;
    if degree(F,vars[1])<=degree(F,vars[2]) then elim:=vars[1] else elim:=vars[2] end if:
    step := `a4b_solver/ResultantStepGCDSplitCyclotomic`(F,vars,elim,N,z):
    keep := step["retained"]: ret := []: subtrees := []:
    degG := add(degree(step["input"],v),v in vars):
    for b in step["branches"] do
        process_resultant:=true:
        if b["resultant"]<>0 and degree(b["resultant"],keep[1])>150 then
            rawdiag:=`a4b_solver/CyclotomicNormNoUnitRoot`(b["resultant"],keep[1],N,z):
            process_resultant:=not rawdiag["no_unit_root"]:
        end if:
        if b["resultant"]<>0 and process_resultant then
            fl := `a4b_solver/FactorUnivariateEvenAware`(b["resultant"],keep[1],N,z):
            for idxf from 1 to nops(fl) do
                RN := `a4b_solver/CyclotomicNorm`(fl[idxf][1],N,z):
                RN := `a4b_solver/PrimitivePolynomial`(RN,keep):
                diag := `a4b_solver/UnitRootDiagnostics`(RN,keep[1]):
                ret := [op(ret),table(["id"=cat(b["id"],"/R",idxf),
                    "source_id"=b["id"],"transform"=b["transform"],
                    "resultant"=fl[idxf][1],"normed_resultant"=RN,
                    "diagnostics"=diag])]:
            end do:
        end if:
        H := b["common_factor"]:
        degH := add(degree(H,v),v in vars):
        if degH>0 and degH<degG then
            subtrees := [op(subtrees),`a4b_solver/RecursiveTorsion2DDelayedNorm`(H,vars,N,z)]:
        end if:
    end do:
    return table(["dimension"=2,"step"=step,
        "univariate_branches"=ret,"common_subtrees"=subtrees,
        "norm_order"=N,"norm_strategy"="gcd_split_strict"]):
end proc:

`a4b_solver/RecursiveTorsion2DAdaptive` := proc(F,vars::list,N,z)
    local raw,rb,d0,ret,RN,fl,j,diag,keep;
    raw := `a4b_solver/RecursiveTorsion2D`(F,vars):
    keep := raw["step"]["retained"]:
    ret := []:
    for rb in raw["univariate_branches"] do
        d0 := rb["diagnostics"]:
        # Any nonzero constant gcd is a unit and is equivalent to 1 here.
        if (degree(d0["gcd_sign"],keep[1])=0 and
            degree(d0["gcd_square"],keep[1])=0) or
           nops(d0["cyclotomic_orders"])=0 then
            ret := [op(ret),rb]:
        else
            RN := `a4b_solver/CyclotomicNorm`(rb["resultant"],N,z):
            RN := `a4b_solver/PrimitivePolynomial`(RN,keep):
            fl := factors(RN)[2]:
            for j from 1 to nops(fl) do
                diag := `a4b_solver/UnitRootDiagnostics`(fl[j][1],keep[1]):
                ret := [op(ret),table(["id"=cat(rb["source_id"],"/NF",j),
                    "source_id"=rb["source_id"],
                    "transform"=rb["transform"],
                    "resultant"=rb["resultant"],
                    "normed_resultant"=fl[j][1],
                    "diagnostics"=diag])]:
            end do:
        end if:
    end do:
    return table(["dimension"=2,"step"=raw["step"],
        "univariate_branches"=ret,"norm_order"=N,
        "norm_strategy"="adaptive_gcd_then_factorwise"]):
end proc:

`a4b_solver/RecursiveTorsion3D` := proc(F,vars::list)
    local step,b,diag2,ret,vars2;
    step := `a4b_solver/ResultantStep`(F,vars,vars[3]):
    vars2 := [vars[1],vars[2]]: ret := []:
    for b in step["branches"] do
        if b["resultant"]<>0 then
            diag2 := `a4b_solver/RecursiveTorsion2D`(b["resultant"],vars2):
            ret := [op(ret),table(["id"=b["id"],"resultant"=b["resultant"],
                                   "two_variable_tree"=diag2])]:
        end if:
    end do:
    return table(["dimension"=3,"step"=step,"two_variable_branches"=ret,
                  "branch_count_3to2"=step["branch_count"],
                  "expected_2to1_per_branch"=2^3-1]):
end proc:

NULL:
