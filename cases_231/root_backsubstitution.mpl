# Root-unit candidate extraction and exact 2D back-substitution.

`a4b_solver/PrimitiveResidues` := proc(n)
    local r,out;
    if n=1 then return [0] end if:
    out := []:
    for r from 0 to n-1 do
        if igcd(r,n)=1 then out := [op(out),r] end if:
    end do:
    return out:
end proc:

`a4b_solver/RootValue` := proc(n,r)
    if n=1 then return 1 else return exp(2*Pi*I*r/n) end if:
end proc:

`a4b_solver/AlgebraicReduce` := proc(F)
    local G;
    G := simplify(F,symbolic):
    try G := evala(G) catch: NULL end try:
    return simplify(G,symbolic):
end proc:

`a4b_solver/SubstituteRoot` := proc(F,var,root,N,z)
    local zv,G;
    if N=1 then zv := 1 else zv := exp(2*Pi*I/N) end if:
    G := subs(z=zv,var=root,F):
    return `a4b_solver/AlgebraicReduce`(G):
end proc:

# Exact specialization in a common cyclotomic field.  If the coefficient
# field is Q(zeta_N) and the specialized root is zeta_n^r, work in
# Q(zeta_L), L=lcm(N,n), and express both roots as powers of one primitive
# L-th root.  This avoids floating-point gcds during branch recovery.
`a4b_solver/ExactSubstituteRoot` := proc(F,var,n,r,N,z)
    local L,t,zr,coeffRoot,rootValue,G;
    L := N*n/igcd(N,n):
    t := '_a4b_exact_common_root':
    zr := RootOf(cyclotomic(L,t)):
    if N=1 then coeffRoot:=1 else coeffRoot:=zr^(L/N) end if:
    if n=1 then rootValue:=1 else rootValue:=zr^(r*(L/n)) end if:
    G := subs(z=coeffRoot,var=rootValue,F):
    return evala(normal(G)):
end proc:

# Compute the exact common factor after fixing one root, and also return its
# coefficient norm over the common cyclotomic field.  The latter is needed
# because the remaining one-variable polynomial has coefficients in
# Q(zeta_lcm(N,n)); applying the rational one-variable cyclotomic test to the
# displayed RootOf coefficients directly is not valid.
`a4b_solver/ExactGCDAtRoot` := proc(F,H,keep,n,r,N,z,elim)
    local L,t,zr,w,coeffRoot,rootValue,Aalg,Balg,C,Cformal,Phi,
          terms,tm,minw,unit,RN;
    L:=N*n/igcd(N,n):
    t:='_a4b_exact_gcd_root': zr:=RootOf(cyclotomic(L,t)):
    w:='_a4b_formal_gcd_root':
    if N=1 then coeffRoot:=1 else coeffRoot:=zr^(L/N) end if:
    if n=1 then rootValue:=1 else rootValue:=zr^(r*(L/n)) end if:
    Aalg:=evala(normal(subs(z=coeffRoot,keep=rootValue,F))):
    Balg:=evala(normal(subs(z=coeffRoot,keep=rootValue,H))):
    C:=gcd(Aalg,Balg):
    if evalb(C=0) then return table(["algebraic"=C,"formal"=0,"normed"=0,"order"=L]) end if:
    Cformal:=subs(zr=w,C):
    Phi:=cyclotomic(L,w):
    terms:=[op(expand(Cformal))]:
    if nops(terms)=0 then return table(["algebraic"=C,"formal"=Cformal,"normed"=0,"order"=L]) end if:
    minw:=min(seq(degree(tm,w),tm in terms)):
    Cformal:=expand(w^(-minw)*Cformal):
    RN:=expand(resultant(Phi,Cformal,w)):
    unit:=(-1)^(degree(Phi,w)*minw):
    RN:=`a4b_solver/PrimitivePolynomial`(expand(unit*RN),[elim]):
    return table(["algebraic"=C,"formal"=Cformal,"normed"=RN,"order"=L])
end proc:

`a4b_solver/ExactSubstituteTwoRoots` := proc(F,var1,n1,r1,var2,n2,r2,N,z)
    local L,t,zr,coeffRoot,root1,root2,G;
    L := N*n1*n2/igcd(N,n1)/igcd(N*n1/igcd(N,n1),n2):
    t := '_a4b_exact_two_common_root':
    zr := RootOf(cyclotomic(L,t)):
    if N=1 then coeffRoot:=1 else coeffRoot:=zr^(L/N) end if:
    if n1=1 then root1:=1 else root1:=zr^(r1*(L/n1)) end if:
    if n2=1 then root2:=1 else root2:=zr^(r2*(L/n2)) end if:
    G := subs(z=coeffRoot,var1=root1,var2=root2,F):
    return evala(normal(G)):
end proc:

# Convert an exponent p modulo M, representing zeta_M^p, to the canonical
# pair (order,residue) with zeta_order^residue = zeta_M^p.  Keeping this
# conversion exact is important: the finite-index lattice lift must not use
# numerical angles or floating-point root comparisons.
`a4b_solver/ExponentToRootData` := proc(p,M)
    local pp,g,n,r;
    pp:=irem(p,M):
    if pp=0 then return table(["order"=1,"residue"=0]) end if:
    g:=igcd(pp,M): n:=M/g: r:=irem(pp/g,n):
    return table(["order"=n,"residue"=r]):
end proc:

# Lift every torsion point of a full-lattice reduction back to the original
# variables.  If U=x^b1x*y^b1y and V=x^b2x*y^b2y, the reduced candidate gives
# two congruences for p,q in x=eta_M^p, y=eta_M^q.  Taking
# M=index*lcm(ord(U),ord(V)) makes the Smith-normal-form kernel visible:
# exactly |det(b1,b2)| solutions occur.  The explicit adjugate formula gives
# one solution modulo L=lcm(ord(U),ord(V)); the complete kernel is obtained
# by enumerating the L-periodic lifts and filtering the two congruences.
`a4b_solver/LiftReducedTorsionCandidates` := proc(candidates,red)
    local out,vars,orig,basis,b1,b2,det,D,idx,cand,ru,rv,nu,nv,L,M,
          a,b,t1,t2,A,B,p0,q0,k,l,p,q,datx,daty,key,seen,tab,
          reduced_retained,reduced_eliminated;
    out:=[]: seen:=table():
    vars:=red["variables"]: orig:=red["original_variables"]:
    if red["status"]="full" then
        for cand in candidates do
            if cand["retained"]=orig[1] then
                tab:=table(["x_order"=cand["retained_order"],
                           "x_residue"=cand["retained_residue"],
                           "y_order"=cand["eliminated_order"],
                           "y_residue"=cand["eliminated_residue"],
                           "branch"=cand["branch"],
                           "reduced_candidate"=cand]):
            else
                tab:=table(["x_order"=cand["eliminated_order"],
                           "x_residue"=cand["eliminated_residue"],
                           "y_order"=cand["retained_order"],
                           "y_residue"=cand["retained_residue"],
                           "branch"=cand["branch"],
                           "reduced_candidate"=cand]):
            end if:
            key:=cat(tab["x_order"],"/",tab["x_residue"],"|",
                     tab["y_order"],"/",tab["y_residue"]):
            if not assigned(seen[key]) then
                seen[key]:=true: out:=[op(out),copy(tab)]
            end if:
        end do:
        return out:
    end if:
    if red["status"]<>"reduced" then return out end if:
    basis:=red["basis"]: b1:=basis[1]: b2:=basis[2]:
    det:=red["signed_determinant"]: D:=abs(det): idx:=red["index"]:
    if D<>idx or D=0 then error "invalid finite-index lattice data" end if:
    for cand in candidates do
        if cand["retained"]=vars[1] then
            nu:=cand["retained_order"]: ru:=cand["retained_residue"]:
            nv:=cand["eliminated_order"]: rv:=cand["eliminated_residue"]:
        else
            nu:=cand["eliminated_order"]: ru:=cand["eliminated_residue"]:
            nv:=cand["retained_order"]: rv:=cand["retained_residue"]:
        end if:
        L:=ilcm(nu,nv): M:=D*L:
        # U=eta_M^(a*D), V=eta_M^(b*D), where eta_M is primitive.
        a:=if nu=1 then 0 else ru*(L/nu) end if:
        b:=if nv=1 then 0 else rv*(L/nv) end if:
        t1:=a*D: t2:=b*D:
        A:=b2[2]*t1-b1[2]*t2:
        B:=-b2[1]*t1+b1[1]*t2:
        if not type(A/det,integer) or not type(B/det,integer) then
            error "lattice lift adjugate is not integral"
        end if:
        p0:=irem(A/det,L): q0:=irem(B/det,L):
        for k from 0 to D-1 do
            for l from 0 to D-1 do
                p:=p0+k*L: q:=q0+l*L:
                if irem(b1[1]*p+b1[2]*q-t1,M)=0 and
                   irem(b2[1]*p+b2[2]*q-t2,M)=0 then
                    datx:=`a4b_solver/ExponentToRootData`(p,M):
                    daty:=`a4b_solver/ExponentToRootData`(q,M):
                    key:=cat(datx["order"],"/",datx["residue"],"|",
                             daty["order"],"/",daty["residue"]):
                    if not assigned(seen[key]) then
                        seen[key]:=true:
                        tab:=table(["x_order"=datx["order"],
                                   "x_residue"=datx["residue"],
                                   "y_order"=daty["order"],
                                   "y_residue"=daty["residue"],
                                   "branch"=cat(cand["branch"],";lift"),
                                   "reduced_candidate"=cand,
                                   "lift_modulus"=M,
                                   "lift_index"=D]):
                        out:=[op(out),copy(tab)]
                    end if:
                end if:
            end do:
        end do:
    end do:
    return out:
end proc:

# Compose two finite-index lifts.  This is used when a coefficient-field
# factor is first reduced from (X,Y) to (U,V), and only then normed over Q.
# The normed polynomial is solved in a further coordinate system if needed;
# the first lift returns (U,V)-torsion data and the second returns (X,Y).
`a4b_solver/LiftReducedTorsionCandidatesComposite` := proc(candidates,redNorm,redBase)
    local mid,mc,out,c,tab,uv1,uv2,seen,key;
    mid:=`a4b_solver/LiftReducedTorsionCandidates`(candidates,redNorm):
    out:=[]: seen:=table():
    for c in mid do
        mc:=table(["retained"=redBase["variables"][1],
                   "retained_order"=c["x_order"],
                   "retained_residue"=c["x_residue"],
                   "eliminated"=redBase["variables"][2],
                   "eliminated_order"=c["y_order"],
                   "eliminated_residue"=c["y_residue"],
                   "branch"=c["branch"]]):
        for tab in `a4b_solver/LiftReducedTorsionCandidates`([mc],redBase) do
            key:=cat(tab["x_order"],"/",tab["x_residue"],"|",
                     tab["y_order"],"/",tab["y_residue"]):
            if not assigned(seen[key]) then
                seen[key]:=true: out:=[op(out),copy(tab)]
            end if:
        end do:
    end do:
    return out:
end proc:

`a4b_solver/Recover2DTorsionCandidates` := proc(tree,N,z)
    local F,vars,elim,keep,b,orders,n,r,rootR,A,B,C,diag2,orders2,
          n2,r2,rootE,check,out,source_id,diag0,gcdinfo,nontrivial;
    F := tree["step"]["input"]:
    vars := tree["step"]["variables"]:
    elim := tree["step"]["eliminated"]:
    keep := tree["step"]["retained"]:
    out := []:
    for b in tree["univariate_branches"] do
      try
        source_id := b["source_id"]:
        # The stored diagnostic can be stale when the branch table was built
        # while factoring several normed factors.  Recompute it from the
        # exact rational norm whenever that field is present.
        try
            diag0 := `a4b_solver/UnitRootDiagnostics`(
                         b["normed_resultant"],keep[1]):
            orders := diag0["cyclotomic_orders"]:
        catch:
            # Some branches have a zero/constant normed resultant, for which
            # the one-variable diagnostic is intentionally undefined.  The
            # historical tree field is only a label in those branches, so it
            # must not be indexed as a table; such a branch contributes no
            # isolated retained root.
            orders := []:
        end try:
        for n in orders do
            for r in `a4b_solver/PrimitiveResidues`(n) do
              try
                gcdinfo := `a4b_solver/ExactGCDAtRoot`(
                              F,b["transform"],keep[1],n,r,N,z,elim):
                C := gcdinfo["algebraic"]:
                # A zero gcd means the two specialized polynomials are
                # identically zero; it is not a finite isolated root
                # candidate and has no univariate factor to diagnose.
                nontrivial:=true:
                if evalb(C=0) then nontrivial:=false end if:
                if nontrivial then if evalb(C=1) then nontrivial:=false end if end if:
                if nontrivial then if evalb(C=-1) then nontrivial:=false end if end if:
                if nontrivial then
                    try
                        diag2 := `a4b_solver/UnitRootDiagnostics`(
                                    gcdinfo["normed"],elim):
                        orders2 := diag2["cyclotomic_orders"]:
                    catch:
                        orders2 := []:
                    end try:
                    for n2 in orders2 do
                        for r2 in `a4b_solver/PrimitiveResidues`(n2) do
                            check := `a4b_solver/ExactSubstituteTwoRoots`(
                                      F,keep[1],n,r,elim,n2,r2,N,z):
                            if evalb(check=0) then
                                rootR:=`a4b_solver/RootValue`(n,r):
                                rootE:=`a4b_solver/RootValue`(n2,r2):
                                out := [op(out),table(["branch"=b["id"],
                                      "source_id"=source_id,"retained"=keep[1],
                                      "retained_order"=n,"retained_residue"=r,
                                      "eliminated"=elim,"eliminated_order"=n2,
                                      "eliminated_residue"=r2,
                                      "retained_root"=rootR,"eliminated_root"=rootE])]:
                            end if:
                        end do:
                    end do:
                end if:
              catch:
                printf("RECOVER_POINT_ERROR %s n=%a r=%a detail=%a\\n",b["id"],n,r,lasterror):
              end try:
            end do:
        end do:
      catch:
        printf("RECOVER_BRANCH_ERROR %s detail=%a\\n",b["id"],lasterror):
      end try:
    end do:
    return out:
end proc:

NULL:
