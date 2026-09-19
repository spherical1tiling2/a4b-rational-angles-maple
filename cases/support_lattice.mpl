with(LinearAlgebra):

`a4b_solver/SupportExponents2D` := proc(F,x,y)
    local G, termsList, term, E, ex, ey;

    G := expand(F):
    if G = 0 then
        return []:
    end if:

    if type(G,`+`) then
        termsList := [op(G)]:
    else
        termsList := [G]:
    end if:

    E := []:
    for term in termsList do
        ex := degree(term,x):
        ey := degree(term,y):
        if not member([ex,ey],E) then
            E := [op(E),[ex,ey]]:
        end if:
    end do:
    return E:
end proc:

`a4b_solver/AnalyzeSupportLattice2D` := proc(F,xyvars)
    local x, y, support, differences, d, i, j, rank,
          index, g, det, basis, T, H, hrow;

    if nops(xyvars) <> 2 then
        error "xyvars must contain exactly two variables":
    end if:
    x := xyvars[1]:
    y := xyvars[2]:
    support := `a4b_solver/SupportExponents2D`(F,x,y):

    T := table():
    T["polynomial"] := expand(F):
    T["support"] := support:

    if nops(support) = 0 then
        T["rank"] := 0:
        T["index"] := 0:
        T["is_full"] := false:
        T["difference_vectors"] := []:
        T["basis"] := []:
        return T
    end if:

    differences := []:
    for i from 2 to nops(support) do
        d := [support[i][1]-support[1][1],
              support[i][2]-support[1][2]]:
        if d <> [0,0] and not member(d,differences) then
            differences := [op(differences),d]:
        end if:
    end do:
    T["difference_vectors"] := differences:

    if nops(differences) = 0 then
        T["rank"] := 0:
        T["index"] := 0:
        T["is_full"] := false:
        T["basis"] := []:
        return T
    end if:

    rank := 1:
    for i from 1 to nops(differences)-1 do
        for j from i+1 to nops(differences) do
            det := differences[i][1]*differences[j][2]-
                   differences[i][2]*differences[j][1]:
            if det <> 0 then
                rank := 2:
                break:
            end if:
        end do:
        if rank = 2 then
            break:
        end if:
    end do:
    T["rank"] := rank:

    if rank = 1 then
        d := differences[1]:
        g := igcd(abs(d[1]),abs(d[2])):
        if g = 0 then g := 1 end if:
        basis := [[d[1]/g,d[2]/g]]:
        T["index"] := infinity:
        T["is_full"] := false:
        T["basis"] := basis:
        return T
    end if:

    # A first pair of independent support differences need not be a basis
    # of the full integer lattice (it may span a proper sublattice).  Use the
    # row Hermite normal form of all differences, which returns a genuine
    # Z-basis of the lattice they generate.  This is the step that prevents
    # false "non-integral lattice coordinate" failures.
    H := HermiteForm(Matrix(differences)):
    basis := []:
    for i from 1 to RowDimension(H) do
        hrow := [H[i,1],H[i,2]]:
        if hrow <> [0,0] then basis := [op(basis),hrow] end if:
    end do:
    if nops(basis)<>2 then error "Hermite form did not produce rank-two basis" end if:
    det := basis[1][1]*basis[2][2]-basis[1][2]*basis[2][1]:
    index := abs(det):
    T["index"] := index:
    T["is_full"] := evalb(index = 1):
    T["basis"] := basis:
    return T
end proc:

`a4b_solver/PrintLattice` := proc(T)
    printf("support = %a\n",T["support"]):
    printf("difference_vectors = %a\n",T["difference_vectors"]):
    printf("rank = %a, index = %a, is_full = %a\n",
        T["rank"],T["index"],T["is_full"]):
    printf("basis = %a\n",T["basis"]):
    return NULL:
end proc:

# Reduce a rank-two, non-full exponent lattice to a full-lattice polynomial.
# If b1,b2 is a lattice basis and e0 is one support exponent, write
#   e=e0+m1*b1+m2*b2,
# and replace the monomial X^e by U^m1 V^m2.  The removed monomial and the
# finite-index map are recorded for the exact torsion lift after solving the
# reduced curve.
`a4b_solver/FullLatticeReduction2D` := proc(F,xyvars::list)
    local x,y,L,b1,b2,e0,det,U,V,G,minm1,minm2,ex,dx,dy,m1,m2,c,
          terms,rank,index,rec;
    if nops(xyvars)<>2 then error "xyvars must contain exactly two variables" end if:
    x:=xyvars[1]: y:=xyvars[2]:
    L:=`a4b_solver/AnalyzeSupportLattice2D`(F,xyvars):
    rank:=L["rank"]: index:=L["index"]:
    if rank=1 then
        return table(["status"="rank1","polynomial"=F,
                      "variables"=xyvars,"original_variables"=xyvars,
                      "index"=infinity,"basis"=L["basis"],
                      "rank"=1,"lattice"=L])
    end if:
    if rank<>2 or index=1 then
        return table(["status"="full","polynomial"=F,
                      "variables"=xyvars,"original_variables"=xyvars,
                      "index"=1,"rank"=2,"basis"=[[1,0],[0,1]],
                      "base_exponent"=[0,0],"coordinate_shift"=[0,0]])
    end if:
    b1:=L["basis"][1]: b2:=L["basis"][2]: e0:=L["support"][1]:
    det:=b1[1]*b2[2]-b2[1]*b1[2]:
    U:='_a4b_lattice_U': V:='_a4b_lattice_V': G:=0:
    minm1:=0: minm2:=0:
    # AnalyzeSupportLattice2D returns a duplicate-free support list.  The
    # duplicate-free convention is important when a coefficient expands into
    # several zeta monomials.
    for ex in L["support"] do
        dx:=ex[1]-e0[1]: dy:=ex[2]-e0[2]:
        m1:=normal((b2[2]*dx-b2[1]*dy)/det):
        m2:=normal((-b1[2]*dx+b1[1]*dy)/det):
        if not type(m1,integer) or not type(m2,integer) then
            error "support exponent is not integral in the lattice basis"
        end if:
        c:=coeff(coeff(F,x,ex[1]),y,ex[2]):
        G:=G+c*U^m1*V^m2:
        if m1<minm1 then minm1:=m1 end if:
        if m2<minm2 then minm2:=m2 end if:
    end do:
    G:=expand(U^(-minm1)*V^(-minm2)*G):
    rec:=table(["status"="reduced","polynomial"=G,
                "variables"=[U,V],"original_variables"=xyvars,
                "index"=abs(det),"basis"=[b1,b2],
                "base_exponent"=e0,"coordinate_shift"=[minm1,minm2],
                "signed_determinant"=det,"lattice"=L]):
    return rec
end proc:

# Rank-one reduction.  If the support differences generate Z*d with d
# primitive, every monomial exponent is e0+m*d and, on the torus,
# F(x,y)=monomial(x,y)*h(x^d1*y^d2).  The monomial is nonzero, so the torsion
# problem is exactly the univariate root-of-unity problem h(t)=0.  The output
# keeps the Laurent shift and the direction explicitly; this represents the
# resulting one-parameter torsion family instead of silently discarding it.
`a4b_solver/Rank1Reduction2D` := proc(F,xyvars::list)
    local x,y,L,d,e0,t,h,minm,ex,dx,dy,m,c,terms;
    x:=xyvars[1]: y:=xyvars[2]:
    L:=`a4b_solver/AnalyzeSupportLattice2D`(F,xyvars):
    if L["rank"]<>1 then error "polynomial does not have rank-one support lattice" end if:
    d:=L["basis"][1]: e0:=L["support"][1]: t:='_a4b_rank1_t': h:=0: minm:=0:
    for ex in L["support"] do
        dx:=ex[1]-e0[1]: dy:=ex[2]-e0[2]:
        if d[1]<>0 then m:=normal(dx/d[1]) else m:=normal(dy/d[2]) end if:
        if not type(m,integer) then error "rank-one support coordinate is not integral" end if:
        c:=coeff(coeff(F,x,ex[1]),y,ex[2]): h:=h+c*t^m:
        if m<minm then minm:=m end if:
    end do:
    h:=expand(t^(-minm)*h):
    return table(["status"="rank1","polynomial"=h,"variable"=t,
                  "direction"=d,"base_exponent"=e0,
                  "coordinate_shift"=minm,"original_variables"=xyvars,
                  "lattice"=L])
end proc:
