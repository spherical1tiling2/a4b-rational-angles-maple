# Formal Laurent-polynomial form of equation (2.6).
# The phase basis is supplied as coefficient vectors in (a,b,c,d,e,q).

`a4b_solver/PhasePolynomial` := proc(T, xphase, yphase)
    local u1,u2,aa,bb,cc,dd,ee,qq,Lx,Ly,lx0,ly0,M,
          E,C,S, z, p, r, c0, L, a1,a2, coeffs, Zp;
    u1 := '_u1': u2 := '_u2':
    aa := T["affine_expressions"][1]:
    bb := T["affine_expressions"][2]:
    cc := T["affine_expressions"][3]:
    dd := T["affine_expressions"][4]:
    ee := T["affine_expressions"][5]:
    qq := T["affine_expressions"][6]:

    Lx := add(xphase[z]*T["affine_expressions"][z],z=1..6):
    Ly := add(yphase[z]*T["affine_expressions"][z],z=1..6):
    lx0 := subs(u1=0,u2=0,Lx):
    ly0 := subs(u1=0,u2=0,Ly):
    M := Matrix([[coeff(Lx,u1),coeff(Ly,u1)],
                 [coeff(Lx,u2),coeff(Ly,u2)]]):
    if LinearAlgebra:-Determinant(M)=0 then
        error "phase basis is singular on the affine model"
    end if:

    # Convert exp(i*pi*L) into x^p*y^r times the exact constant phase.
    E := proc(L)
        local l0,l1,l2,rr,constant;
        l0 := subs(u1=0,u2=0,L):
        l1 := coeff(L,u1): l2 := coeff(L,u2):
        rr := LinearAlgebra:-LinearSolve(M,Vector([l1,l2])):
        constant := simplify(l0-rr[1]*lx0-rr[2]*ly0,symbolic):
        return exp(I*Pi*constant)*x^rr[1]*y^rr[2]
    end proc:
    C := L -> (E(L)+E(-L))/2:
    S := L -> (E(L)-E(-L))/(2*I):

    Zp := ((1-C(bb))*S(dd-aa/2)-(1-C(cc))*S(ee-aa/2))*S((dd-ee)/2)
          -(1-C(bb-cc))*S(aa/2)*S((dd+ee)/2):
    return normal(expand(Zp)):
end proc:

# Fallback for a printed phase basis whose two coordinates become dependent
# after the vertex equations are imposed.  The polynomial is then generated
# in canonical free-parameter coordinates X=exp(i*pi*u1/D1),
# Y=exp(i*pi*u2/D2). The returned table records the coordinate denominators.
`a4b_solver/PhasePolynomialFreeCoordinates` := proc(T)
    local u1,u2,xf,yf,aa,bb,cc,dd,ee, Ls,D1,D2,L,E,C,S,Zp,j,L0,L1,L2;
    u1 := '_u1': u2 := '_u2': xf := 'Xfree': yf := 'Yfree':
    aa := T["affine_expressions"][1]: bb := T["affine_expressions"][2]:
    cc := T["affine_expressions"][3]: dd := T["affine_expressions"][4]:
    ee := T["affine_expressions"][5]:
    Ls := [bb,cc,dd-aa/2,ee-aa/2,(dd-ee)/2,bb-cc,
           aa/2,(dd+ee)/2]:
    D1 := 1: D2 := 1:
    for L in Ls do
        D1 := lcm(D1,denom(normal(coeff(L,u1)))):
        D2 := lcm(D2,denom(normal(coeff(L,u2)))):
    end do:
    E := proc(L)
        local c0,c1,c2;
        c0 := subs(u1=0,u2=0,L):
        c1 := coeff(L,u1): c2 := coeff(L,u2):
        return exp(I*Pi*c0)*xf^(c1*D1)*yf^(c2*D2)
    end proc:
    C := L -> (E(L)+E(-L))/2:
    S := L -> (E(L)-E(-L))/(2*I):
    Zp := ((1-C(bb))*S(dd-aa/2)-(1-C(cc))*S(ee-aa/2))*S((dd-ee)/2)
          -(1-C(bb-cc))*S(aa/2)*S((dd+ee)/2):
    return table(["polynomial"=normal(expand(Zp)),"D1"=D1,"D2"=D2]):
end proc:

NULL:
