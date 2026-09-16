# Canonical Laurent-polynomial generator.
# It never assumes that the printed phase variables are independent.
# X=exp(i*pi*u1/D1), Y=exp(i*pi*u2/D2), so all exponents are integral.

with(LinearAlgebra):

`a4b_solver/CanonicalPhasePolynomial` := proc(T, xphase, yphase)
    local u1,u2,X,Y,z,aa,bb,cc,dd,ee,Ls,D1,D2,N,L,c0,c1,c2,e,
          E,C,S,Zp,terms,t,dx,dy,minx,miny,poly,denoms,j;
    u1 := '_u1': u2 := '_u2': X := 'X': Y := 'Y': z := 'zeta':
    aa := T["affine_expressions"][1]:
    bb := T["affine_expressions"][2]:
    cc := T["affine_expressions"][3]:
    dd := T["affine_expressions"][4]:
    ee := T["affine_expressions"][5]:

    Ls := [bb,cc,dd-aa/2,ee-aa/2,(dd-ee)/2,bb-cc,
           aa/2,(dd+ee)/2]:
    D1 := 1: D2 := 1: N := 1:
    for L in Ls do
        c1 := normal(coeff(L,u1)): c2 := normal(coeff(L,u2)):
        D1 := lcm(D1,denom(c1)):
        D2 := lcm(D2,denom(c2)):
        c0 := normal(subs(u1=0,u2=0,L)):
        N := lcm(N,2*denom(c0)):
    end do:

    # Formal primitive N-th root zeta.  exp(i*pi*c0)=zeta^(N*c0/2).
    E := proc(L)
        local c0i,c1i,c2i,ei;
        c0i := normal(subs(u1=0,u2=0,L)):
        c1i := normal(coeff(L,u1)): c2i := normal(coeff(L,u2)):
        ei := normal(N*c0i/2):
        if not type(ei,integer) then
            error "phase constant was not integral in the chosen cyclotomic field"
        end if:
        return z^ei*X^(D1*c1i)*Y^(D2*c2i)
    end proc:
    C := L -> (E(L)+E(-L))/2:
    S := L -> (E(L)-E(-L))/(2*I):
    Zp := ((1-C(bb))*S(dd-aa/2)-(1-C(cc))*S(ee-aa/2))*S((dd-ee)/2)
          -(1-C(bb-cc))*S(aa/2)*S((dd+ee)/2):

    # Clear negative Laurent exponents; x,y are nonzero on the torus.
    terms := [op(expand(Zp))]:
    minx := 0: miny := 0:
    for t in terms do
        dx := degree(t,X): dy := degree(t,Y):
        if dx < minx then minx := dx end if:
        if dy < miny then miny := dy end if:
    end do:
    poly := expand(X^(-minx)*Y^(-miny)*Zp):
    return table(["polynomial"=poly,"raw_laurent"=Zp,
                  "N"=N,"D1"=D1,"D2"=D2,
                  "variables"=[X,Y],
                  "phase_coordinates"=[
                    add(xphase[j]*T["affine_expressions"][j],j=1..6),
                    add(yphase[j]*T["affine_expressions"][j],j=1..6)] ]):
end proc:

NULL:
