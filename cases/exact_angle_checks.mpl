# Exact checks for angles measured in units of Pi.
`a4b_solver/ExactEquation26` := proc(angles)
    local aa,bb,cc,dd,ee,L,t,S,C,r,P;
    aa,bb,cc,dd,ee:=op(angles):
    L:=ilcm(seq(denom(r),r in angles)):
    S:=r->(t^(2*L*r)-t^(-2*L*r))/2:
    C:=r->(t^(2*L*r)+t^(-2*L*r))/2:
    # Each term has two sine factors; their common factor I^(-2) cancels.
    r:=((1-C(bb))*S(dd-aa/2)-(1-C(cc))*S(ee-aa/2))*S((dd-ee)/2)
       -(1-C(bb-cc))*S(aa/2)*S((dd+ee)/2):
    P:=numtheory:-cyclotomic(4*L,t):
    return evalb(rem(numer(normal(r)),P,t)=0)
end proc:

`a4b_solver/FamilyEquation26` := proc(angles,q)
    local aa,bb,cc,dd,ee,t,E,S,C,r;
    aa,bb,cc,dd,ee:=op(angles):
    E:=r->simplify(evalc(exp(I*Pi*subs(q=0,r))))*t^coeff(r,q):
    S:=r->(E(r)-E(-r))/2:
    C:=r->(E(r)+E(-r))/2:
    r:=((1-C(bb))*S(dd-aa/2)-(1-C(cc))*S(ee-aa/2))*S((dd-ee)/2)
       -(1-C(bb-cc))*S(aa/2)*S((dd+ee)/2):
    return evalb(normal(r)=0)
end proc:
