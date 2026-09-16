# Exact four-variable torus polynomial for the experimental system
#
#   alpha+beta+gamma = 2,
#   alpha+beta+gamma+delta+epsilon = 3+4/f,
#   formula (2.6) = 0.
#
# Output coordinates are
#   X = exp(I*Pi*beta),
#   Y = exp(I*Pi*gamma),
#   W = exp(2*I*Pi*delta),
#   T = -exp(4*I*Pi/f).
# The output support lattice is Z^4 (index one).

if not assigned(a4b_solver) then a4b_solver := table() end if:

`a4b_solver/Formula26FourVariablePolynomial` := proc(X,Y,W,T)
    local x0,y0,z0,t0,C0,S0,Flaurent,F,terms,term,ex,ey,ez,et,
          coeff0,G,eX,eY,eZ,eT,coeff1,H;

    x0 := '_a4b_half_beta':
    y0 := '_a4b_half_gamma':
    z0 := '_a4b_half_delta':
    t0 := '_a4b_half_fphase':
    C0 := P -> expand(2-P-1/P):
    S0 := P -> expand(P-1/P):

    # Common nonzero scalar -1/8 in formula (2.6) is omitted.
    Flaurent := expand(
        (C0(x0^2)*S0(-x0*y0*z0^2)
         -C0(y0^2)*S0(-x0*y0*t0^2/z0^2))
          *S0(z0^2/t0)
        -C0(x0^2/y0^2)*S0(-1/(x0*y0))*S0(t0)
    ):
    F := `a4b_solver/PrimitivePolynomial`(Flaurent,[x0,y0,z0,t0]):

    # First compression: X=x0^2, Y=y0^2, Z=z0^2, T=t0^2.
    G := 0:
    for term in [op(expand(F))] do
        ex := degree(term,x0): ey := degree(term,y0):
        ez := degree(term,z0): et := degree(term,t0):
        if irem(ex,2)<>0 or irem(ey,2)<>0 or
           irem(ez,2)<>0 or irem(et,2)<>0 then
            error "formula (2.6) did not have the expected even support"
        end if:
        coeff0 := normal(term/(x0^ex*y0^ey*z0^ez*t0^et)):
        # Keep z0^2 temporarily under the formal name _a4b_Z.
        G := G+coeff0*X^(ex/2)*Y^(ey/2)
                    *'_a4b_Z'^(ez/2)*T^(et/2):
    end do:
    G := expand(G):

    # The residual _a4b_Z support is even; W=_a4b_Z^2=z0^4.
    H := 0:
    for term in [op(G)] do
        eX := degree(term,X): eY := degree(term,Y):
        eZ := degree(term,'_a4b_Z'): eT := degree(term,T):
        if irem(eZ,2)<>0 then
            error "residual delta support was not even"
        end if:
        coeff1 := normal(term/(X^eX*Y^eY*'_a4b_Z'^eZ*T^eT)):
        H := H+coeff1*X^eX*Y^eY*W^(eZ/2)*T^eT:
    end do:
    return `a4b_solver/PrimitivePolynomial`(expand(H),[T,X,Y,W]):
end proc:

NULL:
