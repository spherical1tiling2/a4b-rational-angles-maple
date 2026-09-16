restart:
zeta4 := I: zeta8 := exp(I*Pi/4):
x := 'x': y := 'y':

# Reference polynomials copied from the six-case Appendix.
P1 := zeta4*x^4*y^8+x^4*y^7+2*zeta4*x^4*y^6-x^3*y^7-x^4*y^5
 +2*zeta4*x^3*y^6+zeta4*x^4*y^4+2*x^3*y^5-2*zeta4*x^3*y^4
 -2*x^2*y^5+zeta4^2*x^3*y^3+8*zeta4*x^2*y^4+x*y^5
 -2*zeta4^2*x^2*y^3-2*zeta4*x*y^4+2*zeta4^2*x*y^3+zeta4*y^4
 +2*zeta4*x*y^2-zeta4^2*y^3-zeta4^2*x*y+2*zeta4*y^2
 +zeta4^2*y+zeta4:

P2 := x^6*y^8-x^8*y^3-x^7*y^4+2*x^6*y^5-x^5*y^6+x^4*y^7
 -x^6*y^2-2*x^5*y^3-4*x^4*y^4-2*x^3*y^5-x^2*y^6+x^4*y
 -x^3*y^2+2*x^2*y^3-x*y^4-y^5+x^2:

P3 := x^10*y^7+x^10*y^6-x^9*y^6+x^9*y^5+x^8*y^6-2*x^8*y^5
 -x^8*y^4-2*x^7*y^5+2*x^7*y^4-3*x^6*y^4-x^6*y^3
 -2*x^5*y^4+2*x^5*y^3+x^4*y^4+3*x^4*y^3-2*x^3*y^3
 +2*x^3*y^2+x^2*y^3+2*x^2*y^2-x^2*y-x*y^2+x*y-y-1:

P4 := x^3*y^7+zeta4*x^2*y^6+x^3*y^5-zeta4*x^3*y^4+x^2*y^5
 +zeta4*x^2*y^4+x^3*y^3+x*y^5+2*zeta4*x*y^4+2*x^2*y^3
 +zeta4*x^2*y^2+zeta4*y^4+x*y^3+zeta4*x*y^2-y^3
 +zeta4*y^2+x*y+zeta4:

P5 := x^10*y^7+x^9*y^8-x^8*y^9+x^10*y^4-x^9*y^5-2*x^8*y^6
 +2*x^7*y^7-x^6*y^8+x^8*y^3-2*x^7*y^4-3*x^6*y^5
 +2*x^5*y^6-2*x^5*y^3+3*x^4*y^4+2*x^3*y^5-x^2*y^6
 +x^4*y-2*x^3*y^2+2*x^2*y^3+x*y^4-y^5+x^2-x*y-y^2:

P6 := zeta8^2*x^3*y^4-zeta8*x^4*y^2+zeta8^3*x^3*y^3
 +zeta8*x^2*y^4+x^5+zeta8^2*x^4*y-x^3*y^2
 +2*zeta8^2*x^2*y^3-x*y^4-zeta8^3*x^4+2*zeta8*x^3*y
 -zeta8^3*x^2*y^2+zeta8*x*y^3+zeta8^3*y^4+zeta8^2*x^3
 +x^2*y-zeta8^2*x*y^2+zeta8*x^2:

appendix_3v_reference_polynomials := [P1,P2,P3,P4,P5,P6]:
