# Corrections confirmed during manual audit.
# APP-076 omitted the nonzero-coordinate factor (x-1) in the entered reference polynomial.
appendix_reference_corrections := table():
appendix_reference_corrections["APP-076"] := x-1:
appendix_reference_corrections["APP-230"] :=
  (zeta3^2-zeta3)*x*y^7+(zeta3^2-2*zeta3+1)*y^8-zeta3*x^2*y^5
  -(zeta3^2+zeta3-2)*x*y^6-2*zeta3^2*x^2*y^4
  -(zeta3^2-zeta3-2)*x*y^5-x^2*y^3+4*zeta3^2*x*y^4
  -zeta3*y^5-(zeta3^2-2*zeta3-1)*x*y^3-2*zeta3^2*y^4
  -(zeta3^2-2*zeta3+1)*x*y^2-y^3
  +(zeta3^2+zeta3-2)*x^2+(zeta3^2-1)*x*y:
NULL:
