restart:
read "recursive_torsion_engine.mpl":
read "formula26_four_variable_core.mpl":

T:='T': X:='X': Y:='Y': W:='W':
F:=`a4b_solver/Formula26FourVariablePolynomial`(X,Y,W,T):

printf("=== BACK-SUBSTITUTION OF RANK-DEFICIENT 4->3 PROJECTIONS ===\n"):
printf("Y=1 (gamma=0): %a\n",factor(subs(Y=1,F))):
printf("X=1 (beta=0): %a\n",factor(subs(X=1,F))):
printf("T=1 (forces excluded f=4): %a\n",factor(subs(T=1,F))):
printf("T=-1 (no finite f): %a\n",factor(subs(T=-1,F))):
printf("X=Y (beta=gamma): %a\n",factor(subs(Y=X,F))):
printf("X*Y=1 (boundary alpha=0): %a\n",factor(expand(subs(Y=1/X,F)*X^3))):

# T*X*Y+1=0 is equivalent, in the admissible angle range, to alpha=4/f.
Falpha4 := `a4b_solver/PrimitivePolynomial`(subs(T=-1/(X*Y),F),[X,Y,W]):
printf("T*X*Y+1=0 (alpha=4/f): %a\n",factor(Falpha4)):

quit:
