restart:
with(LinearAlgebra):
read "recursive_torsion_engine.mpl":
read "generic_recursive_torsion_nd.mpl":

printf("=== FOUR FREE ANGLES: ONE VERTEX EQUATION PLUS THE FIVE-ANGLE TOTAL ===\n"):
aa:='_aa': bb:='_bb': cc:='_cc': dd:='_dd': ee:='_ee': q:='_q':
# The paper uses five angles and has already divided out Pi:
# alpha+beta+gamma+delta+epsilon = 3+4*q, q=1/f.
# We keep q as a variable, so six unknowns and two equations leave four free variables.
A4 := Matrix([[1,1,1,0,0,0],[1,1,1,1,1,-4]]):
b4 := Vector([2,3]):
R4 := ReducedRowEchelonForm(<A4|b4>):
printf("equations: alpha+beta+gamma=2; alpha+beta+gamma+delta+epsilon=3+4*q\n"):
printf("RREF=%a\n",R4):
printf("rank=%a, unknowns=%a, free_dimension=%a\n",Rank(A4),ColumnDimension(A4),ColumnDimension(A4)-Rank(A4)):
printf("solution: alpha=2-beta-gamma, epsilon=1+4*q-delta, q=1/f\n"):
printf("with f fixed, three angle variables remain free; with q included, four free variables remain.\n"):

# Use X=e^(i*Pi*alpha), ..., V=e^(i*Pi*epsilon).  If T=e^(i*Pi*4/f),
# the five-angle total becomes XYZWV+T=0 because exp(3*i*Pi)=-1.
X:='X': Y:='Y': Z:='Z': W:='W': V:='V': T:='T':
V4 := X*Y*Z-1:
S4 := X*Y*Z*W*V+T:
G4 := Groebner:-Basis([V4,S4],plex(X,Y,Z,W,V,T)):
printf("torus equations: %a=0, %a=0\n",V4,S4):
printf("Groebner projection of the two equations: %a\n",G4):
printf("after XYZ=1, the second equation is W*V+T=0.\n"):
printf("\n"):

printf("=== SIX ANGLES: ONLY THE TOTAL ANGLE EQUATION ===\n"):
a1:='_a1': a2:='_a2': a3:='_a3': a4:='_a4': a5:='_a5': a6:='_a6': q:='_q':
A := Matrix([[1,1,1,1,1,1,-4]]):
b := Vector([3]):
Aug := <A|b>:
R := ReducedRowEchelonForm(Aug):
printf("equation: a1+a2+a3+a4+a5+a6=3+4*q, q=1/f\n"):
printf("RREF=%a\n",R):
printf("rank=%a, unknowns=%a, free_dimension=%a\n",Rank(A),ColumnDimension(A),ColumnDimension(A)-Rank(A)):
printf("with f fixed, q is fixed and five angle variables remain free.\n"):
printf("torus form for fixed q: X1*X2*X3*X4*X5*X6 + exp(i*Pi*4/f)=0 up to the sign exp(3*i*Pi)=-1.\n"):
quit:
