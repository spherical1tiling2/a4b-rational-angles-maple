restart:
read "recursive_torsion_engine.mpl":
read "four_to_three_unique.m":

T:='T': X:='X': Y:='Y':
P:=uniqueFactors[17]:
step:=`a4b_solver/ResultantStep`(P,[X,Y],X):
orders:=[]:

printf("=== TORSION PROJECTION U4_017 ===\n"):
printf("P=%a\n",P):
printf("branches=%d\n",step["branch_count"]):
for b in step["branches"] do
  if b["resultant"]<>0 then
    for entry in factors(b["resultant"])[2] do
      diag:=`a4b_solver/UnitRootDiagnostics`(entry[1],Y):
      printf("branch=%a factor=%a orders=%a gcdSign=%a gcdSquare=%a\n",
             b["id"],entry[1],diag["cyclotomic_orders"],
             diag["gcd_sign"],diag["gcd_square"]):
      for n in diag["cyclotomic_orders"] do
          if not member(n,orders) then orders:=[op(orders),n] end if
      end do:
    end do:
  end if:
end do:
printf("candidate X orders=%a\n",sort(orders)):
quit:
