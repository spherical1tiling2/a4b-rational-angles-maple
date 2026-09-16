restart:
read "recursive_torsion_engine.mpl":
read "root_backsubstitution.mpl":
read "exact_torus_backsub_nd.mpl":

U:='U': V:='V': x:='x': y:='y':
red:=table(["rank"=2,"status"="reduced","index"=2,
            "basis"=[[1,1],[1,-1]],"variables"=[U,V],
            "original_variables"=[x,y]]):
ru:=table(["variable"=U,"order"=1,"residue"=0]):
rv:=table(["variable"=V,"order"=1,"residue"=0]):
lifts:=`a4b_solver/LiftTorsionND`([ru,rv],red):
printf("LIFT_COUNT=%d ROOTS=%a\n",nops(lifts),
       [seq([[seq(lifts[j]["roots"][k]["order"],k=1..2)],
              [seq(lifts[j]["roots"][k]["residue"],k=1..2)]],
             j=1..nops(lifts))]):

fixed:=[table(["variable"=x,"order"=3,"residue"=1])]:
rec:=`a4b_solver/RecoverRemainingUnitRoots`(x+y,fixed,y):
printf("REC_STATUS=%s COUNT=%d DATA=%a\n",rec["status"],nops(rec["roots"]),
       [seq([rec["roots"][j]["order"],rec["roots"][j]["residue"]],
            j=1..nops(rec["roots"]))]):
quit:
