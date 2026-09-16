restart:
read "recursive_torsion_engine.mpl":

# The paper's three-variable example after u=yz, w=y/z.
x:='x': u:='u': w:='w':
L := u^2*w*x^4+u^2*w^2*x^2-u^2*w*x^3-2*u*w*x^4-u^2*x^3
     +4*u*w*x^3-w^2*x^3+w*x^4-u^2*w*x+2*u^2*x^2-6*u*w*x^2
     +2*w^2*x^2-w*x^3+u^2*w-u^2*x+4*u*w*x-w^2*x-2*u*w-w*x+x^2+w:

tree3 := `a4b_solver/RecursiveTorsion3D`(L,[x,u,w]):
printf("3-variable to 2-variable branches = %a\n",tree3["branch_count_3to2"]):
printf("each 2-variable branch has nominally %a one-variable transforms\n",
       tree3["expected_2to1_per_branch"]):
printf("nominal total 3-to-1 branches = %a\n",
       tree3["branch_count_3to2"]*tree3["expected_2to1_per_branch"]):

# After a root-unit candidate x=zeta_N^r is identified, use this routine
# on the original L(x,u,w), leaving (u,w) as the two angle variables.
if assigned(_candidate_x) then
  Lfixed := expand(subs(x=_candidate_x,L)):
  fixed_tree := `a4b_solver/RecursiveTorsion2D`(Lfixed,[u,w]):
  printf("fixed x=%a; remaining two-variable branches=%a\n",
         _candidate_x,fixed_tree["step"]["branch_count"]):
end if:
quit:
