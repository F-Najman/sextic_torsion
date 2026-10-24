/* Alternative final [6] test for X1(49). Published J(Q) has a cyclic order-two 2-primary part. The exact
   JOne49 certificate proves T2 is nonzero on it, hence A2 kills all J(Q).
   The completed degree1..5 prefix remains a separate certificate. */
SetColumns(0);
SetMemoryLimit(8*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke.m");
Attach("fast_hecke_11.m");
Attach("fast_hecke_cyclic_commonzero.m");
Attach("hecke_two_commonzero.m");
t0:=Cputime();X:=MDX1(49,GF(3));
printf "C49_A2_SETUP genus=%o time=%o\n",Genus(X),Cputime(t0);
xs:=NoncuspidalPlaces(X,6);rs:=PlacesUpToDiamondFast(X,xs);
assert #xs eq 154 and #rs eq 8;
printf "C49_A2_PLACES total=%o representatives=%o time=%o\n",#xs,#rs,Cputime(t0);
survivors:=[];
for i in [1..#rs] do
 x:=rs[i];D:=Divisor(x);
 printf "C49_A2_START i=%o time=%o\n",i,Cputime(t0);
 T:=HeckeTwoCommonZeros(X,x);
 printf "C49_A2_HECKE i=%o time=%o\n",i,Cputime(t0);
 E:=EllipticCurve(X,x);P:=LevelStructure(X,x)`P;
 diamond:=ModuliPointCyclicCommonZeros(X,E,rec<recformat<P:PtEll>|P:=2*P>);
 printf "C49_A2_DIAMOND i=%o time=%o\n",i,Cputime(t0);
 A:=T-D-2*Divisor(diamond);
 principal:=IsPrincipal(A);
 printf "C49_A2_TEST i=%o principal=%o time=%o\n",i,principal,Cputime(t0);
 if principal then Append(~survivors,i);end if;
end for;
printf "C49_A2_DONE total=%o survivors=%o data=%o time=%o\n",#rs,#survivors,survivors,Cputime(t0);
quit;
