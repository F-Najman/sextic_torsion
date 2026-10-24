SetColumns(0);
SetMemoryLimit(4*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke_11.m");
t0:=Cputime();X:=MDX11(4,16,GF(9));
assert #NoncuspidalPlaces(X,1) eq 0 and #NoncuspidalPlaces(X,2) eq 0;
xs:=NoncuspidalPlaces(X,3);assert #xs eq 320;
left:=Seqset(xs);rs:=[];
while #left gt 0 do
 x:=Rep(left);Append(~rs,x);
 for a in [1,3,5,7] do
  y:=DiamondOperator(X,a,x);assert y in xs;Exclude(~left,y);
 end for;
end while;
assert #rs eq 80;
print "SIEVE_4_16_P3_SETUP",Genus(X),#xs,#rs,Cputime(t0);
survivors:=[];
for i in [1..#rs] do
 x:=rs[i];D:=Divisor(x);killed:=false;
 for q in [5,13,17] do
  T:=HeckeOperatorFast(X,q,x);
  A:=T-D-q*Divisor(DiamondOperator(X,q,x));
  principal:=IsPrincipal(A);
  print "SIEVE_4_16_P3_TEST",i,q,principal,Cputime(t0);
  if not principal then killed:=true;break;end if;
 end for;
 if not killed then Append(~survivors,i);end if;
end for;
print "SIEVE_4_16_P3_DONE",#rs,#survivors,survivors,Cputime(t0);
quit;
