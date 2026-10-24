SetColumns(0);
SetMemoryLimit(1800000000);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke_11.m");
t0:=Cputime();
X:=MDX11(3,24,GF(7));
xs:=NoncuspidalPlaces(X,3);
reps:=PlacesUpToDiamondFast(X,xs);
assert #xs eq 48 and #reps eq 12;
print "SIEVE_3_24_SETUP",Genus(X),#xs,#reps,Cputime(t0);
survivors:=[];
for i in [1..#reps] do
 x:=reps[i];D:=Divisor(x);killed:=false;
 for q in [13,19] do
  ti:=Cputime();T:=HeckeOperatorFast(X,q,x);
  print "HECKE_STAGE",i,q,Cputime(ti);
  ti:=Cputime();A:=T-q*Divisor(DiamondOperator(X,q,x))-D;
  print "DIAMOND_STAGE",i,q,Cputime(ti);
  ti:=Cputime();principal:=IsPrincipal(A);
  print "PRINCIPAL_STAGE",i,q,principal,Cputime(ti);
  if not principal then killed:=true;break;end if;
 end for;
 if not killed then
  Append(~survivors,i);
  print "SIEVE_3_24_SURVIVOR",i;
 end if;
 print "SIEVE_3_24_PROGRESS",i,#reps,#survivors,Cputime(t0);
end for;
print "SIEVE_3_24_DONE",#reps,#survivors,Cputime(t0);
quit;
