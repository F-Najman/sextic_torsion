SetColumns(0);SetMemoryLimit(2*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke_11.m");
Attach("fast_diamond_hash_m2.m");
t0:=Cputime();X:=MDX11(2,24,GF(5));
for dd in [3..6] do
 xs:=NoncuspidalPlaces(X,dd);rs:=PlacesUpToDiamondM2HashDirect(X,xs);
 print "M2_DIRECT_HASH_COUNT",dd,#xs,#rs,Cputime(t0);
 if dd in [3,4] then
  K:=BaseRing(EllipticCurve(X,xs[1]));
  // Compare full orbit partitions, including degree4 short orbits.
  left:=Seqset(xs);old:={};
  while #left gt 0 do
   x:=Rep(left);orbit:={DiamondOperator(X,a,x):a in [1,5,7,11]};
   for y in orbit do assert y in xs;Exclude(~left,y);end for;
   Include(~old,orbit);
  end while;
  ids:=AssociativeArray();
  for x in xs do
   id:=Sort(MDM2DiamondKeysDirect(X,x,K))[1];
   if IsDefined(ids,id) then Include(~ids[id],x);else ids[id]:={x};end if;
  end for;
  new:={ids[id]:id in Keys(ids)};
  assert new eq old and #new eq #rs;
  print "M2_DIRECT_HASH_FULL_PARTITION_MATCH",dd,#old,Cputime(t0);
 end if;
 if dd eq 3 then assert #rs eq 4;end if;
 if dd eq 4 then assert #rs eq 32;end if;
 if dd eq 5 then assert #rs eq 168;end if;
end for;
print "M2_DIRECT_HASH_TEST_DONE complete=true",Cputime(t0);
quit;
