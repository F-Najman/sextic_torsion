/* Resume only the five mixed support types in the retained genus9
   quotient-filter session. All five single-place types completed in
   m2_30_quotient_prefix.log before a parser error between outer loops.
   No completed single-place test is fabricated or rerun here. */
assert [#places[d]:d in [1..6]] eq [0,24,192,732,2736,20000];
assert [#reps[d]:d in [1..6]] eq [0,6,48,186,684,5006];
assert Invariants(GY) eq [4,4,4,136,16320,0];
assert Invariants(HQY) eq [4,8160];
assert X`_E[1] eq 0 and X`_E[3] eq 0;
assert #survivors eq 20;
assert #[r:r in survivors|r[1] eq [5]] eq 2;
assert #[r:r in survivors|r[1] eq [6]] eq 18;
assert &and[r[1] in [[5],[6]]:r in survivors];
print "M2_30_MIXED_RESUME_SETUP retained_single_survivors=20 single_types_complete=5";
for ds in [[2,2],[3,2],[4,2],[3,3]] do
 da,db:=Explode(ds);passed:=0;total:=0;
 for x in reps[da] do
  for y in places[db] do
   a:=smallY[x]+smallY[y];total+:=1;
   if qY(a) in complementsY[6-da-db] then
    passed+:=1;Append(~survivors,[*ds,DescribePart(Divisor(x)+Divisor(y))*]);
    print "M2_30_QUOTIENT_SURVIVOR",ds,total,Eltseq(a),Cputime(t0);
   end if;
  end for;
 end for;
 assert total eq #reps[da]*#places[db];
 print "M2_30_QUOTIENT_TYPE",ds,total,passed,Cputime(t0);
end for;
passed:=0;total:=0;
for x in reps[2] do
 for i in [1..#places[2]] do
  for j in [i..#places[2]] do
   y:=places[2][i];z:=places[2][j];a:=smallY[x]+smallY[y]+smallY[z];total+:=1;
   if qY(a) eq QY!0 then
    passed+:=1;Append(~survivors,[*[2,2,2],DescribePart(Divisor(x)+Divisor(y)+Divisor(z))*]);
    print "M2_30_QUOTIENT_SURVIVOR",[2,2,2],total,Eltseq(a),Cputime(t0);
   end if;
  end for;
 end for;
end for;
assert total eq #reps[2]*Binomial(#places[2]+1,2);
print "M2_30_QUOTIENT_TYPE",[2,2,2],total,passed,Cputime(t0);
Write("data/m2_30_quotient_survivors.m","M2_30_QuotientSurvivors := " cat Sprint(survivors) cat ";");
print "M2_30_QUOTIENT_DATA_BEGIN";print survivors;print "M2_30_QUOTIENT_DATA_END";
printf "M2_30_QUOTIENT_RESUME_DONE survivors=%o class_maps=%o single_types=5 mixed_types=5 complete=true time=%o\n",#survivors,number,Cputime(t0);
quit;
