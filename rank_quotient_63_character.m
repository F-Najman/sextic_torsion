/* Exact replay of the character kernel recorded in the quintic proof. */
SetColumns(0);
AttachSpec("mdmagma/Magma/magma.spec");
units:=[a:a in [1..62]|GCD(a,63) eq 1];
ker:={1,5,25,38,58,62};matching:=[];
for c in units do
 if {a:a in units|ConreyCharacterValue(63,c,a) eq 1} eq ker then
  assert ConreyCharacterValue(63,c,8) eq -1;Append(~matching,c);
 end if;
end for;
assert #matching eq 2;
print "C63_POSITIVE_CHARACTER_QUOTIENT_DONE",matching,"diamond8=-1";
quit;
