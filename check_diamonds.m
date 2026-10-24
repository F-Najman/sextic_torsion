/* Lightweight exact character checks; no modular symbols computation. */
SetColumns(0);
AttachSpec("mdmagma/Magma/magma.spec");
for t in [ <91, [9,30], [3]>, <95, [8], [11]>,
           <133, [39,31], [64]>, <247, [8], [12,56]> ] do
    n, chars, gens := Explode(t);
    for c in chars do
        for a in gens do
            v := ConreyCharacterValue(n,c,a);
            printf "N=%o chi=%o a=%o value=%o\n",n,c,a,v;
            assert v eq 1;
        end for;
    end for;
    z := Integers(n);
    norm := func<a | Minimum(Integers()!(z!a), Integers()!(-z!a))>;
    H := {1};
    repeat
        old := H;
        H join:= { norm(a*b) : a in H, b in gens };
    until H eq old;
    printf "N=%o generators=%o subgroup_mod_sign=%o order=%o\n", n,gens,Sort(Setseq(H)),#H;
end for;
n := 133;
units := [a : a in [1..n-1] | GCD(a,n) eq 1];
ordchar := func<v | Minimum([j : j in [1..12] | v^j eq 1])>;
gs := [a : a in units | ConreyCharacterValue(n,39,a) eq 1 and ordchar(ConreyCharacterValue(n,31,a)) eq 6];
hs := [a : a in units | ConreyCharacterValue(n,31,a) eq 1 and ordchar(ConreyCharacterValue(n,39,a)) eq 3];
printf "N=133 g choices=%o\n", gs;
printf "N=133 h choices=%o\n", hs;
for a in [gs[1],hs[1]] do
    printf "N=133 a=%o chi39=%o chi31=%o\n", a,ConreyCharacterValue(n,39,a),ConreyCharacterValue(n,31,a);
end for;
ts := [a : a in [1..94] | GCD(a,95) eq 1 and ConreyCharacterValue(95,8,a) eq -1];
printf "N=95 t choices with chi8=-1: %o\n", ts;
print "DIAMOND_CHECKS_DONE";
quit;
