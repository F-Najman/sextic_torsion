/* independent_counts.m -- independent recomputation of the numbers of non-cuspidal closed points of
   the curves of Table 4 (Section 11 of the paper), using Sutherland's GL2 package (https://github.com/AndrewVSutherland/Magma, file gl2points.m).

   That package counts #X_H(F_q) and #Y_H(F_q) from Frobenius similarity classes and weighted
   class numbers, without using a model of the modular curve, so the counts are independent of
   the models used in the sieves.

   X_1(m,n), m | n, is X_H for the pointwise stabiliser in GL(2,Z/n) of the level structure
   (P,Q) = ((n/m)e_1, e_2) of type Z/m x Z/n on E[n], namely
        H = { [a,0;c,1] : a = 1 mod m,  c = 0 mod m }.
   The genus of X_H is printed as a check on this identification.

   Usage: magma -n -b independent_counts.m   (uses gl2/magma.spec, see README.md) */
SetColumns(0);
SetMemoryLimit(2*10^9);
AttachSpec("gl2/magma.spec");

function LevelSubgroup(m,n)
    G := GL(2,Integers(n));
    gens := [ G![1,0,m,1] ]
        cat [ G![a,0,0,1] : a in [x : x in [1..n] | GCD(x,n) eq 1 and (x-1) mod m eq 0] ];
    return sub<G | gens>;
end function;

// <m, n, p, expected #Y_1(m,n)(F_{p^d}) for d = 1..5>, read off sextic model enumeration of the paper
pairs := [
<1,32,3,[0,0,24,32,80,864]>,
<1,34,3,[0,0,0,104,160,480]>,
<1,38,3,[0,0,9,0,360,531]>,
<1,33,5,[0,20,180,940,3010,16160]>,
<2,22,3,[0,0,0,0,300,600]>,
<1,39,2,[0,0,0,0,0,72]>,
<1,50,3,[0,0,0,0,0,360]>,
<1,40,3,[0,0,0,128,360,480]>,
<1,49,3,[0,0,0,84,315,924]>,
<2,26,3,[0,0,0,0,180,732]>
];

allok := true;
for t in pairs do
    m,n,p,exp := Explode(t);
    H := LevelSubgroup(m,n);
    assert GL2Level(H) eq n and #H eq (n div m)*#[x : x in [1..n] | GCD(x,n) eq 1 and (x-1) mod m eq 0];
    ycnt := [ &+GL2jCounts(H,p^d) : d in [1..6] ];
    ok := ycnt eq exp; allok and:= ok;
    printf "GL2COUNT X_1(%o,%o) p=%o: genus %o | #Y(F_{p^d}) = %o | sextic model enumeration = %o | %o\n",
        m, n, p, GL2Genus(H), ycnt, exp, ok select "agree" else "MISMATCH";
end for;
printf "GL2_POINT_COUNTS_DONE all_agree=%o\n", allok;
quit;
