SetColumns(0);
SetMemoryLimit(2*10^9);
F:=GF(64);
for n in [57,65] do
    qualifying:=0;
    for j in F do
        E0:=EllipticCurveFromjInvariant(j);
        for E in Twists(E0) do
            o:=#E;
            if o mod n ne 0 then continue; end if;
            qualifying+:=1;
            assert o eq n and j eq 0;
            A,mp:=AbelianGroup(E);
            assert Invariants(A) eq [n];
            P:=mp(A.1);
            auts:=Automorphisms(E);
            stabilizers:=[a:a in [1..n-1]|GCD(a,n) eq 1 and &or[u(P) eq a*P:u in auts]];
            assert 8 in stabilizers;
            printf "GOOD2 N=%o j=%o curve=%o order=%o aut_count=%o stabilizer_multipliers=%o\n",n,j,aInvariants(E),o,#auts,stabilizers;
        end for;
    end for;
    printf "GOOD2_ALL_CURVES N=%o qualifying_twists=%o\n",n,qualifying;
end for;
Z<t>:=PolynomialRing(Integers());
f:=t^3-t^2+t-2;
assert Discriminant(f) eq -83;
printf "C65_CUBIC_FIELD disc=%o factor_mod2=%o\n",Discriminant(f),Factorization(PolynomialRing(GF(2))!f);
print "GOOD2_57_65_CHECK_DONE";
quit;
