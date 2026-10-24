/* Noncyclic extension of fast_hecke.m.
   The kernel-orbit algorithm is identical; it carries both level points.
   Attach this after the v2 mdmagma spec. */

intrinsic MDIsogeniesFastPoints(E::CrvEll[FldFin], points::SeqEnum[PtEll], q::RngIntElt) -> List
{ For each Frobenius-orbit of cyclic subgroups G of order q of E (over the base field K of E), return
  <E/G, images of points, orbit size e>, with E/G and the images defined over the field of definition L of G }
    require IsPrime(q) and IsOdd(q) and q ne Characteristic(BaseRing(E)) : "q must be an odd prime different from the characteristic";
    K := BaseRing(E);
    fq := DivisionPolynomial(E, q);
    fac := [t[1] : t in Factorization(fq)];
    assert &+[Degree(h) : h in fac] eq (q^2-1) div 2;   // separable: distinct x-coordinates
    a1, a2, a3, a4, a6 := Explode(aInvariants(E));
    unassigned := [true : h in fac];
    result := [* *];
    total := 0;
    for i in [1..#fac] do
        if not unassigned[i] then continue; end if;
        h := fac[i]; d := Degree(h);
        if d eq 1 then F1 := K; x1 := -Coefficient(h, 0)/Coefficient(h, 1);
        else F1 := ext<K | h>; x1 := F1.1; end if;
        // y-coordinate: solve y^2 + a1 x y + a3 y = x^3 + a2 x^2 + a4 x + a6 over F1 or a quadratic extension
        PY<Y> := PolynomialRing(F1);
        ypol := Y^2 + (F1!a1*x1 + F1!a3)*Y - (x1^3 + F1!a2*x1^2 + F1!a4*x1 + F1!a6);
        rts := Roots(ypol);
        if #rts gt 0 then F2 := F1; y1 := rts[1][1];
        else F2 := ext<F1 | 2>; y1 := Roots(PolynomialRing(F2)!ypol)[1][1]; end if;
        E2 := EllipticCurve([F2!a1, F2!a2, F2!a3, F2!a4, F2!a6]);
        P1 := E2![F2!x1, y1];
        assert q*P1 eq E2!0 and P1 ne E2!0;
        xs := [ F1!((k*P1)[1]) : k in [1..(q-1) div 2] ];          // x-coordinates lie in F1
        PX<X> := PolynomialRing(F1);
        g := &*[X - xk : xk in xs];                                  // kernel polynomial of G = <P1>
        cs := Coefficients(g);
        L := sub< F1 | cs cat [F1!K.1] >;                            // field of definition of G, contains K
        e := Degree(L) div Degree(K);
        // mark the factors belonging to this Frobenius-orbit: divisors of Norm_{L/K}(g)
        qK := #K;
        Ng := PX!1; gg := g;
        for j in [1..e] do
            Ng *:= gg;
            gg := PX![c^qK : c in Coefficients(gg)];
        end for;
        assert gg eq g;   // Frobenius^e fixes g
        NgK := PolynomialRing(K)![K!c : c in Coefficients(Ng)];
        deg_assigned := 0;
        for j in [1..#fac] do
            if unassigned[j] and IsDivisibleBy(NgK, fac[j]) then unassigned[j] := false; deg_assigned +:= Degree(fac[j]); end if;
        end for;
        assert deg_assigned eq e*(q-1) div 2;
        // the isogeny over L
        PL<XL> := PolynomialRing(L);
        gL := PL![L!c : c in cs];
        EL := EllipticCurve([L!(F1!a) : a in aInvariants(E)]);
        EG, phi := IsogenyFromKernel(EL, gL);
        images := [phi(EL![L!(F1!(K!P[1])), L!(F1!(K!P[2]))]) : P in points];
        Append(~result, <EG, images, e>);
        total +:= e;
    end for;
    assert total eq q+1;
    return result;
end intrinsic;

intrinsic HeckeOperatorFast(X::MDCrvMod11, q::RngIntElt, x::PlcCrvElt) -> DivCrvElt
{ T_q(x) on X_1(2,N), carrying both level points through each kernel orbit }
    require X`M eq 2 : "ModuliPoint is currently implemented only for M=2";
    require IsPrime(q) and IsOdd(q) : "q must be an odd prime";
    require GCD(Characteristic(BaseRing(X)),q) eq 1 and GCD(Level(X),q) eq 1 : "q must be coprime to the characteristic and level";
    E := EllipticCurve(X,x);
    Ls := LevelStructure(X,x);
    isog := MDIsogeniesFastPoints(E,[Ls`P,Ls`Q],q);
    D := DivisorGroup(Curve(X))!0;
    for t in isog do
        EG, points, e := Explode(t);
        y := ModuliPoint(X,EG,rec<recformat<P:PtEll,Q:PtEll> | P:=points[1],Q:=points[2]>);
        mult := e*Degree(x)/Degree(y);
        assert Denominator(mult) eq 1;
        D +:= (Integers()!mult)*Divisor(y);
    end for;
    assert Degree(D) eq (q+1)*Degree(x);
    return D;
end intrinsic;

intrinsic HeckeOperatorFast(X::MDCrvMod11, q::RngIntElt, D::DivCrvElt) -> DivCrvElt
{ T_q(D) on X_1(2,N), for divisors supported on noncuspidal places }
    supp,mults := Support(D);
    answer := DivisorGroup(Curve(X))!0;
    for i in [1..#supp] do
        answer +:= mults[i]*HeckeOperatorFast(X,q,supp[i]);
    end for;
    return answer;
end intrinsic;

intrinsic PlacesUpToDiamondFast(X::MDCrvMod11, S::SeqEnum[PlcCrvElt]) -> SeqEnum[PlcCrvElt]
{ One representative per diamond orbit, using simultaneous isomorphism tests on P and Q }
    N := Level(X); F := BaseRing(X); pF := #F;
    data := [* *];
    for x in S do
        E := EllipticCurve(X,x); Ls := LevelStructure(X,x);
        Append(~data,<E,Ls`P,Ls`Q,MinimalPolynomial(jInvariant(E),F),Degree(x)>);
    end for;
    groups := AssociativeArray();
    for i in [1..#S] do
        key := <data[i][4],data[i][5]>;
        if IsDefined(groups,key) then Append(~groups[key],i); else groups[key]:=[i]; end if;
    end for;
    reps := [];
    for key in Keys(groups) do
        idx := groups[key]; done := [false : i in idx];
        for a in [1..#idx] do
            if done[a] then continue; end if;
            done[a] := true; Append(~reps,S[idx[a]]);
            E1 := data[idx[a]][1]; P1:=data[idx[a]][2]; Q1:=data[idx[a]][3];
            K1 := BaseRing(E1); k := data[idx[a]][5];
            mults := {@ <c*P1,c*Q1> : c in [1..N-1] | GCD(c,N) eq 1 @};
            for b in [a+1..#idx] do
                if done[b] then continue; end if;
                E2:=data[idx[b]][1]; P2:=data[idx[b]][2]; Q2:=data[idx[b]][3]; K2:=BaseRing(E2);
                if K2 cmpne K1 then
                    ok := true;
                    try
                        ainv := [K1!c : c in aInvariants(E2)];
                        pc := [K1!c : c in Eltseq(P2)]; qc := [K1!c : c in Eltseq(Q2)];
                    catch e
                        ok := false;
                    end try;
                    if not ok then
                        Embed(K2,K1);
                        ainv := [K1!c : c in aInvariants(E2)];
                        pc := [K1!c : c in Eltseq(P2)]; qc := [K1!c : c in Eltseq(Q2)];
                    end if;
                else
                    ainv:=aInvariants(E2); pc:=Eltseq(P2); qc:=Eltseq(Q2);
                end if;
                found := false;
                for i in [0..k-1] do
                    pe := pF^i;
                    E2i := EllipticCurve([c^pe : c in ainv]);
                    P2i := E2i![c^pe : c in pc]; Q2i := E2i![c^pe : c in qc];
                    iso,phi := IsIsomorphic(E2i,E1);
                    if not iso then continue; end if;
                    for alpha in Automorphisms(E1) do
                        if <alpha(phi(P2i)),alpha(phi(Q2i))> in mults then found:=true; break; end if;
                    end for;
                    if found then break; end if;
                end for;
                if found then done[b]:=true; end if;
            end for;
        end for;
    end for;
    return reps;
end intrinsic;

