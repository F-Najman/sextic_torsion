/* TTY continuation: preserve18 complete jet blocks, diagnose basis19+. */
assert N eq 63 and gx eq 49 and #cs eq 22 and IsDefined(localjets,<18,1>);
printf "Y63_DIAGNOSTIC_CONTINUATION preserved_blocks=18 time=%o\n",Cputime(t0);
for i in [19..#cs] do
    c:=cs[i]; maxm:=D div ds[i]; k:=Valuation(omega_ref,c);
    prec:=maxm+Maximum(0,-k)+2;
    repeat
        printf "RANK0_CARTIER_JET_START N=%o cusp=%o/%o degree=%o reference_order=%o precision=%o time=%o\n",N,i,#cs,ds[i],k,prec,Cputime(t0);
        S,ex:=Completion(AF,c:Precision:=prec); t:=S.1;
        printf "RANK0_CARTIER_COMPLETION N=%o cusp=%o time=%o\n",N,i,Cputime(t0);
        Kp,kmap:=ResidueClassField(c); Ks:=BaseRing(S);
        printf "Y63_RESIDUE_FIELDS cusp=%o same_parent=%o degree_place=%o degree_series=%o defining_place=%o defining_series=%o\n",i,Kp cmpeq Ks,Degree(Kp),Degree(Ks),DefiningPolynomial(Kp),DefiningPolynomial(Ks);
        expansions:=[];
        for bidx in [1..#ratios] do
            printf "Y63_BASIS_START cusp=%o basis=%o/%o time=%o\n",i,bidx,#ratios,Cputime(t0);
            sb:=t^k*ex(ratios[bidx]);
            Append(~expansions,sb);
            printf "Y63_BASIS_DONE cusp=%o basis=%o precision=%o time=%o\n",i,bidx,AbsolutePrecision(sb),Cputime(t0);
        end for;
        goodprec:=&and[AbsolutePrecision(s) ge maxm:s in expansions];
        printf "RANK0_CARTIER_EXPANDED N=%o cusp=%o sufficient_precision=%o time=%o\n",N,i,goodprec,Cputime(t0);
        if not goodprec then prec:=2*prec; end if;
    until goodprec;
    assert &and[IsZero(s) or Valuation(s) ge 0:s in expansions];
    rows:=[];
    for s in expansions do
        row:=[];
        for j in [0..maxm-1] do
            co:=Eltseq(Coefficient(s,j),GF(2)); assert #co eq ds[i]; row cat:=[GF(2)!e:e in co];
        end for;
        Append(~rows,row);
    end for;
    A:=Matrix(GF(2),rows);
    if #JetFile gt 0 then Write(JetFile,Sprintf("basisjets[%o]:=Matrix(GF(2),%o,%o,%o);",i,Nrows(A),Ncols(A),Eltseq(A))); end if;
    for m in [1..maxm] do
        J:=BasisMatrix(U)*Submatrix(A,1,1,gx,m*ds[i]); localjets[<i,m>]:=J;
        if #JetFile gt 0 then Write(JetFile,Sprintf("localjets[<%o,%o>]:=Matrix(GF(2),%o,%o,%o);",i,m,Nrows(J),Ncols(J),Eltseq(J))); end if;
    end for;
    printf "RANK0_CARTIER_JET N=%o cusp=%o/%o degree=%o reference_order=%o time=%o\n",N,i,#cs,ds[i],k,Cputime(t0);
end for;
if #JetFile gt 0 then Write(JetFile,"print \"FI_JET_DATA_READY\";"); end if;
VO:=VectorSpace(GF(2),gx); bad:=[]; tested:=0;
for v in vs do
    if Limit gt 0 and tested ge Limit then break; end if;
    A:=ZeroMatrix(GF(2),udim,0);
    for i in [1..#cs] do if v[i] gt 0 then A:=HorizontalJoin(A,localjets[<i,v[i]>]); end if; end for;
    assert Ncols(A) eq D; r:=Rank(A);
    if tested lt DirectCheck or r ne D then
        Ediv:=&+[v[i]*Divisor(cs[i]):i in [1..#cs]]; assert Degree(Ediv) eq D;
        V,vm:=DifferentialSpace(Ediv);
        V_original:=sub<VO|[VO!Eltseq(vm(b) @@ om):b in Basis(V)]>;
        direct_rank:=udim-Dimension(U meet V_original); assert direct_rank eq r;
        printf "RANK0_CARTIER_DIRECT N=%o index=%o jet_rank=%o direct_rank=%o time=%o\n",N,tested+1,r,direct_rank,Cputime(t0);
    end if;
    tested+:=1;
    if r ne D then Append(~bad,v); printf "RANK0_CARTIER_EXCEPTION N=%o index=%o rank=%o vector=%o time=%o\n",N,tested,r,v,Cputime(t0); end if;
    if tested mod 1000 eq 0 then printf "RANK0_CARTIER_PROGRESS N=%o tested=%o total=%o exceptions=%o time=%o\n",N,tested,#vs,#bad,Cputime(t0); end if;
    if StopAfterExceptions gt 0 and #bad ge StopAfterExceptions then print "RANK0_CARTIER_EXCEPTION_LIMIT_REACHED"; break; end if;
end for;
printf "RANK0_CARTIER_DONE N=%o D=%o tested=%o total=%o exceptions=%o complete=%o time=%o\n",N,D,tested,#vs,#bad,tested eq #vs,Cputime(t0);
quit;
