/* TTY continuation after slow generic Poles(j) call, same AF33 model. */
assert N eq 51 and gx eq 33 and Dimension(U) eq 30;
ps17:=Poles(jInvariant(E17));
assert {*Degree(p):p in ps17*} eq {*1^^8,4^^2*};
printf "Y51_BASE_CUSPS count=%o time=%o\n",#ps17,Cputime(t0);
cs_all:=[];
for idx in [1..#ps17] do
    qs:=Decomposition(AF,ps17[idx]);
    for q in qs do
        if q notin cs_all and Valuation(jaf,q) lt 0 then Append(~cs_all,q); end if;
    end for;
    printf "Y51_CUSP_DECOMPOSITION index=%o/%o raw_signature=%o total_unique_poles=%o time=%o\n",idx,#ps17,{*Degree(q):q in qs*},#cs_all,Cputime(t0);
end for;
assert #Seqset(cs_all) eq #cs_all;
assert {*Degree(c):c in cs_all*} eq {*1^^16,4^^4*};
pole_degree:=&+[-Valuation(jaf,c)*Degree(c):c in cs_all];
assert pole_degree eq 576;
printf "Y51_ALL_CUSPS pole_degree=%o signature=%o time=%o\n",pole_degree,{*Degree(c):c in cs_all*},Cputime(t0);
cs:=[c:c in cs_all|Degree(c) le D]; ds:=[Degree(c):c in cs];
printf "RANK0_CARTIER_CUSPS N=%o signature=%o time=%o\n",N,{*d:d in ds*},Cputime(t0);
if N eq 51 and D eq 6 then assert gx eq 33 and {*d:d in ds*} eq {*1^^16,4^^4*}; end if;
function WeightedVectors(ds,d)
    if #ds eq 0 then return d eq 0 select [[]] else []; end if;
    result:=[];
    for k in [0..d div ds[1]] do
        for v in WeightedVectors(ds[2..#ds],d-k*ds[1]) do Append(~result,[k] cat v); end for;
    end for;
    return result;
end function;
vs:=WeightedVectors(ds,D);
printf "RANK0_CARTIER_COUNT N=%o total=%o time=%o\n",N,#vs,Cputime(t0);
if N eq 51 and D eq 6 then assert #vs eq 54808; end if;
if #JetFile gt 0 then
    Write(JetFile,Sprintf("// Rank0 Cartier FI jets: fresh path.\nN:=%o; As:=[]; D:=%o; gx:=%o; o:=1; udim:=%o; JetFactor:=1; ds:=%o; perms:=[];\nlocaljets:=AssociativeArray();",N,D,gx,udim,ds));
    Write(JetFile,Sprintf("CartierMatrix:=Matrix(GF(2),%o,%o,%o);",gx,gx,Eltseq(CM)));
    Write(JetFile,Sprintf("ImageBasis:=Matrix(GF(2),%o,%o,%o);",udim,gx,Eltseq(BasisMatrix(U))));
    Write(JetFile,"basisjets:=AssociativeArray();");
end if;
etaAF:=Differential(SeparatingElement(AF)); omega_ref:=etaAF;
ratios:=[w/etaAF:w in B];
localjets:=AssociativeArray();
for i in [1..#cs] do
    c:=cs[i]; maxm:=D div ds[i]; k:=Valuation(omega_ref,c);
    prec:=maxm+Maximum(0,-k)+2;
    repeat
        printf "RANK0_CARTIER_JET_START N=%o cusp=%o/%o degree=%o reference_order=%o precision=%o time=%o\n",N,i,#cs,ds[i],k,prec,Cputime(t0);
        S,ex:=Completion(AF,c:Precision:=prec); t:=S.1;
        printf "RANK0_CARTIER_COMPLETION N=%o cusp=%o time=%o\n",N,i,Cputime(t0);
        expansions:=[t^k*ex(f):f in ratios];
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
