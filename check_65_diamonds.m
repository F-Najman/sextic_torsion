/* Lemma 9.4: the diamond operator <2> on X_1(5) in the coordinate t of the
   Tate normal form b=c=t.  Over Q it is t -> -1/t, so over F_2 it is t -> 1/t.
   As <3>=<-2>, the operator <3> acts in the same way.  Together with the action
   z -> z+1 of <2> on X_1(13)/<3> over F_2 (marker FACTOR13_AF_INVOLUTION in
   logs/cartier/X1_13_twist_parameter.log), this identifies the automorphisms
   tau_T and tau_U of Y = X_1(65)/<4> as diamond operators. */
SetColumns(0);
AttachSpec("mdmagma/v2/mdmagma.spec");
for K in [* Rationals(), GF(2) *] do
    F<t>:=FunctionField(K);
    E:=EllipticCurve([1-t,-t,-t,0,0]);
    P:=E![0,0];
    assert P ne E!0 and 5*P eq E!0;
    bc2:=TateNormalForm_bc(E,2*P);
    bc3:=TateNormalForm_bc(E,3*P);
    assert bc2 eq [-1/t,-1/t] and bc3 eq bc2;
    printf "X1_5_DIAMOND char=%o <2>,<3>: t -> %o\n",Characteristic(K),bc2[1];
end for;
print "CHECK65_DIAMONDS_DONE";
quit;
