SetColumns(0);SetMemoryLimit(10^9);t0:=Cputime();
A:=JOne(36);assert Dimension(A) eq 17;
f:=FrobeniusPolynomial(A,5);assert Degree(f) eq 34;
z:=Integers()!Evaluate(f,1);assert z mod (12*252*7812) eq 0;
print "JAC36_F5_ORDER",z,Factorization(z),"index_Q",z div (12*252*7812),Cputime(t0);
print "JAC36_F5_DONE";quit;
