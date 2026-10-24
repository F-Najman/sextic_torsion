/* Independent global identity guard for the already-running N30 filter.
   The assertions hold in the universal function field, so at every
   nonsingular noncuspidal specialization in that run. */
SetColumns(0);SetMemoryLimit(2*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
X:=MDX11(2,30,GF(7));
assert X`_E[1] eq 0 and X`_E[3] eq 0;
print "M2_30_NORMALIZATION_GUARD_DONE global_a1_zero=true global_a3_zero=true",Cputime(t0);
quit;
