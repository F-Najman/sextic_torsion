# sextic_torsion

Code and logs accompanying the paper *Classification of torsion of elliptic curves over sextic
fields*, which proves

  Φ(6) = Φ^∞(6) ∪ {Z/25, Z/37, Z/4 × Z/12, Z/2 × Z/24}.

The table in [Which script proves what](#which-script-proves-what) is the authoritative map from
each statement of the paper that rests on a computation to the script that verifies it and the log
that script produced. In the digital version of the paper, the words reporting a computation
("we compute", "we check", "we verify") are links to the script that carries it out. Most logs can
be regenerated with [`verify_all.sh`](#reproducing-the-computations). A few cannot be regenerated
by a single command, because the run was interrupted and then continued by hand inside the same
Magma session; they are marked "archived" in `jobs.txt`, and the section
[Archive-only logs](#archive-only-logs) explains each of them. The logs are those of the original
runs, byte for byte; see [Provenance of the committed logs](#provenance-of-the-committed-logs).

## Getting the code

The Magma code depends on [mdmagma](https://github.com/koffie/mdmagma), which supplies the
moduli-theoretic routines (`MDX1`, `MDX11`, `HeckeOperator`, `DiamondOperator`, `ModuliPoint`,
`CuspOrbitsQ`) and, through its own `Magma` submodule, Sutherland's optimised models of
`X_1(N)` and `X_1(M,N)`. The point counts use Sutherland's
[GL2 package](https://github.com/AndrewVSutherland/Magma). Both are pinned as submodules.

```bash
git clone https://github.com/F-Najman/sextic_torsion
cd sextic_torsion
git submodule update --init                       # mdmagma and gl2
git -C mdmagma submodule update --init Magma      # Sutherland's models, used by mdmagma
git -C mdmagma apply ../mdmagma-magma2.30.patch   # needed with Magma V2.30, see below
```

Do not use `git clone --recurse-submodules`: mdmagma has a further submodule
`tests/magma-unittest` registered under an SSH URL, so a full recursive clone can fail without a
GitHub SSH key. It is not used here.

The pinned commits are mdmagma `d1c019c6` (with its `Magma` submodule at `6a5a68a8`) and the GL2
package at `e70690a`, exactly as in
[quintic-torsion](https://github.com/F-Najman/quintic-torsion). The latter two are different
checkouts of the same repository on purpose: mdmagma pins a commit that predates `gl2points.m`,
so the point-count scripts attach the separate copy in `gl2/`, and every other script attaches
mdmagma's.

**The mdmagma patch.** Magma V2.30 rejects two statements of the pinned mdmagma: an
`end intrinsic` without its semicolon in `v2/curves/divisors.m` and in
`v2/modular_curves/base_curve.m`. All computations here were run with a copy of mdmagma in which
these two semicolons had been added, and nothing else changed. `mdmagma-magma2.30.patch` makes
exactly this change; apply it with the command above. `verify_all.sh` checks whether it is
applied and applies it if not. With Magma V2.29 the patch is harmless.

## Running the computations

All scripts are run from the repository root; the `Attach`, `AttachSpec` and `load` paths and the
paths `data/...` are resolved from there. Everything was run with Magma V2.30-1. The Python
scripts need Python 3; `candidates.py`, `verify_c3x18_arithmetic.py` and
`verify_a3_n30_arithmetic.py` also need SymPy (`python3 -m pip install sympy`), the others only
the standard library. The shell driver requires Bash 4.4 or later and GNU command-line tools.
Magma commands use `-n` to suppress personal startup files.

Several drivers take parameters, passed as `magma -n -b Name:=value script.m`:

* `hecke_sieve_deg5_v6.m` (the driver of the quintic paper, here with `D:=6`),
  `hecke_sieve_cached.m`, `hecke_sieve_cyclic_commonzero.m` — the Hecke sieve on `X_1(N)`, with
  parameters `N`, `p`, `Qs`, `D`, `M`, `MemGB`, `Shard`, `NShard`. The last two also accept
  `FromDegree` and `MaxPart`, which restrict the run to the partition types of degree at least
  `FromDegree` and largest part at most `MaxPart`; such a run ends in `SIEVE_RESUME_DONE`
  rather than `SIEVE_DONE`. The three differ only in caching and in the Hecke routine.
* `sieve_M2_strong.m`, `sieve_M2_fast_cached.m`, `sieve_M2_tate_cached.m` — the same sieve on
  `X_1(2,N)` (`M:=2`), with a congruence condition on `q` computed respectively from the
  strong cyclotomic condition, the same with fast routines, and the exact Tate-cusp
  permutations in `data/cusp_conditions_m2.m`.
* `analyze_survivors.m` — parameters `N`, `p`, `D`, `SURV`.
* `torsion_over_Fq.m` — parameters `Ns`, `Ps`, `K`.
* `cusp_double_rr_local_ref_lazy.m`, `cusp_double_rr_3n.m`, `cusp_double_rr_3n_localref.m` —
  parameters `N`, `D`, `DirectCheck`, `UseJets`, `Limit`.
* `sieve_42_fixed_shard.m` — parameter `Shard` in `0..3`.
* `rational_two_torsion.m`, `two_torsion_49_T2.m` — `MatrixFile` (and `N`, `Qs`) to export the
  integral matrices that the Python replays read.
* `cusp63_direct_cartier.m`, `cusp51_direct_cartier.m`, `cusp51_direct_rr_kernels_full.m` —
  parameters `D`, `DirectCheck`, `StopAfterExceptions`, `MemGB`, `Limit`, `JetFile`.

Everything else runs without parameters. The commands that produced the committed logs are in
`jobs.txt`. For example:

```bash
magma -n -b N:=35 p:=3 Qs:=11,13 D:=6 M:=1 MemGB:=2 hecke_sieve_deg5_v6.m
magma -n -b N:=45 D:=6 DirectCheck:=3 UseJets:=1 cusp_double_rr_local_ref_lazy.m
magma -n -b Shard:=0 sieve_42_fixed_shard.m
python3 replay_fi_jets.py data/fi63_direct_cartier_jets.m
python3 candidates.py
```

Magma may buffer its output, so a terminated run can leave an empty or incomplete log, and some
Magma errors return exit status 0. Check both the exit status and the closing marker: every
script prints a marker when it finishes, and the third column of `jobs.txt` gives it for each
log. Some scripts end without `quit;`; run them with standard input at end of file
(`< /dev/null`), as `verify_all.sh` does, so that Magma exits when the script ends.

## Reproducing the computations

`jobs.txt` lists each log file, the command that produces it, and its closing marker. An optional
fourth column gives an expected result as an extended regular expression.

```bash
./verify_all.sh                 # every job, Magma jobs one at a time
./verify_all.sh -j 4            # Magma jobs four at a time
./verify_all.sh X1_42           # only the jobs whose log path contains X1_42
./verify_all.sh -j 4 counts     # both
python3 -m unittest discover -s tests   # the driver's tests; no Magma needed
git diff                        # what the rerun changed
```

The driver overwrites the selected logs in place; `git diff` compares the new output with the
committed run. It runs the selected Magma jobs first, then the Python jobs one at a time and in
the order of `jobs.txt`, because the Python jobs read the logs and data files written by the
Magma jobs and by one another (`candidates.py` reads almost everything). Python-only selections
run without Magma or its submodules. When Magma jobs are selected, the driver first applies the
mdmagma patch if necessary and then attaches every package once, in separate processes, so that
parallel jobs do not compile the same signature files at the same time.

A job passes only if its command exits 0, its closing line is the expected marker (optionally
followed by a space and details), and any result expression matches. The closing line is the
last non-blank line of the log, ignoring carriage returns, a trailing `Script done on ...` line
written by `script(1)` (four original runs were recorded on a terminal), and trailing lines
that consist of a CPU time only (Magma prints the time on a line of its own after a sequence).
Failures are listed at the end and give a nonzero exit status. Missing or malformed job lists,
invalid options, and filters matching no runnable job also give a nonzero status.

Rerunning all jobs costs about 6 hours of CPU time, most of it on the shards of
`X_1(42)`, the resumed sieves on `X_1(33)`, `X_1(2,22)` and `X_1(49)`, the doubled-cuspidal test
on `X_1(45)` and the class sieve on `X_1(36)`. Some jobs rewrite data files in `data/` (the
matrix exports, the survivor lists of the `Z/2 × Z/30` filter, the coverage records); a correct
rerun rewrites them with the same content.

What a successful reproduction should change, and nothing else:

* timing figures, including the `time=...` fields, the `Time:` lines and the bare CPU times;
* the transcript lines of `logs/tests/probe_42_hecke_hash.log` and `logs/n30/other_quotients.log`,
  which were recorded on a terminal (`Script started on`, `Script done on`, carriage returns),
  and the path in Magma's `Loading "..."` lines;
* whitespace, and whatever a different version of Magma prints differently, such as the printed
  form of a model or an ideal, or the order of the elements in a printed set;
* in `data/sieve_42_coverage.json`, the field `shard_source_sha256` (see
  [Provenance](#provenance-of-the-committed-logs)).

Everything else is a genuine discrepancy: the place counts on the `NONCUSP` and `REPS` lines, the
Hecke primes on the `QLIST`, `CC` and `TATE_CUSP_CONDITION` lines, every `RESULT`, `SURVIVOR`
and `SUMMARY` line, the group invariants and class-membership results, the dimensions and ranks
of the Riemann–Roch and Cartier tests, and the presence of the closing marker.

### Archive-only logs

These logs are kept because the proof uses them, but no single command regenerates them. Each is
listed in `jobs.txt` with the command `archived: ...`, and is not rerun by `verify_all.sh`.

* `logs/sieve/X1_33_p5_prefix.log`, `logs/sieve/X1_2_22_p3_prefix.log`,
  `logs/sieve/X1_49_p3_prefix.log` — sieve runs interrupted during degree 6. Each partition type
  completed in the prefix is recorded by a `RESULT` or `TYPE ...: no divisors` line; the rest were
  run separately (`..._resume.log`, or `X1_49_p3_A2.log`). `combine_sieve33.py`,
  `combine_sieve2_22.py` and `combine_sieve49.py` check that the union covers all 29 partition
  types of degree at most 6 and write the combined logs `X1_33_p5.log`, `X1_2_22_p3.log`,
  `X1_49_p3.log`. To redo the whole computation in one run, use the command of the prefix
  without the interruption, writing to a new log, e.g.
  `magma -n -b N:=33 p:=5 Qs:=7 D:=6 M:=1 hecke_sieve_cached.m > X1_33_p5_full.log`; it should end
  in `SIEVE_DONE` with `total_survivors=0`.
* `logs/sieve/X1_42_p5_prefix_session.log` — the terminal transcript of the `X_1(42)` sieve
  (`hecke_sieve_42_roots.m`), which completed the types (3), (4), (5) and was interrupted during
  (6); `export_42_fixed_reps.m` was then loaded into the same session. It wrote the 392 sextic
  diamond representatives to `data/c42_fixed_representatives.m`, checking that each one
  reconstructs its place, and tested the 8 divisors of type (3,3). The four shards
  `sieve_42_fixed_shard.m` test the 392 representatives from that file, and
  `combine_sieve42.py` checks that every index has a witness and writes `X1_42_p5.log`.
* `logs/cartier/Y63_session.log`, `logs/cartier/Y63_continuation.log` — the `Z/63` certificate.
  `cusp63_direct_cartier.m` computed the model and the first 18 of the 22 cusp blocks and was
  interrupted in block 19; `cartier_direct_63_diagnostic_continuation.m`, which repeats the code
  of the script with additional diagnostic output, computed blocks 19–22 and tested all 6542
  divisors in the same session. The standalone
  `magma -n -b JetFile:=data/fi63_new.m cusp63_direct_cartier.m` does the whole computation.
* `logs/cartier/Y51_session.log`, `logs/cartier/Y51_cusp_union.log`,
  `logs/cartier/Y51_rr_kernels.log` — the `Z/51` certificate, a chain of three steps in one
  session (model from `cusp51_direct_cartier.m`, cusps from
  `cartier_direct_51_union_continuation.m`, restriction kernels from
  `cusp51_retained_rr_kernels_continuation.m`). Only the first four lines of
  `Y51_session.log` are used; that log also records exploratory continuations which ended in
  errors (an assertion that each base cusp has exactly two places above it, which is false, and
  an undeclared identifier). They supply nothing to the certificate. The standalone script
  `cusp51_direct_rr_kernels_full.m` combines the three steps; it was started once but stopped
  before completion, so no log of it exists. The independent script
  `independent/cartier_51_63.m` repeats the whole computation in one run (see
  [Independent recomputations](#independent-recomputations)).
* `logs/n30/first_quotient_prefix.log`, `logs/n30/first_quotient_mixed_continuation.log` — the
  first quotient filter for `Z/2 × Z/30`. The run of `m2_30_quotient_prefilter.m` completed all
  five single-place types; two assertions were then added to the source while Magma was still
  reading it, which shifted its read position and produced the parser error
  `User error: bad syntax` at the end of the prefix log. `m2_30_mixed_continuation.m`, loaded into
  the same session, ran the five mixed types and wrote `data/m2_30_quotient_survivors.m`.
  `audit_m2_30_quotients.py` checks that the parser error is the only error and that the two
  logs together cover all ten support types. The script in this repository is the corrected
  source, with the two assertions; `magma -n -b m2_30_quotient_prefilter.m` runs the whole filter.
* `logs/finite_fields/C4xC24_hasse.log` — produced by a short Python computation that was not
  archived: the Hasse intervals over `F_5`, `F_25`, `F_125` contain no multiple of 96. The paper
  proves this by hand.

## Which script proves what

Numbering refers to the manuscript (version of 24 October 2026), and the rows follow its order.

| Statement | Script | Log |
|---|---|---|
| **Section 2** (after Lemma 2.3), as used in Propositions 4.3, 5.3, 5.5 and 8.5. No elliptic curve over `F_{2^e}`, `e ≤ 6`, has a point of order 45, 49, 51, 55, 75, 77, 81 or 85 (order 35 does occur over `F_64`) | `torsion_over_Fq.m` | `logs/finite_fields/torsion_F2_49_51_55_75_77_85.log`, `logs/finite_fields/torsion_F2_35_45_49_51_81.log` |
| **Section 2 and Sections 4–9.** Every value of Abramovich's bound and every order of a diamond operator quoted in the paper; the subgroup generated by 12 and 56 in `Δ(247)` has order 9 | `global_table.py` | `logs/global_table.log` |
| **Lemma 3.1.** The 66 minimal groups, by a short direct enumeration | `group_list.py` | `logs/group_list.log` |
| **Lemma 3.1 and Table 1.** The bookkeeping of the proof: 66 minimal groups outside the 48 known ones, 57 after the `Z/37` case, each excluded by one of the computations below | `candidates.py` (full subgroup enumeration and the ledger of exclusions, which reads the logs below); `audit_classification_glue.py` (independent enumeration through maximal subgroups) | `logs/candidates.log`, `logs/classification_glue_audit.log`; `data/candidates.json` |
| **Section 4; Propositions 4.1, 4.3, 8.5, 9.1, 9.3 and 9.5.** The weight 2 newforms of positive analytic rank at every level used (all ranks recorded as proved), and the values of their characters at the elements used | `lmfdb_ranks.py` (offline, from the archived LMFDB records; `--fetch` refetches them) | `logs/lmfdb_ranks.log`; `data/lmfdb_newforms_cache.json` |
| **Proposition 4.1 and Table 2.** The global method for `n` = 119, 121, 125, 143, 169, 187, 209, 221, 289, 323, 361: the inequality `6k_3 < 325·i(n)/2^15`, the order of `⟨a⟩` not dividing 6, and the least CM degree exceeding 6 | `global_table.py`; the same checks, in exact arithmetic, in `candidates.py`. The finiteness of `(⟨a⟩−1)J_1(n)(Q)` is taken from the quintic paper | `logs/global_table.log`, `logs/candidates.log`; `data/lmfdb_ranks.json`, `data/lmfdb_ranks_summary.txt`, `data/char_kernels_all.log` |
| **Proposition 4.3 and Table 3.** `Z/91`, `Z/95`, `Z/133`, `Z/247`: the character values used with several diamond operators (`χ(31) = −1` at 95; at 133, `χ_f(6) = χ_s(3) = 1`, with the orders 9 and 18 of `⟨3⟩` and `⟨6⟩` in `Δ(133)` from `global_table.py`; `⟨12, 56⟩` of order 9 at 247) | `check_diamonds.m`, `global_table.py` | `logs/check_diamonds.log`, `logs/global_table.log` |
| **Lemma 5.2.** `Z/7 × Z/7`, `Z/9 × Z/9`: every simple new factor of `⊕ S_2(Γ_0(m^2), ε)`, `ε` even mod `m`, has winding dimension equal to its dimension (m = 7, 9) | `full_level_rank.m` | `logs/full_level_rank.log` |
| **Proposition 5.3.** `Z/3 × Z/27`: no elliptic curve over `F_4`, `F_16`, `F_64` contains `Z/3 × Z/27` (proved by hand in the paper) | `noncyclic_finite_fields.m` | `logs/finite_fields/C3xC27_F4_F16_F64.log` |
| **Propositions 5.3 and 5.5.** The same for `Z/3 × Z/15`, `Z/3 × Z/21`, `Z/3 × Z/27` | `check_local_3n.m` | `logs/finite_fields/C3xCn_F4_F16_F64.log` |
| **Proposition 5.3.** `Z/4 × Z/24`: no good reduction above 5 | none: the paper proves this by hand from the Hasse bound over `F_5`, `F_25` and `F_125`; the archived log, whose script was not kept, records the same check | `logs/finite_fields/C4xC24_hasse.log` |
| **Proposition 5.5.** `Z/45`: the closed cusps of `X_1(45)` over `F_2` of degree at most 6 have degrees `1^12, 2^4, 4^5, 6^2`, and `dim L(2C) = 1` for all 19048 effective cuspidal divisors `C` of degree 6 | `cusp_double_rr_local_ref_lazy.m` | `logs/rr/X1_45_F2_double_cusp.log` |
| **Proposition 5.5** (validation of the method). The jet computation agrees with ordinary Riemann–Roch on all 24 degree-two cuspidal divisors of `X_1(21)` over `F_2` | `cusp_double_rr_local_ref_lazy.m` | `logs/rr/X1_21_F2_double_cusp_validation.log` |
| **Proposition 5.5.** `Z/3 × Z/21`: `dim L(2C) = 1` for all 2608 effective cuspidal divisors of degree 3 on `X_1(3,21)` over `F_4` | `cusp_double_rr_3n.m` | `logs/rr/X1_3_21_F4_double_cusp.log` |
| **Proposition 5.5.** `Z/3 × Z/15`: of the 944 effective cuspidal divisors of degree 3 on `X_1(3,15)` over `F_4`, 912 have `dim L(2C) = 1` and the other 32 are supported on the cusps of degree 1, with `dim L(C) = 1` | `cusp_double_rr_3n_localref.m` | `logs/rr/X1_3_15_F4_double_cusp.log` |
| **Proposition 5.5.** `Z/3 × Z/15`: with `σ_a(ζ_15) = ζ_15^a`, the generator `σ_7` of `Gal(Q(ζ_15)/Q(ζ_3))` and the Frobenius `σ_4 = σ_7^2` at 2 each fix exactly 16 of the 32 cusps (so the same 16) | `cusp_orbits.py` | `logs/classes/X1_3_15_cusp_orbits.log` |
| **Proposition 5.5.** `Z/3 × Z/15`: reduction at 2 is injective on `J(Q(ζ_3))_tors` (cusp group `[3,3,24,96,96]` over `F_4`, `J(F_7) = [3,3,3,24,480,480]`, equal 2-primary orders 8192) | `verify_3_15_reduction_injective.m` | `logs/classes/X1_3_15_reduction_injective.log` |
| **Proposition 5.7.** `Z/54`: on `X_1(27)` over `F_2` the non-cuspidal place counts are (0,0,0,0,0,9), one diamond orbit; `dim L(B) = 1` at 2 and at 5 for the divisor `B` of the proof; the curve at `T_2` has `E(F_25) ≅ Z/27` | `certify_54_via27.m`; independent count `count_27_gl2.m` | `logs/rr/X1_27_C54_certificate.log`, `logs/counts/X1_27_p2.log` |
| **Proposition 6.1 and Table 4.** `Z/34`, `Z/38`: the Hecke sieve at 3 (q = 5, 7, 11) leaves no survivors | `hecke_sieve_deg5_v6.m` | `logs/sieve/X1_34_p3.log`, `logs/sieve/X1_38_p3.log` |
| **Proposition 6.1 and Table 4.** `Z/35`: the Hecke sieve at 3 (q = 11, 13) leaves no survivors among the 1, 3, 11, 12 divisors of types [3], [5], [6], [3,3] | `hecke_sieve_deg5_v6.m`; place counts and fast/reference Hecke divisors checked by `verify_35.m` | `logs/sieve/X1_35_p3.log`, `logs/counts/X1_35_p3_verify.log` |
| **Proposition 6.1 and Table 4.** `Z/39`: the Hecke sieve at 2 (q = 5, 7) leaves no survivors; the fast Hecke divisors agree with the reference ones in characteristic 2 | `hecke_sieve_deg5_v6.m`; `check_char2_hecke.m` | `logs/sieve/X1_39_p2.log`, `logs/tests/check_char2_hecke.log` |
| **Proposition 6.1 and Table 4.** `Z/50`: the Hecke sieve at 3 (q = 7, 11) leaves no survivors | `hecke_sieve_deg5_v6.m` | `logs/sieve/X1_50_p3.log` |
| **Proposition 6.1 and Table 4.** `Z/40`: the Hecke sieve at 3 (q = 7, 11) leaves no survivors | `hecke_sieve_cached.m` | `logs/sieve/X1_40_p3.log` |
| **Proposition 6.1 and Table 4.** `Z/33`: the Hecke sieve at 5 (q = 7) leaves no survivors; the interrupted run and the resumed run together cover all 29 partition types | `hecke_sieve_deg5_v6.m` (prefix), `hecke_sieve_cached.m` (resume), `combine_sieve33.py` | `logs/sieve/X1_33_p5_prefix.log`, `logs/sieve/X1_33_p5_resume.log`, `logs/sieve/X1_33_p5.log`; `data/sieve_33_coverage.json` |
| **Proposition 6.1, Table 4 and Remark 6.2.** `Z/42`: the Hecke sieve at 5 (q = 11, 13, 17) leaves no survivors among the 2, 33, 108, 8 and 392 divisors of types (3), (4), (5), (3,3), (6) | `hecke_sieve_42_roots.m` and `export_42_fixed_reps.m` (types (3), (4), (5), (3,3)), `sieve_42_fixed_shard.m` (four shards of the 392 sextic representatives), `combine_sieve42.py`; orbit hashing checked by `probe_42_hecke_hash.m` | `logs/sieve/X1_42_p5_prefix_session.log`, `logs/sieve/X1_42_p5_shard{0,1,2,3}.log`, `logs/sieve/X1_42_p5.log`, `logs/tests/probe_42_hecke_hash.log`; `data/c42_fixed_representatives.m`, `data/sieve_42_coverage.json` |
| **Proposition 6.1 and Table 4.** `Z/2 × Z/22` (hence `Z/44`): the Hecke sieve at 3 (q = 5, 7) leaves no survivors; prefix and resume cover all 29 types | `sieve_M2_strong.m` (prefix), `sieve_M2_fast_cached.m` (resume), `combine_sieve2_22.py`; the `X_1(2,N)` routines of `fast_hecke_11.m` are checked in `test_fast_hecke_11.m`, `test_fast_hecke_22.m`, `test_fast_diamonds_22.m` | `logs/sieve/X1_2_22_p3_prefix.log`, `logs/sieve/X1_2_22_p3_resume.log`, `logs/sieve/X1_2_22_p3.log`, `logs/tests/test_fast_hecke_{11,22}.log`, `logs/tests/test_fast_diamonds_22.log`; `data/sieve_2_22_coverage.json` |
| **Proposition 6.1, Table 4 and Remark 6.2.** `Z/2 × Z/28` (hence `Z/56`): the Hecke sieve at 3 (q = 5, 11) leaves no survivors among the 8 and 16 representatives of degrees 5 and 6; both primes satisfy the cuspidal-complement condition | `sieve_M2_tate_cached.m`; `cusp_conditions_m2.py` (Tate-cusp permutations) | `logs/sieve/X1_2_28_p3.log`, `logs/sieve/cusp_conditions_m2.log`; `data/cusp_conditions_m2.m` |
| **Remark 6.2.** The congruence condition on `q`, per curve | the sieve drivers above | the `CUSPFIELDS`, `CC` and `TATE_CUSP_CONDITION` lines of the sieve logs |
| **Proposition 6.3 and Table 4.** `Z/49`: at 3, `A_5` excludes the 1 and 3 representatives of degrees 4 and 5 (prefix run) and the 8 of degree 6 (resumed run); all other types are empty | `hecke_sieve_cached.m` (prefix), `hecke_sieve_cyclic_commonzero.m` (resumed run), `combine_sieve49.py` | `logs/sieve/X1_49_p3_prefix.log`, `logs/sieve/X1_49_p3_resume.log`, `logs/sieve/X1_49_p3.log`; `data/sieve_49_coverage.json` |
| **Remark 6.4.** `A_2 = T_2 − 1 − 2⟨2⟩` annihilates `J_1(49)(Q)`: on the integral homology of `JOne(49)` (rank 138), `ker A_5` has dimension 2 mod 2 and meets `ker T_2` trivially; `A_2` then also excludes the 8 representatives of degree 6 | `two_torsion_49_T2.m` (independent binary replay `verify_two49_T2.py`); `sieve_49_A2.m` with `hecke_two_commonzero.m` (checked by `test_hecke_two_commonzero.m`) | `logs/homology/J1_49_T2.log`, `logs/homology/J1_49_T2_replay.log`, `logs/sieve/X1_49_p3_A2.log`, `logs/tests/test_hecke_two_commonzero.log`; `data/two_torsion_49_T2_matrices.m` |
| **Proposition 6.6 and Table 5.** `Z/3 × Z/24` at 7: `A_13` is non-principal on all 12 diamond representatives of the 48 cubic places | `sieve_3_24.m`; the `M = 3` routines are checked in `test_moduli_point3.m`, `test_fast_hecke_3.m` | `logs/sieve/X1_3_24_p7.log`, `logs/tests/test_moduli_point3.log`, `logs/tests/test_fast_hecke_3.log` |
| **Proposition 6.6 and Table 5.** `Z/6 × Z/12` at 7: `A_13 = T_13 − 14` is non-principal on all 48 representatives of the 96 cubic places | `sieve_6_12.m` with `model_6_12.m`; model and Hecke checks `test_moduli_612.m`, `test_hecke_612.m` | `logs/sieve/X1_6_12_p7.log`, `logs/tests/test_moduli_612.log`, `logs/tests/test_hecke_612.log` |
| **Proposition 6.6 and Table 5.** `Z/4 × Z/16` at the prime above 3 of `Q(i)`: `A_5` is non-principal on all 80 representatives of the 320 cubic places over `F_9` | `sieve_4_16_p3.m`, `preflight_4_16_p3.m` | `logs/sieve/X1_4_16_p3.log`, `logs/tests/preflight_4_16_p3.log` |
| **Proposition 6.6 and Table 5.** `Z/4 × Z/20` at 3: `T_29` excludes the 24 quadratic representatives and `T_13` the 80 cubic ones | `sieve_4_20_p3.m`, `preflight_4_20_p3.m` | `logs/sieve/X1_4_20_p3.log`, `logs/tests/preflight_4_20_p3.log` |
| **Remark 6.7.** The place counts of Table 5, verified independently of the models: `#Y_1(3,24)(F_{7^d}) = 0, 0, 144`, `#Y_1(6,12)(F_{7^d}) = 0, 0, 288`, `#Y_1(4,16)(F_{9^d}) = 0, 0, 960`, `#Y_1(4,20)(F_{9^d}) = 0, 192, 960` | `independent_m3_counts.m`, `independent_m6_counts.m`, `independent_m4_p3_counts.m`, `independent_m4_20_p3_counts.m` | `logs/counts/X1_3_18_and_3_24_p7.log`, `logs/counts/X1_6_12_p7.log`, `logs/counts/X1_4_16_p3.log`, `logs/counts/X1_4_20_p3.log` |
| **Table 4 and Proposition 7.2.** `Z/32`: the Hecke sieve at 3 (q = 5, 7, 11) leaves 4 divisors of degree 6, none of whose classes lies in the reduction of `J_1(32)(Q) = [2,2,2,4,120,11640]` | `hecke_sieve_deg5_v6.m`; `analyze_survivors.m` | `logs/sieve/X1_32_p3.log`, `logs/classes/X1_32_p3_class_membership.log`; `data/survivors_32_p3.txt` |
| **Table 4 and Proposition 7.3.** `Z/2 × Z/26` (hence `Z/52`): the Hecke sieve at 3 (q = 5, 7) leaves two sextic representatives; the 18 rational cusps generate `[2,14,266,3990,11970] = J(Q)` in `J(F_3)`, and neither survivor's class lies in it | `sieve_M2_tate_cached.m`, `inspect_2_26_classes.m`, `match_2_26_survivors.m` | `logs/sieve/X1_2_26_p3.log`, `logs/classes/X1_2_26_p3_classes.log`, `logs/classes/X1_2_26_p3_survivor_match.log` |
| **Proposition 7.4.** `Z/36`: in `J(F_5) = [4,4,12,12,36,7812,7812]`, the cuspidal subgroup `G` is `[3,36,7812,7812]` (index 768), and no non-cuspidal part of type [3], [4], [5], [6], [3,3] (48, 144, 660, 2856, 1176 tests) has its class in `G` | `class_cuspidal_36_certificate.m`; `#J_1(36)(F_5) = 5061852020736` from the Frobenius polynomial in `jacobian36_5.m` | `logs/classes/X1_36_p5_class_sieve.log`, `logs/classes/J1_36_F5_order.log` |
| **Proposition 7.5.** `X_1(2,24)` at 5: the cuspidal subgroup generated by the 19 unit-signature blocks is `[2,4,4,120,240]`; no non-cuspidal reduction of type [3], [4], [5], [3,3] survives (4, 32, 168, 64 tests), and of type [6] exactly two of 684 survive, with the cuspidal classes `L_1`, `L_2` | `fine_2_24_cusp_reduction.m` | `logs/n24/X1_2_24_p5_cusp_reduction.log`; `data/N24survivor_cusp_lifts.txt` |
| **Proposition 7.5** (diamond orbits). The `X_1(2,N)` orbit hashing agrees with the reference diamond orbits in degrees 3 and 4 and gives 168 and 684 orbits in degrees 5 and 6 on `X_1(2,24)` over `F_5` | `test_diamond_hash_m2_direct.m` | `logs/tests/test_diamond_hash_m2_explicit.log` |
| **Proposition 7.7.** `Z/3 × Z/18`: at 7 and 13 the Hecke tests leave 36 cubic places each; the 24 unit-signature cusp blocks give the same cusp group (order 20253807) at both primes, and none of the 1296 pairs has compatible 9-multiples | `diagonal_3_18_cusp_classes_certificate.m`; independent integer replay `verify_c3x18_arithmetic.py`; counts `independent_m3_counts.m`, `independent_m3_18_p13_counts.m` | `logs/classes/X1_3_18_diagonal_classes.log`, `logs/classes/X1_3_18_arithmetic_replay.log`, `logs/counts/X1_3_18_and_3_24_p7.log`, `logs/counts/X1_3_18_p13.log`; `data/C3x18_arithmetic_certificate.txt` |
| **Proposition 7.8.** `Z/2 × Z/30` (hence `Z/60`): `J(X_1(30))(F_7) = [4,4,4,136,16320]`; of the 22706 non-cuspidal reductions at 7 (10 support types), the map to `X_1(30)` keeps 68 and the two other maps keep 10; the quotient `W` has genus 9, and the class obstruction on `W` excludes all 10 (`J(W)(F_7) = [2^6,102,10608]`, order 69249024) | `m2_30_quotient_prefilter.m` with `m2_30_mixed_continuation.m`, `m2_30_other_quotients.m`, `filter_a3_m2_30.m` with `a3_m2_30_model.m`, `a3_m2_30_frobenius.m`, `test_a3_m2_30_map.m`, `m2_30_normalization_guard.m`; audits and replays `audit_m2_30_quotients.py`, `verify_a3_n30_arithmetic.py`, `replay_a3_projection.py`, `audit_n30_final.py` | `logs/n30/*.log`; `data/m2_30_*.m`, `data/A3_m2_30_arithmetic.txt`, `data/A3_n30_arithmetic_replay_root.json`, `data/n30_final_coverage.json` |
| **Proposition 8.5.** `Z/51`: on `Y = X_1(51)/⟨16⟩` (genus 33), the closed cusps over `F_2` have degrees `1^16, 4^4`, the image of `C − 1` has dimension 30, and it restricts with rank 6 to each of the 54808 effective cuspidal divisors of degree 6; the 100 restriction kernels and all ranks are replayed independently | `cusp51_direct_cartier.m`, `cartier_direct_51_union_continuation.m`, `cusp51_retained_rr_kernels_continuation.m` (standalone: `cusp51_direct_rr_kernels_full.m`); replays `verify_rr_transport.py`, `replay_rr_kernels.py`, `replay_cartier_projection.py`, `replay_fi_jets.py`; cusp count `audit_c51_cusp_orbits.py` | `logs/cartier/Y51_session.log` (first four lines), `logs/cartier/Y51_cusp_union.log`, `logs/cartier/Y51_rr_kernels.log`, `logs/cartier/Y51_*_replay.log`, `logs/cartier/Y51_cusp_orbits.log`; `data/fi51_retained_rr_kernels.m` |
| **Proposition 8.5.** `Z/63`: on `Y = X_1(63)/⟨8⟩` (genus 49), the closed cusps over `F_2` have degrees `1^9, 2^3, 3^9, 6`, the image of `C − 1` has dimension 45, and it restricts with rank 6 to each of the 6542 effective cuspidal divisors of degree 6 | `cusp63_direct_cartier.m` with `cartier_direct_63_diagnostic_continuation.m`; replays `replay_cartier_projection.py`, `replay_fi_jets.py` | `logs/cartier/Y63_session.log`, `logs/cartier/Y63_continuation.log`, `logs/cartier/Y63_projection_replay.log`, `logs/cartier/Y63_rank_replay.log`; `data/fi63_direct_cartier_jets.m` |
| **Proposition 8.5.** `Z/63`: the character of the only newform of positive analytic rank takes the value −1 at 8; no curve over `F_64` has order divisible by 63 | `rank_quotient_63_character.m`, `finite_field_63.m` | `logs/cartier/X1_63_character.log`, `logs/finite_fields/C63_F64.log` |
| **Proposition 9.1.** `Z/57`, `Z/65`: no good reduction above 2 (over `F_64` only `j = 0`, with the stated multipliers; the cubic field of the `X_0(65)` points has two primes above 2) | `check_good2_57_65.m` | `logs/finite_fields/good2_57_65_F64.log` |
| **Lemma 9.2.** `Z/57`: `⟨2⟩` acts trivially on `J_1(57)(Q)[2]`: on the integral homology of `JOne(57)` (rank 170), `ker A_5` mod 2 has dimension 2 and is fixed by `⟨2⟩` | `rational_two_torsion.m`; independent binary replay `verify_two_torsion.py` | `logs/homology/J1_57_two_torsion.log`, `logs/homology/J1_57_two_torsion_replay.log`; `data/two_torsion_57_matrices.m`, `data/two_torsion_57_replay.json` |
| **Section 9, before Lemma 9.4.** The coordinate `z` on `X_1(13)/⟨3⟩`, the `j`-map `j_13(z)` and the invariant `A_13` over `F_2`; the diamond operator `⟨2⟩` acts on `z` by `z ↦ z+1` (marker `FACTOR13_AF_INVOLUTION`), which is used in the proof of Lemma 9.4 | `setup13_twist_parameter_certified.m` | `logs/cartier/X1_13_twist_parameter.log` |
| **Lemma 9.4 (identification of `τ_T`, `τ_U`).** On `X_1(5)`, with the Tate normal form `b = c = t`, the diamond operators `⟨2⟩ = ⟨3⟩` act by `t ↦ -1/t` over `Q` and by `t ↦ 1/t` over `F_2` | `check_65_diamonds.m` | `logs/cartier/X1_5_diamond.log` |
| **Lemma 9.4, and Proposition 9.5 (second check).** `Z/65`: `Y = X_1(65)/⟨4⟩` has genus 41 and `τ_T`, `τ_U` are commuting automorphisms of order 4; the image of `C − 1` on the summed diamond images (dimensions 35 and 28) restricts with rank 6 to each of the 1792 cuspidal divisors of degree 6 over `F_2` that are compatible with the fields of definition of the cusps of `X_1(65)` (the proof of Proposition 9.5 uses the computation on all 3576 divisors in the next row) | `cusp65_genus41_rr_transport_full.m`; replays as for `Z/51`; cusp count `audit_c65_cusp_orbits.py` | `logs/cartier/Y65_rr_transport.log`, `logs/cartier/Y65_*_replay.log`, `logs/cartier/Y65_cusp_orbits.log`; `data/fi65_genus41_rr_transport.m` |
| **Proposition 9.5.** `Z/65`: from the equations of Lemma 9.4, the summed diamond images have dimension 35, the image `U` of `C − 1` has dimension 28, and `U` restricts with rank 6 to each of the 3576 effective cuspidal divisors of degree 6 over `F_2` (the computation used in the proof) | `independent/cartier_65.m` | `independent/logs/cartier_65.log` |
| **Proposition 10.1.** `Z/37`: the 37-division polynomial of `E_0` (j = −9317) has factor degrees 6, 6, 6, 666; for each sextic factor the point on the corresponding twist has order 37, and the reductions at degree-one primes above 29 and 41 have 37 points, so the torsion is `Z/37`; the three sextic factors define one cyclic field and give one curve up to isomorphism | `verify_37.m`; `verify_sporadic_curves.m` | `logs/verify_37.log`, `logs/verify_sporadic_curves.log` |
| **Proposition 10.2.** `Z/25` occurs: van Hoeij's sextic point gives `(0,0)` of exact order 25, and good reductions above 17 and 19 have order 25. `Z/4 × Z/12` occurs: the base change of 162.d2 (Cremona label 162d1) to `Q(α, i)`, `α^3 − 3α − 4 = 0`, has torsion `[4, 12]` | `verify_known_examples.m`; `verify_sporadic_curves.m` | `logs/verify_known_examples.log`, `logs/verify_sporadic_curves.log` |
| **Proposition 10.3(1) and the last sentence.** `Z/2 × Z/24` occurs, with the stated field, `j`-invariant and conductor; `#Aut(K) = 2` and `E^σ ≅ E`; both points of the Derickx–Sutherland model give this curve, with exact torsion `Z/2 × Z/24` | `verify_sporadic_curves.m`; `example_2_24_direct.m`, `example_2_24_second_direct.m` | `logs/verify_sporadic_curves.log`, `logs/n24/example_first.log`, `logs/n24/example_second.log` |
| **Proposition 10.3(2)** (consistency check). The examples reduce mod 5 exactly to the two surviving places of Proposition 7.5 (representatives 87 and 577) | `match_n24_exact_points.m` | `logs/n24/match_exact_points.log` |
| **Proposition 10.4.** `Z/48`: the three 2-isogenous curves of `E_{2,24}` have torsion `Z/24`, `Z/2 × Z/12`, `Z/24` | `verify_sporadic_curves.m`; `example_2_24_two_isogenies.m`, `example_2_24_second_two_isogenies.m`; independent upper bounds `reduction_bounds_n24_isogenies.m` | `logs/verify_sporadic_curves.log`, `logs/n24/example_first_isogenies.log`, `logs/n24/example_second_isogenies.log`, `logs/n24/reduction_bounds_isogenies.log` |
| **Section 11.** The place counts of Table 4 (and of `X_1(36)`, `X_1(2,24)`, `X_1(27)`), verified independently of the models by counting elliptic curves with level structure (Sutherland's GL2 package) | `independent_counts.m` (`X_1(32)`, `X_1(33)`, `X_1(34)`, `X_1(38)`, `X_1(39)`, `X_1(40)`, `X_1(49)`, `X_1(50)`, `X_1(2,22)`, `X_1(2,26)`), `count_short_sieves.m` (`X_1(2,26)`, `X_1(2,28)`, `X_1(36)`, `X_1(40)`, `X_1(54)`), `verify_35.m`, `count_42_gl2.m`, `independent_m2_24_counts.m`, `count_27_gl2.m` | `logs/counts/*.log` |

The place counts of `X_1(3,18)`, `X_1(3,24)`, `X_1(6,12)`, `X_1(4,16)` and `X_1(4,20)` are checked in the
rows of Propositions 6.6 and 7.7 and Remark 6.7. The counts of `X_1(2,30)` over `F_7` used in
Proposition 7.8 are not checked independently.

## Independent recomputations

The folder `independent/` contains a second implementation of the computations on which the proof
depends most heavily. It was written independently of the code above, and it constructs its own
models of the curves from scratch: `X_1(51)/⟨16⟩` from `X_1(17)` (itself recomputed from the
Tate normal form) and the 3-division polynomial, `X_1(63)/⟨8⟩` from `X_1(9)` and the 7-division
polynomial, the genus 41 model of `X_1(65)/⟨4⟩` from the equations printed in Lemma 9.4,
`X_1(45)` from `X_1(9)` and a point of order 5, and `X_1(3,21)` over `F_4` as the 3-division
field of the universal curve over `X_1(7)`. The kernels of the restriction maps are computed
from Riemann–Roch spaces of differentials, without local expansions.

| Statement | Script | Log |
|---|---|---|
| **Proposition 5.5** (`Z/45`): `dim L(2C) = 1` for all 19048 divisors | `independent/doubled_cuspidal_45.m` | `independent/logs/doubled_cuspidal_45.log` |
| **Proposition 5.5** (`Z/3 × Z/21`): `dim L(2C) = 1` for all 2608 divisors | `independent/doubled_cuspidal_3_21.m` | `independent/logs/doubled_cuspidal_3_21.log` |
| **Remark 6.4 and Lemma 9.2**, in both conventions for the Hecke action on homology | `independent/two_torsion_57_49.m` | `independent/logs/two_torsion_57_49.log` |
| **Proposition 8.5** (`Z/51`): rank 6 on all 54808 divisors | `magma -n -b Case:=51 independent/cartier_51_63.m` | `independent/logs/cartier_51.log` |
| **Proposition 8.5** (`Z/63`): rank 6 on all 6542 divisors | `magma -n -b Case:=63 independent/cartier_51_63.m` | `independent/logs/cartier_63.log` |
| **Lemma 9.4 and Proposition 9.5** (`Z/65`): the model from the equations of Lemma 9.4, and rank 6 on all 3576 cuspidal divisors of degree 6 (this is the computation used in the proof of Proposition 9.5; the original one covers the 1792 compatible divisors) | `independent/cartier_65.m` | `independent/logs/cartier_65.log` |
| Genera, indices and numbers of cusps of the curves `X_1(m,n)` used | `independent/genera.sage`, `independent/cusps.sage` (SageMath) | none |

These logs were produced with the scripts as they stand here, run in a scratch directory (some
under earlier file names). The Magma runs for `Z/63` and `Z/65` take about 4 and 10 hours of CPU
time, respectively, almost all of it in the Riemann–Roch spaces at the cusps. Note that
in a tower of function fields Magma's `Decomposition` of the infinite places can return the same
place more than once, so the scripts remove duplicates.

## The archived data

`data/` holds the inputs that are not produced here and the outputs that later steps read:

* `lmfdb_newforms_cache.json` — the raw LMFDB records read by `lmfdb_ranks.py`.
* `lmfdb_ranks.json`, `lmfdb_ranks_summary.txt`, `char_kernels_all.log` — the rank and
  character-kernel data of the quintic paper, copied from its working directory. The first two
  are identical to `data/` of quintic-torsion; `char_kernels_all.log` differs from
  `logs/char_kernels_all.log` there only in wording (`paper's a` instead of
  `arXiv quartic choice a`), with the same numbers.
* `candidates.json`, `candidates.txt` — written by `candidates.py`, read by
  `audit_classification_glue.py`.
* `survivors_32_p3.txt` — the `SURVIVOR_DATA` expressions of `logs/sieve/X1_32_p3.log`, one per
  line, the input of `analyze_survivors.m`.
* `cusp_conditions_m2.m`, `cusp_conditions_m2.json` — written by `cusp_conditions_m2.py`,
  loaded by `sieve_M2_tate_cached.m`.
* `c42_fixed_representatives.m` — the 392 sextic diamond representatives on `X_1(42)` over
  `F_5` (SHA-256 `fa129d97…`, checked by `combine_sieve42.py`).
* `two_torsion_49_T2_matrices.m`, `two_torsion_49_T2_root_matrices.m` (two identical exports
  from two runs; `verify_two49_T2.py` checks they agree), `two_torsion_57_matrices.m` — integral
  homology matrices exported by `two_torsion_49_T2.m` and `rational_two_torsion.m`.
* `fi63_direct_cartier_jets.m`, `fi51_retained_rr_kernels.m`, `fi65_genus41_rr_transport.m` —
  the Cartier matrices, image bases, restriction kernels and jets exported by the `Z/63`,
  `Z/51`, `Z/65` certificates; the replays read them.
* `N24survivor_cusp_lifts.txt`, `C3x18_arithmetic_certificate.txt` — the lines between
  `N24_LIFT_CODE_BEGIN` and `N24_LIFT_CODE_END` of `logs/n24/X1_2_24_p5_cusp_reduction.log`,
  and between `M3N18_ARITHMETIC_CODE_BEGIN` and `M3N18_ARITHMETIC_CODE_END` of
  `logs/classes/X1_3_18_diagonal_classes.log`, respectively.
* `m2_30_quotient_survivors.m`, `m2_30_three_quotient_survivors.m`, `m2_30_a3_survivors.m`,
  `A3_m2_30_arithmetic.txt` — the divisors retained by the successive `Z/2 × Z/30` filters and
  the exported class arithmetic.
* `sieve_33_coverage.json`, `sieve_2_22_coverage.json`, `sieve_42_coverage.json`,
  `sieve_49_coverage.json`, `n30_final_coverage.json`, `C3x18_arithmetic_replay.json`,
  `A3_n30_arithmetic_replay_root.json`, `two_torsion_57_replay.json` — coverage and replay
  records read by `candidates.py` and `audit_n30_final.py`.

`original_names.tsv` lists, for every file of the original computations, its name in the
original computational record.

## Provenance of the committed logs

The computations were run in several working directories, with absolute or
directory-relative paths. The scripts here are the copies kept in the final record of the
computations, in which those paths had already been rewritten relative to one directory; apart
from paths, the scripts are identical to the versions that ran. The committed logs are those of
the original runs. The following scripts were changed further for this repository, only in file
paths, so their logs were produced by the unedited version:

* `gl2/magma.spec` instead of `sutherland-magma-current/magma.spec`: `independent_counts.m`,
  `count_short_sieves.m`, `count_27_gl2.m`, `count_42_gl2.m`, `independent_m2_24_counts.m`,
  `verify_35.m`;
* data files moved to `data/`: `sieve_M2_tate_cached.m` (`cusp_conditions_m2.m`),
  `sieve_42_fixed_shard.m` and `export_42_fixed_reps.m` (`c42_fixed_representatives.m`),
  `cusp65_genus41_rr_transport_full.m`, `cusp51_retained_rr_kernels_continuation.m`,
  `cusp51_direct_rr_kernels_full.m` (the exported jet files), `m2_30_quotient_prefilter.m`,
  `m2_30_mixed_continuation.m`, `m2_30_other_quotients.m`, `filter_a3_m2_30.m` (the survivor
  lists and `A3_m2_30_arithmetic.txt`);
* the paths of the logs and data files they read and write: `candidates.py`,
  `combine_sieve33.py`, `combine_sieve2_22.py`, `combine_sieve42.py`, `combine_sieve49.py`,
  `cusp_conditions_m2.py`, `verify_two49_T2.py`, `verify_c3x18_arithmetic.py`,
  `verify_a3_n30_arithmetic.py`, `replay_a3_projection.py`, `audit_m2_30_quotients.py`,
  `audit_n30_final.py` (which keeps the original names as the keys of the hashes it records, so
  that `data/n30_final_coverage.json` stays valid).

The scripts written together with the manuscript (`group_list.py`, `global_table.py`,
`lmfdb_ranks.py`, `verify_sporadic_curves.m`, `check_65_diamonds.m` and the folder
`independent/`) and `cusp_orbits.py`,
to which a closing line was added (its original run left no log), were run as they stand here.

Magma's `Write` appends to an existing file. So that a rerun replaces the exported data instead of
appending to the committed copy, the first write to each exported file now has the option
`Overwrite:=true`, in `cusp51_direct_cartier.m`, `cusp51_direct_rr_kernels_full.m`,
`cusp63_direct_cartier.m`, `cusp65_genus41_rr_transport_full.m`, `export_42_fixed_reps.m`,
`filter_a3_m2_30.m`, `m2_30_other_quotients.m`, `m2_30_quotient_prefilter.m`,
`rational_two_torsion.m` and `two_torsion_49_T2.m`. The original runs wrote to new files, so this
does not change their output.

Further points a reproduction should be aware of:

* `data/sieve_42_coverage.json` records in `shard_source_sha256` the SHA-256 of the copy of
  `sieve_42_fixed_shard.m` that was run, whose `Attach` and `load` paths differ from the
  script here. `candidates.py` compares this file with a fresh audit, so run
  `combine_sieve42.py` before `candidates.py` (as `jobs.txt` does); it rewrites the field.
* `fast_hecke_11.m` was extended during the computations; the archived version is its final
  state. The logs `logs/tests/test_fast_hecke_11.log`, `test_fast_diamonds_22.log`,
  `test_moduli_point3.log`, `test_fast_hecke_3.log`, `logs/sieve/X1_3_24_p7.log` and
  `logs/sieve/X1_2_22_p3_resume.log` were produced with earlier states of that file. The
  `X_1(2,N)` part used in the resumed `X_1(2,22)` run is preserved as
  `fast_hecke_11_m2_original.m`; `logs/tests/test_fast_hecke_22.log` is a rerun with the final
  version.
* The commands of the logs listed at the top of `jobs.txt` are not recorded in the original
  material; their parameters were read off the logs, and the drivers of the sieves on
  `X_1(32)`, `X_1(34)`, `X_1(38)`, `X_1(39)`, `X_1(50)`, the interrupted and resumed runs on
  `X_1(33)`, `X_1(2,22)`, `X_1(49)` were inferred from the file times. The candidate drivers
  print identical output.
* Several computations were run twice, the second time while auditing (the original record
  calls these `..._root.log`). Where both exist, the log kept here is the one that
  `candidates.py` reads; the two differ only in timings. This concerns `X1_6_12_p7.log`,
  `logs/counts/X1_6_12_p7.log`, `X1_4_16_p3.log`, `logs/counts/X1_4_16_p3.log`,
  `logs/counts/X1_4_20_p3.log`, `logs/tests/preflight_4_20_p3.log`,
  `X1_27_C54_certificate.log`, `example_first_isogenies.log`, `example_second_isogenies.log`,
  `match_exact_points.log`, `reduction_bounds_isogenies.log` and `J1_49_T2.log`; the kept
  `J1_57_two_torsion.log` is the audit run with `Qs:=5`, the first run (with the default
  `Qs`) stopped after `q = 5` with the same result.

## Hardware

The computations were run on the server Mordell at the University of Zagreb, with an AMD EPYC
9175F CPU (16 cores, 32 threads) and 384 GB of RAM, using Magma V2.30-1. The committed runs took
about 10.5 hours of CPU time in total (the CPU times recorded in the logs, counting interrupted
runs up to the interruption). The longest single computations were the `Z/51` certificate
(3787 s, including its setup), the resumed `X_1(49)` sieve (2772 s), the `Z/63` certificate
(2728 s), the `Z/2 × Z/30` first filter (2708 s), the `Z/45` test (2678 s), the interrupted
`X_1(33)` and `X_1(42)` sieves (2365 s, 2348 s) and the four `X_1(42)` shards (about 1100 s
each). Runs were given memory limits of 1-24 GB. The independent recomputations in
`independent/` are not included in these figures.
