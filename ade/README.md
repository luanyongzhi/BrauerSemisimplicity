# Exact ADE certificates

From the cloned repository root, run:

```sh
sage -python ade/run.py
```

The default command reconstructs all three types and fails on an assertion,
nonzero child exit status, changed source checksum, or unexpected result.
Assertions must be enabled. It needs SageMath, NumPy and SymPy; no network,
GAP, Julia or separate python-flint installation is used. SageMath 10.1 was
tested. Expected total runtime on the test computer is about one minute.

The runner resolves paths from its own location, so the current working
directory is immaterial. All constructors and their outputs run in disposable
scratch space. No supplied expected result or source file is rewritten. To
retain a concise actual report, choose a separate path:

```sh
sage -python ade/run.py --report ade-results.json
sage -python ade/run.py --only e6
sage -python ade/run.py --only e7
sage -python ade/run.py --only e8
```

Do not run the legacy constructors directly in `ade/code`: some write result
files next to themselves or in their working directory. `run.py` supplies the
required order and isolation. Expected values in `expected/results.json` are
comparison targets only; constructors never read them. Their hashes of full
branching/projector data use canonical JSON (sorted keys, no extra whitespace).

## What is certified

| Paper input | Constructor/check | Certificate |
| --- | --- | --- |
| E6 singleton spectra | `e6_spectrum.py` | 72-root geometry; all 25 class fingerprints; the actual perpendicular A5 embedding; 25 by 11 branching and all eleven determinant factorizations |
| E6 higher layers | `e6_spectrum.py` | Fresh polynomial matrices of orders 270, 270, 540, 135; formal relations and root conjugations; exact rational ranks; order-two jet of size 1080; monic determinant degree exhaustion |
| E7 singleton spectra | `e7_singleton.py`, `e7_fusion_check.py` | Geometric source classes and all 60 target-class identifications; 60 by 34 aggregated branching; exceptional set |
| E7 higher collision inputs | `e7_pair_branching.py`, `e7_higher_branching.py` | Actual stabilizers, coefficient actions, annihilator checks and branching |
| E7 remaining parameter -7 | Three `e7_*` residual scripts | Full ranks 945, 2835, 945 modulo 1009, plus formal relation checks; full modular rank proves characteristic-zero invertibility |
| E8 singleton/pair inputs | `e8_singleton.py`, `e8_pair_branching.py`, `e8_exact_inputs.py` | Actual 240-root class fusion; all 112 fingerprints distinct; branching; unbounded-integer character-product cross-checks |
| E8 four-root layer and screen | `e8_four_branching.py`, `e8_final_screen.py` | 9450-element orbit, exact normalizer order 73728, coefficient quotient and branching; residual values -3 and 5 |
| E8 pair action | `e8_relations.py` | 90,720 bridge comparisons, 158 formal defining relations and 960 root-conjugation identities |
| E8 six remaining channels | `e8_annihilator.py` | Four integral projector certificates covering six W-types; nonzero values -5, -20, -20 at -3 and -27, -12, -12 at 5 |

For E6, exact nullities at candidate roots give lower bounds on determinant
multiplicities. The standard-pair matrix has nullity 406 at 1 and order-two
jet nullity 466. Together with the other roots these bounds exhaust its degree
1080. The other three monic higher determinants likewise exhaust their degrees.
This is a characteristic-zero factorization certificate, not a deduction from
modular rank deficiency.

E7 combines the three split D6 character pairs only where their central
reflection scalars agree. This suffices for the obstruction union; it is not
a fully separated 60 by 37 table. E8 normalizer generation uses a fixed-seed
Schreier sample only to propose generators, then verifies exact closure order.
The projector construction divides by positive integer gcds after each factor;
its final integral blocks are primitive. The four base blocks are respectively
128 I9, 64 I9, -528 I16, -528 I16.

## Scope and provenance

The programs certify the finite inputs used by the paper. The paper supplies
the normal-form/cellular results, specialized-layer reduction, central-element
argument, perpendicular contraction and treatment of the maximal layers. A
successful run does not replace those proofs. In particular, the E8 screen is
conditional on the algebraic reductions stated there.

The three `weyle*.jl` files are pinned CHEVIE source data. The scripts evaluate
their character values at q=1, check orthogonality and dimensions, and construct
the actual subgroup fusion geometrically. They do not independently derive the
irreducible character tables from first principles. Julia is not needed.

`PROVENANCE.json` records source and distributed SHA256 hashes and every
adaptation. E6 uses Sage rational-matrix rank in place of python-flint integer
rank. E7's standard residual matrix uses Sage finite-field rank. The E8 checker
retains the E8 part of an archived mixed checker, removing unrelated Type-C
calculations and host-specific wrappers. Geometric and certificate assertions
are retained. No previous paper, previous output table, or remote source is read
as mathematical input at runtime.

CHEVIE is distributed under the retained `code/CHEVIE_LICENSE` (MIT). The
author-supplied research scripts arrived without an explicit license; this
package does not assign them a new one. The authors should choose the license
for their research code before publication. That choice must preserve the
separate CHEVIE notice. `TESTED.json` records the actual packaged checks.
