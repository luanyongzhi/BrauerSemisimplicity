Folded finite computations for Brauer-algebra semisimplicity

Run from this directory with SageMath 10.1:

    sage -python verify_all.py

To run selected types or display every computed certificate:

    sage -python verify_all.py --only F4 H4
    sage -python verify_all.py --details

The scripts also work when invoked by their path from another directory.
Sage supplies Python, NumPy and SymPy; no extra packages or downloaded data
are required. The tested versions are Sage 10.1, Python 3.11.1, NumPy 1.23.5
and SymPy 1.12. The complete run took about 8 seconds on the authors' machine.
Do not enable Python's -O option: the mathematical checks use assertions.

A successful run exits with code 0 and prints four PASS messages to stderr
and a JSON summary to stdout. An assertion or summary mismatch causes a
nonzero exit. All calculations run before comparison with expected.json.
The entrypoint does not write files or update expected.json. Redirect stdout
to a separate file if a local record is wanted. Floating-point arithmetic
is not used.

What is verified

G2: The two middle cells and the minimal cell, all of dimension 3, are
constructed explicitly over Q[delta]. The script checks 45 presentation
identities, 12 contravariance identities and all three determinants.
The indexing is nu_1=1, nu_2=3; the ordered triple-edge pair is (i,j)=(2,1).

H3: Starting with A5 permutations, the script constructs four 15-dimensional
cells and one 5-dimensional cell. It checks 145 presentation identities,
30 contravariance identities and five symbolic determinants. Here t=delta^2
and e_i^2=t*e_i. The source documents the corrections to the final letters
in the two five-bond relations of the cited preprint.

F4: f4_model.py reconstructs the positive E6 roots, the fixed boundary sets,
the coefficient actions and all 14 nongroup Gram forms. The 868 polynomial
presentation/contravariance identities have entry degree at most 6. Checking
each at the seven integers -2,...,4 therefore proves them. A conservative
bound of 55,296,000 on intermediate integer products proves that the NumPy
int64 arithmetic in these checks does not overflow. Each determinant has
degree at most N*d, where N is its matrix order and d its entry-degree bound.
Its exact Sage integer determinant is compared with the displayed formula
at 0,...,N*d: 650 evaluations in total. This proves all 14 polynomial
identities, including constants. In the standard S3 coefficient module,
the metric is [[2,1],[1,2]]; its determinant contributes the factors 3^12
for the two standard singleton forms and 3^30 for the pair standard form.
Labels 1, 2, 4 indicate boundary sizes. Singleton labels are B3 bipartitions;
the remaining labels are the trivial, sign and standard S3 modules, followed
by the scalar maximal layer. The variable z in this check is delta.

H4: The scripts rebuild doubled E8 roots, the full 14,400-element signed
root action, 60 root pairs and 75 frames. They check the basepair normalizer
and its A5 quotient, all choices in the transport and maximal corner rules,
42 presentation identities on each layer (5,670 column identities), and
73,800 group-ring contravariance entries. The coefficient representations
are constructed from A5 permutations, its standard module, its exterior
square, and its action on the six Sylow 5-subgroups. The pair-form spectra
are proved by squarefree annihilators and exact Lagrange projector traces.
A commuting operator with square 5I verifies equal multiplicities on the
two 3-dimensional A5 constituents. The maximal graph spectrum and the
resulting degree-150 Gram determinant are verified exactly. All matrix
products used in the spectral certification have checked integer norm
bounds. Here t=delta^2. Pair determinant factors are given up to a nonzero
coefficient-metric constant; the maximal determinant is given exactly.

Scope

These calculations establish the finite identities used in the paper.
The published basis theorems, the cellular normal forms, and the identification
of the constructed pairings with the cellular forms are mathematical inputs
explained in the paper; they do not follow merely from checking a collection
of representation matrices. The classical and general dihedral arguments
are analytic and do not require additional exploratory scripts here.

Files

verify_all.py          portable entrypoint and expected-summary comparison
verify_g2.py           G2 symbolic verification
verify_h3.py           H3 permutation construction and symbolic verification
f4_model.py            F4 cells from E6 roots
verify_f4.py           F4 exact polynomial certificates
verify_h4_relations.py H4 geometry, normalizer, actions and group-ring checks
verify_h4_spectrum.py  H4 exact spectra and determinant consequences
expected.json         concise expected summaries, fixed before the rerun
PROVENANCE.txt        source digests, adaptations and license status
README.txt            these instructions

No stored matrices, character tables, machine-specific paths, execution logs,
historical manuscripts, or generated caches are needed by this subpackage.
