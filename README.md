# Verification of semisimplicity criteria for Brauer algebras

This is the finite verification package accompanying *Semisimplicity of
Brauer algebras associated with finite Coxeter types*, by Shoumin Liu and
Yongzhi Luan. It contains the computations used in the exceptional-type
proofs. The structural reductions and the classical and dihedral arguments
are in the paper.

## Run

From an extracted copy or a clone, use SageMath 10.1:

```sh
sage -python verify_all.py
```

The command checks supplied-file hashes, reconstructs every finite model,
performs the exact checks, compares the independently computed results with
the expected results, and exits with a nonzero status on a failed check.
Do not use Python's `-O` option or set `PYTHONOPTIMIZE`.

Tested environment: SageMath 10.1, Python 3.11.1, NumPy 1.23.5, SymPy 1.12,
Debian under WSL. These libraries are included in the tested Sage installation.
No network access, Julia, GAP, or separately installed python-flint is needed.
The complete run takes approximately one minute on the test machine;
performance depends on hardware. The command also works from another current
directory if given the absolute path of `verify_all.py`.

Expected terminal result: checksum PASS, one PASS for each of E6, E7, E8,
G2, H3, F4, H4, and a final PASS. To print the computed certificates, use
`sage -python verify_all.py --details`. These are actual computed results,
not copies of the expected JSON. Runtime and environment are printed separately.

The default command leaves the supplied files unchanged. Intermediate outputs
are kept in temporary directories and removed after the run. The expected
results are comparison targets, never mathematical inputs to the constructors.
For selective runs use `--only ade` or `--only folded`; the component READMEs
give commands for individual types and optional retained reports.

## Contents and proof coverage

| Files | Purpose |
| --- | --- |
| `verify_all.py`, `SHA256SUMS.json` | One entry point; integrity of the distributed files |
| `ade/code/*.py`, `ade/run.py` | E6 exact ranks and jet, E7 fusion/branching and residual ranks, E8 fusion/branching and residual projector checks |
| `ade/code/weyle6.jl`, `weyle7.jl`, `weyle8.jl` | Fixed CHEVIE character-table inputs, evaluated at q=1 |
| `ade/expected/results.json` | Expected ADE spectra, determinant factors, ranks, branching digests, and residual certificates |
| `folded/*.py`, `folded/expected.json` | G2/H3 symbolic identities; F4 exact determinant evaluation; H4 root geometry and exact spectra |
| Component READMEs and provenance files | Precise proof-to-program map, original source hashes, adaptations, and arithmetic conventions |
| `ade/code/CHEVIE_LICENSE` | Required third-party MIT notice |
| `ade/TESTED.json`, `VALIDATION.txt` | Recorded execution environment and verification scope |

The package deliberately contains no manuscript, old research archive, chat
history, generated matrix dumps, caches, or optional low-rank experiments.
All necessary geometry and matrices are reconstructed. The three character
tables are external inputs pinned by file hash; orthogonality and branching
checks do not constitute an independent derivation of those tables.

A successful execution certifies the stated finite calculations. It does not
replace the paper's identification of the matrices with cellular pairings,
normal-form basis theorems, or reductions to the listed residual cases.
Full modular rank at one prime is used only for the E7 integer matrices,
not to infer an entire complex polynomial root set.

## Provenance, permissions, and citation

CHEVIE source data are redistributed with their existing MIT notice. The
research scripts supplied for this project did not carry an explicit license;
this package does not assign a license on behalf of the authors. The authors
can add their chosen research-code license before public release, retaining
the CHEVIE notice. Provenance files identify the adaptations.

When publishing the package, create a fixed release and cite that release or
its commit in the paper. The manuscript has an empty `\verificationurl` macro
for that URL. No repository URL or DOI has been invented here. A DOI archive
of a fixed release is optional; it can be created through Zenodo following
[GitHub's documentation](https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content).
The READMEs already identify the paper and authors; no additional citation
metadata is needed to run the verification.
