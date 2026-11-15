#!/usr/bin/env python3
"""Reconstruct and verify the folded finite calculations in the paper.

Run: sage -python verify_all.py [--only G2 H3 F4 H4] [--details]
All arithmetic assertions run before comparison with the expected summary.
No result files, caches or expected snapshots are written by this entrypoint.
"""
import argparse
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
if not __debug__:
    raise RuntimeError("Run without -O: assertions are the certificate checks.")


def summary(name, result):
    if name == "G2":
        return result
    if name == "H3":
        return {
            "status": "PASS",
            "dimensions": [v["dimension"] for v in result["modules"]],
            "defining_relation_identities": result["total_defining_relation_checks"],
            "contravariance_identities": result["contravariance_identities"],
            "determinants": [v["determinant"] for v in result["modules"]],
            "dimension_square_sum": result["dimension_square_sum"],
            "bad_t": result["bad_t"],
        }
    if name == "F4":
        return {key: value for key, value in result.items() if key != "determinants"} | {
            "cell_orders": {v["label"]: v["order"] for v in result["determinants"]}
        }
    relations, spectra = result
    return {"relations": relations, "spectra": spectra["spectra"],
            "reflection_class_sum_offdiagonal_matches": spectra["reflection_class_sum_offdiagonal_matches"],
            "constant_terms": spectra["constant_terms"],
            "maximal_determinant": spectra["maximal_determinant"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", nargs="+", choices=["G2", "H3", "F4", "H4"])
    parser.add_argument("--details", action="store_true", help="print full computed certificates")
    args = parser.parse_args()
    expected = json.loads((Path(__file__).parent / "expected.json").read_text(encoding="utf-8"))
    selected = args.only or ["G2", "H3", "F4", "H4"]
    outputs = {}
    for name in selected:
        print("Reconstructing and checking " + name + " ...", file=sys.stderr, flush=True)
        if name == "G2":
            from verify_g2 import verify
            result = verify()
        elif name == "H3":
            from verify_h3 import verify
            result = verify()
        elif name == "F4":
            from verify_f4 import verify
            result = verify()
        else:
            from verify_h4_relations import geom, verify as verify_relations
            from verify_h4_spectrum import verify as verify_spectrum
            geometry = geom()
            result = (verify_relations(geometry), verify_spectrum(geometry))
        # Normalize integer JSON keys, e.g. spectral roots, for comparison.
        actual = json.loads(json.dumps(summary(name, result)))
        assert actual == expected[name], name + " differs from the expected summary"
        outputs[name] = result if args.details else actual
        print(name + ": PASS", file=sys.stderr, flush=True)
    print(json.dumps(outputs, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
