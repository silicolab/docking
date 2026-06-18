# Golden ground-truth values

These files are reference outputs from **official AutoDock Vina v1.2.7** (the official
precompiled Windows binary), used to validate the `docking` crate. The Rust engine must
reproduce these numbers within the documented tolerances.

**Do not edit these by hand.** Regenerate them only intentionally, via:

```bash
bash tests/golden/regenerate.sh        # set $VINA to the official Vina binary
```

and review the resulting diff carefully — these are the reference the tests check against.

## Provenance

- Engine: AutoDock Vina v1.2.7, run single-threaded (`--cpu 1`), verified bit-for-bit
  reproducible for a fixed `--seed`.
- Inputs (`receptor.pdbqt`, `ligand.pdbqt`, `config.txt`) are the prepared PDBQT and box
  taken verbatim from Vina's own `example/basic_docking/solution` (system `1iep`:
  Abl kinase + STI-571 / imatinib).

## Files per system

| File                   | Produced by                          | Validates    |
|------------------------|--------------------------------------|--------------|
| `score_only.out.txt`   | `--score_only`                       | scoring      |
| `local_only.out.txt`   | `--local_only`                       | BFGS energy  |
| `local_only.pdbqt`     | `--local_only --out`                 | BFGS pose    |
| `dock_seed42.out.txt`  | `--seed 42 --exhaustiveness 8`       | search       |
| `dock_seed42.pdbqt`    | `--seed 42 --exhaustiveness 8 --out` | search poses |

The `tests/fixtures/pdbqt/` files are a separate, feature-agnostic set used only for
PDBQT parser round-trip tests; they are not scoring ground truth.
