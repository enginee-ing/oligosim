# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`oligosim` is a forward process model for solid-phase oligonucleotide synthesis.
Given a sequence, its modification pattern, and process conditions, it predicts
the **distribution of species** in the crude product (full-length, truncations,
deletions attributed to the cycle that caused them) rather than a single yield
scalar. Pure Python, standard library only for the core package; `numpy`/`matplotlib`/
`jupyter` are notebook-only extras.

Status: v0.1. The synthesis engine is implemented and tested, but **no kinetic
parameter in this repository is calibrated** — see "Parameter status" below
before changing or adding any default value.

## Commands

```bash
pip install -e ".[dev]"        # install package + pytest
pytest                          # run the full test suite
pytest tests/test_synthesis.py::test_name   # run a single test
pip install -e ".[notebooks]"  # matplotlib/numpy/jupyter/nbformat, for notebook work
python build_notebook.py       # regenerate notebooks/01_sensitivity_analysis.ipynb from source
```

There is no separate lint/format config in this repo — match existing style.

## Architecture

Three modules, each with a single responsibility, wired together by `simulate()`:

- **`oligosim/chemistry.py`** — `Residue` (base + 2' sugar + 3' linkage) and
  `Oligo` (an ordered chain of residues). Masses are computed from elemental
  composition, not tabulated, so they're independently checkable.
- **`oligosim/conditions.py`** — `ProcessConditions`: per-cycle coupling/capping/
  sulfurization efficiencies, with optional per-position overrides. All access
  to coupling efficiency goes through `coupling_at()` so that a future `amidite`
  module (v0.5) can derive it from raw-material quality attributes without
  touching the engine.
- **`oligosim/synthesis.py`** — `simulate()`: propagates the support-bound
  population cycle by cycle and resolves it into `Species` by failure pattern.
  Returns a `SynthesisResult`.

### Key design decisions to preserve

**Synthesis-order indexing.** `Oligo` stores residues 3'→5' internally
(position 1 = the support-loaded 3'-terminal residue), because that's the order
cycles actually run in. `Oligo.sequence_5to3` / `Oligo.from_string` give the
human-facing 5'→3' view. Don't conflate the two orderings when touching
`chemistry.py` or `synthesis.py`.

**State space is species keyed by skipped positions.** `simulate()` tracks
`frozenset[int]` (deleted positions) → mole fraction, exactly up to
`max_deletions` (default 3). Anything beyond that is summed into
`unresolved_fraction` rather than dropped, so `mass_balance` (species fractions
+ unresolved) always closes to 1.0. Any change to the population bookkeeping
must preserve this invariant.

**The capping/deletion mechanism is the point of the model, and is
tested as an analytic invariant** (`tests/test_synthesis.py`,
`test_full_length_equals_product_of_coupling_regardless_of_capping` and
`test_capping_shifts_failures_from_deletions_to_truncations`):
- Full-length fraction == `prod(c_i)` over all cycles, **independent of capping
  efficiency**. A chain is full-length only if it coupled on every cycle;
  capping only acts on chains that already failed.
- What capping efficiency controls is whether a failure becomes a **truncation**
  (capped — short, separates easily) or a **deletion sequence** (escapes
  capping, couples later, differs from product by one internal residue —
  co-elutes closely, hard to separate). This is why capping efficiency matters
  more than its magnitude suggests, even though it doesn't move yield.

If you change the engine, run the full test suite — several tests exist
specifically to catch a regression in one of these invariants, not just to
check example values.

**Nucleoside vs. residue mass.** A "residue" mass includes its 3'-phosphate;
the 3'-terminal residue (position 1, after cleavage from the support) has a
free 3'-OH and must use `nucleoside_mass` instead (residue mass − `HPO3`). An
n-mer carries n−1 phosphates. This is pinned by
`test_oligo_mass_is_residues_plus_water_minus_terminal_phosphate` and the
nusinersen end-to-end mass check — do not "fix" the phosphate count without
re-deriving against those tests.

**Sulfurization is an independent marginal, not joint with the deletion
state** (`SynthesisResult.expected_po_mismatches`,
`fully_sulfurized_fraction`, `po_mismatch_distribution`) — this is a known,
documented approximation in v0.1, not an oversight. Folding it into the joint
population is the first v0.2 item; don't silently "improve" it without
updating the README's "Known gap" section and roadmap to match.

## Parameter status (read before touching defaults)

None of the kinetic parameters (`coupling_efficiency`, `capping_efficiency`,
`sulfurization_efficiency` defaults, `_SUGAR_RELATIVE_REACTIVITY`) are
calibrated values — they are plausible round numbers or, for sugar reactivity,
ordinal placeholders encoding a qualitative ordering only (disabled by default
via `apply_sugar_reactivity=False`). See "Data provenance" below before
adding or changing any parameter.

## Data provenance

Every kinetic parameter added to this repository must carry a public
citation. Never use vendor certificates of analysis, batch records,
customer process parameters, or any data obtained under a commercial
relationship — regardless of how routine the document is. A model whose
provenance cannot be stated publicly is one that cannot be shown to
anyone. When adding a parameter, put the citation in a comment on the
same line.

## Notebooks

`notebooks/01_sensitivity_analysis.ipynb` is generated by `build_notebook.py`
— edit the script, not the notebook file directly, then regenerate with
`python build_notebook.py`.
