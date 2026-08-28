# oligosim

A forward process model for solid-phase oligonucleotide synthesis. Given a
sequence, its modification pattern, and process conditions, it predicts the
**distribution of species** in the crude product — not just a yield number.

> **Status: v0.2.dev0.** v0.1's engine is implemented and tested; joint
> sulfurization, incomplete detritylation, and deamination — the first three
> v0.2 items — have landed ahead of the rest of v0.2. Kinetic parameters are
> placeholders, not calibrated values. See
> [Parameter status](#parameter-status) before using any output quantitatively.

## The problem

The standard estimate of oligonucleotide synthesis yield is

```
full-length fraction = (coupling efficiency) ** (n - 1)
```

This is a scalar, and it answers only one question. It cannot tell you *which*
impurities you have, so it cannot tell you whether they will separate — which
is the question that actually determines whether a batch meets spec.

Two crude mixtures with identical 87% full-length content can behave completely
differently downstream. If the missing 13% is short truncations, purification is
routine. If it is n−1 deletion sequences, it is hard, because a chain missing
one internal residue co-elutes closely with product.

`oligosim` propagates the support-bound population cycle by cycle and resolves
it by failure pattern, so truncations and deletions are distinguished and
attributed to the cycle that produced them.

## The mechanism

Per cycle *i*:

1. **Detritylation** deblocks a fraction `d_i` of active chains, exposing a
   free 5'-OH; the rest stay DMT-protected. A chain that fails to deblock
   can't couple this cycle, and — critically — can't be capped either, since
   capping acetylates free hydroxyls. It becomes a **deletion sequence**
   unconditionally: a second route to that outcome, and one capping
   efficiency has no purchase on at all.
2. **Coupling** extends a fraction `c_i` of the chains that did deblock. The
   rest keep a free 5'-OH.
3. **Capping** acetylates a fraction `k_i` of those coupling failures,
   removing them permanently as **truncations**.
4. Coupling failures that **escape capping** stay active and couple in a
   later cycle, becoming **deletion sequences** by the same outcome as an
   ordinary failed coupling — full-length-minus-one, missing an internal
   residue.

Steps 1 and 4 are why capping only gets you so far. Capping does not improve
yield; the correct-product fraction is `prod(d_i * c_i) * s^n_ps` (folding in
sulfurization) regardless of capping efficiency. What capping does is convert
coupling failures that *did* deblock into truncations instead of deletions,
trading a hard separation problem for an easy one — it has no effect on chains
that never deblocked in the first place. The model reproduces this exactly,
and it is enforced as a test invariant.

**Deprotection is a separate, post-synthesis stage.** Everything above
happens inside the cycle loop, tracked per synthesis cycle. Deprotection does
not: it is base-catalysed, heat-accelerated hydrolysis during the ammonia
deprotection step, applied once to the *finished* species distribution
rather than folded into the cycle-by-cycle state. It deaminates cytosine (C)
to uracil and 5-methylcytosine (mC) to thymine — literally thymine, the same
base every unmodified position already carries, not an analogue.

`deprotect()` expands each synthesis species into sub-species by deamination
count, binomial over `deamination_probability_per_residue` (`p`) and that
species' own count of eligible C/mC residues (`n_eligible`). `n_eligible` is
computed **per species, not per oligo**: a deleted position that happened to
carry a C or mC isn't present to deaminate, and a truncated species only
carries residues up to its truncation point, so two species from the same
oligo can have different eligible counts. Folding this into the invariant
above:

```
correct product fraction = prod(d_i * c_i) * s^n_ps * (1 - p)^n_eligible
```

`deamination_probability_per_residue` defaults to 0.0 — a direct probability,
not a time/temperature model, since no calibrated public rate constant
exists for this reaction; a t/T model is future work pending a literature
value.

## Quick start

```python
from oligosim import Oligo, ProcessConditions, Sugar, Linkage, simulate

# Nusinersen: 18-mer, uniform 2'-O-MOE, full phosphorothioate, 5-methyl-C.
oligo = Oligo.from_string(
    "TmCAmCTTTmCATAATGmCTGG", sugar=Sugar.MOE, linkage=Linkage.PS
)

result = simulate(
    oligo,
    ProcessConditions(
        coupling_efficiency=0.992,
        capping_efficiency=0.95,
        sulfurization_efficiency=0.995,
    ),
)
print(result.summary())
```

```
Oligo                         : TCACTTTCATAATGCTGG [2'-MOE, PS]  (n=18)
Couplings                     : 17

Correct product (FLP)         : 80.1108%
Correct length (any backbone) : 87.2364%
  naive c^(n-1)               : 87.2365%
Deletions                     : 0.5999%
Truncations                   : 12.1635%
Unresolved                    : 0.0001%  (>3 deletions)
Mass balance                  : 1.0000000000

PS linkages                   : 17
Fully sulfurized (of full-length): 91.8318%
E[PO mismatches]              : 0.0794  (resolved to 3)

... (top species table omitted)
```

Attribute impurity to the cycle that caused it:

```python
result.deletions_by_position()   # {2: 0.000352, 3: 0.000352, ...}
```

Mixed chemistry is supported per position — gapmers, siRNA strands, chimeras:

```python
sugars = [Sugar.MOE] * 5 + [Sugar.DNA] * 10 + [Sugar.MOE] * 5
gapmer = Oligo.from_string("A" * 20, sugar=sugars, linkage=Linkage.PS)
```

## Design notes

**Modification-aware from the start.** Therapeutic oligos are essentially never
plain DNA, so `Residue` carries base, 2' sugar and 3' linkage. A model built on
ACGT strings has to be rewritten the moment it meets a real drug substance.

**State space.** Species are keyed by the set of skipped positions and the
count of PO-for-PS sulfurization mismatches — not which linkages mismatched,
since count (not position) is what drives charge and hydrophobicity, and
therefore chromatographic behaviour. Exact tracking to `max_deletions` and
`max_mismatches` (both default 3); anything beyond either cap is reported as
`unresolved_fraction` rather than silently dropped, so mass balance always
closes to 1.0. For a 20-mer at `max_deletions=3` that is 1160 deletion states,
further multiplied by up to 4 mismatch states, and runs in milliseconds.

**Masses are computed, not tabulated**, from elemental composition, so they are
independently checkable. The 3'-terminal residue contributes a *nucleoside*
mass — after cleavage from support it has a free 3'-OH — so an n-mer carries
n−1 phosphates. Getting this wrong adds a spurious HPO₃ to every mass, and the
test suite pins it.

**Validation.** The nusinersen monoisotopic mass is checked end-to-end against
C₂₃₄H₃₄₀N₆₁O₁₂₈P₁₇S₁₇, derived from nucleoside composition. As a cross-check,
the same formula gives an average MW of 7127.2 for the free acid and 7500.9 for
the 17-fold sodium salt, matching the cited value for nusinersen sodium.

## Known gap: this model is optimistic

For an 18-mer 2'-OMe phosphorothioate oligonucleotide, published work on
membrane-enabled liquid-phase synthesis reports crude purity around **72%**.
This model predicts **80.1%** correct product (right length, correctly
sulfurized backbone) at 99.2% coupling — down from the 87.2% you'd get by
checking length alone, since a genuinely correct molecule needs both. Folding
PO-for-PS sulfurization shortfall into the joint population state closed
roughly half the ~15 point gap to the published figure.

What remains, ~8 points, is the set of impurity classes still not modelled:

| Missing | Effect |
|---|---|
| Depurination | Acid-catalysed, dA-dominated, accumulates with cycle count. |
| Cyanoethyl adducts | Acrylonitrile released during deprotection is a Michael acceptor; adducts form preferentially on thymine. |
| n+1 insertions | From *premature* detritylation of the incoming amidite — a different failure mode from *incomplete* detritylation, modelled as a second, capping-immune route to deletions. |
| Cleavage losses | Physical yield loss during support cleavage — distinct from cytosine/5-methylcytosine deamination during the same deprotection step, which is now modelled (see "The mechanism"). |

## Notebooks

[`notebooks/01_sensitivity_analysis.ipynb`](notebooks/01_sensitivity_analysis.ipynb)
— what the impurity profile depends on, with the charts rendered inline.

The strongest result is §5, detritylation vs. capping, swept over their full
ranges on the same axes: capping efficiency doesn't move correct product at
all (flat at 77.4%, holding a representative detritylation loss fixed), while
detritylation efficiency falls by about 12.6 points over a range of just
*one* percentage point (99% to 100%, correct product 67.5% → 80.1%). Capping
only relocates impurity between truncations and deletions; a failed
detritylation is an unconditional, capping-immune deletion that destroys
product outright — no amount of capping buys it back.

The rest of the notebook covers coupling efficiency, length, capping alone
(§3: correct product is flat at 80.111% across the *entire* capping range,
not approximately — exactly, since capping never touches a chain that was
already going to be correct product), and per-cycle attribution.

## Parameter status

**None of the kinetic parameters in this repository are calibrated.** Defaults
are plausible round numbers; the sugar reactivity factors are ordinal
placeholders encoding a documented qualitative ordering (bulky 2' substituents
and bridged sugars couple more slowly than DNA) and are disabled by default.

Use the model for **sensitivity analysis** — how does crude purity respond to a
range of coupling efficiency? — rather than point prediction, until the
parameters are replaced with literature values.

Every parameter added to this repository must carry a public citation. No
vendor certificates of analysis, batch records, or other data obtained under a
commercial relationship, regardless of how routine the document is. A model
whose provenance cannot be stated is a model that cannot be shown to anyone.

## Prior art

This is not a novel modelling approach. It is, as far as I can find, the first
open implementation, and the first to connect the stages end to end.

- **Synthesis kinetics.** *Kinetic Modeling of Solid-Phase Oligonucleotide
  Synthesis: Mechanistic Insights and Reaction Dynamics*, Org. Process Res.
  Dev. 2025, 29(9), 2298–2309. A mechanistic kinetic treatment of coupling,
  capping, oxidation and detritylation. Not available as code.
- **Error propagation.** Earlier work relating constant coupling and capping
  efficiencies to the product length distribution.
- **Purification.** Mechanistic ion-exchange models predicting purity and yield
  in collected fractions, including N−1 and N+1 content. Commercial software.
- **Scale-up.** CFD treatment of packed-bed synthesis reactors, Biotechnol.
  Prog. 2014.

The gap this fills: those exist as separate closed silos. Nothing connects raw
material specification → synthesis → purification → cost of goods in one open
model, and none of them treat **amidite quality attributes** as inputs.

## Roadmap

| Version | Scope |
|---|---|
| **v0.1** | Positional failure propagation; modification-aware chemistry; mass assignment. ✅ |
| v0.2 | Joint sulfurization state ✅. Incomplete detritylation ✅. Deamination ✅. Depurination, cyanoethyl adducts, n+1 insertions still open. 🚧 |
| v0.3 | Full predicted mass spectrum with isotope envelopes. |
| v0.4 | Chromatography: retention model, resolution, pool cut points, purity/yield trade-off curve. |
| v0.5 | `amidite` module — derive per-cycle coupling efficiency from water content, ³¹P purity, free acid, related substances, activator and excess. Inverts into spec setting. |
| v0.6 | Cost of goods: amidite consumption, solvent volume, waste, scale. |

v0.5 is the point of the project. Everything before it is the substrate that
makes the question answerable:

> *What amidite specification do I need to hit 95% crude purity on a 20-mer
> gapmer?*

## Install

```bash
git clone <repo>
cd oligosim
pip install -e ".[dev]"
pytest
```

Pure Python, standard library only. `pytest` for tests.

## Licence

MIT.
