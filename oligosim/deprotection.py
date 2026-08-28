"""Base modification during cleavage/deprotection.

Ammonia deprotection removes the base- and phosphate-protecting groups and
cleaves the oligo from the solid support. It is also base-catalysed,
heat-accelerated hydrolysis that DEAMINATES cytosine residues: cytosine (C)
loses an amine and gains a carbonyl oxygen to become uracil (U);
5-methylcytosine (mC) undergoes the identical chemistry to become thymine
(T) -- not an analogue, the same base every unmodified synthesis position
already carries. This happens during the deprotection incubation, not
during a synthesis cycle, so it is modelled as a single post-synthesis
stage rather than folded into `synthesis.simulate()`'s per-cycle state.

Why this class of impurity matters more than its abundance suggests
---------------------------------------------------------------------
A deamination event is a net loss of NH and gain of O: +0.98402 Da
monoisotopic (`chemistry.DEAMINATION_DELTA`) on a molecule of order 7 kDa --
about a 0.014% mass shift, well inside the width of the natural isotope
envelope. The deaminated species is still full length, still carries the
same net charge (deamination doesn't touch the phosphate backbone), and
elutes with product. Unlike a truncation or an n-1 deletion, there is no
obvious separation handle: it isn't shorter, it isn't missing a residue,
and it isn't a resolvable mass. It ships.

Architecture
------------
Deprotection is a single post-synthesis stage: it does not add a dimension
to `simulate()`'s per-cycle state (deamination timing has nothing to do
with which cycle a residue was added in, and doing that would multiply the
state space for no resolvable benefit, same reasoning as the deletion/
mismatch state design in `synthesis.py`). Given a finished `SynthesisResult`,
`deprotect()` expands each of its species into sub-species by deamination
count, binomial over *that species'* count of C and 5-methyl-C residues --
which depends on which positions were deleted (a deleted C isn't there to
deaminate) and, for a truncated species, on how much of the sequence was
even synthesized. So the eligible-residue count is computed per species,
not once for the parent oligo.

Parameter status
-----------------
`deamination_probability_per_residue` is used directly, not derived from a
time/temperature Arrhenius model. No calibrated, publicly citable rate
constant exists in the literature for this reaction in the deprotection-
cocktail context, and inventing one would violate this repository's
provenance rule (see CLAUDE.md "Data provenance"): every kinetic parameter
must carry a public citation, and a plausible-sounding rate constant is not
one. A t/T model is future work, pending a citable source. The default of
0.0 makes this stage a no-op, matching behaviour before this module existed.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import comb

from .chemistry import Base, DEAMINATION_DELTA, Oligo
from .synthesis import Species, SynthesisResult


def _check_fraction(name: str, value: float) -> None:
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must lie in [0, 1], got {value}")


@dataclass
class DeprotectionConditions:
    """Deprotection-stage parameters.

    Parameters
    ----------
    deamination_probability_per_residue
        Probability that a given cytosine (C) or 5-methylcytosine (mC)
        residue deaminates during ammonia deprotection (C -> U; mC -> T).
        A direct probability, not derived from time/temperature -- see the
        module docstring's "Parameter status" for why. Defaults to 0.0
        (no-op).
    """

    deamination_probability_per_residue: float = 0.0

    def __post_init__(self) -> None:
        _check_fraction(
            "deamination_probability_per_residue",
            self.deamination_probability_per_residue,
        )


@dataclass(frozen=True)
class DeprotectedSpecies:
    """One species after cleavage/deprotection: a synthesis `Species`
    expanded by deamination count.

    Attributes
    ----------
    parent
        The pre-deprotection `Species` this was expanded from. Kept as a
        reference rather than flattened onto this class, so
        `deleted`/`truncated`/`mismatches`/`length` have one source of
        truth.
    deaminations
        Count of C->U / 5-methyl-C->T events. Not attributed to a specific
        residue -- same rationale as `Species.mismatches`: count is what
        determines mass, and tracking position would multiply the state
        space for no resolvable benefit.
    fraction
        Mole fraction of the total support-bound population (same basis as
        `SynthesisResult.species` fractions).
    mass
        Monoisotopic mass, accounting for `deaminations`.
    """

    parent: Species
    deaminations: int
    fraction: float
    mass: float

    @property
    def label(self) -> str:
        """Short identifier, extending `parent.label` with a `+da{n}`
        suffix when deaminated (mirrors `Species.label`'s `+mm{n}`)."""
        return f"{self.parent.label}+da{self.deaminations}" if self.deaminations else self.parent.label


@dataclass
class DeprotectionResult:
    """Outcome of applying deprotection-stage chemistry to a
    `SynthesisResult`."""

    synthesis_result: SynthesisResult
    conditions: DeprotectionConditions
    species: list[DeprotectedSpecies]
    unresolved_fraction: float
    max_deaminations: int

    @property
    def mass_balance(self) -> float:
        """Total accounted mole fraction. Should be 1.0 within float error."""
        return sum(s.fraction for s in self.species) + self.unresolved_fraction

    @property
    def correct_product_fraction(self) -> float:
        """Mole fraction of the fully correct, fully deprotected molecule:
        `parent.deleted` empty, not truncated, zero PS/PO mismatches, *and*
        zero deaminations. The narrowest possible "product" definition --
        narrower than `SynthesisResult.correct_product_fraction`, which
        this stage can only shrink further, never grow.
        """
        return sum(
            s.fraction
            for s in self.species
            if not s.parent.deleted
            and not s.parent.truncated
            and s.parent.mismatches == 0
            and s.deaminations == 0
        )


def _deamination_eligible_count(oligo: Oligo, species: Species) -> int:
    """Count of C / 5-methyl-C residues actually present in `species`.

    Computed per species, not once for the parent oligo: a deleted position
    that happened to carry a C or mC isn't present to deaminate, and a
    truncated species only carries residues up to its truncation point.

    `species.length + len(species.deleted)` recovers the upper bound of
    "positions that would exist absent deletion/truncation" for both cases
    with one formula: for a truncated species this is (cycle - 1), since
    `length = (cycle - 1) - len(deleted)`; for a full-length species it
    reduces to `oligo.n`, since `length = oligo.n - len(deleted)`.
    """
    upper = species.length + len(species.deleted)
    return sum(
        1
        for position in range(1, upper + 1)
        if position not in species.deleted
        and oligo.residue_at(position).base in (Base.C, Base.mC)
    )


def deprotect(
    synthesis_result: SynthesisResult,
    conditions: DeprotectionConditions | None = None,
    max_deaminations: int = 2,
) -> DeprotectionResult:
    """Expand a `SynthesisResult` by deamination count.

    Parameters
    ----------
    max_deaminations
        Deamination states are tracked exactly up to this many events per
        species. Beyond it, probability mass is lumped into
        `unresolved_fraction` so mass balance is preserved and the
        truncation is visible rather than silent -- the same mechanism
        `simulate()` uses for `max_deletions` / `max_mismatches`.

    Returns
    -------
    DeprotectionResult
    """
    if conditions is None:
        conditions = DeprotectionConditions()
    if max_deaminations < 0:
        raise ValueError("max_deaminations must be >= 0")

    oligo = synthesis_result.oligo
    p = conditions.deamination_probability_per_residue

    expanded: list[DeprotectedSpecies] = []
    unresolved = synthesis_result.unresolved_fraction

    for parent in synthesis_result.species:
        n_eligible = _deamination_eligible_count(oligo, parent)

        # Nothing to deaminate, or deamination is switched off: reproduce
        # the parent exactly rather than running a trivial binomial (also
        # guarantees bit-identical output at p=0.0, not just approx-equal).
        if n_eligible == 0 or p == 0.0:
            expanded.append(
                DeprotectedSpecies(
                    parent=parent,
                    deaminations=0,
                    fraction=parent.fraction,
                    mass=parent.mass,
                )
            )
            continue

        for j in range(0, n_eligible + 1):
            prob = comb(n_eligible, j) * (p**j) * ((1.0 - p) ** (n_eligible - j))
            amount = parent.fraction * prob
            if amount == 0.0:
                continue
            if j > max_deaminations:
                unresolved += amount
            else:
                expanded.append(
                    DeprotectedSpecies(
                        parent=parent,
                        deaminations=j,
                        fraction=amount,
                        mass=parent.mass + j * DEAMINATION_DELTA,
                    )
                )

    return DeprotectionResult(
        synthesis_result=synthesis_result,
        conditions=conditions,
        species=expanded,
        unresolved_fraction=unresolved,
        max_deaminations=max_deaminations,
    )


__all__ = [
    "DeprotectionConditions",
    "DeprotectedSpecies",
    "DeprotectionResult",
    "deprotect",
]
