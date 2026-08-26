"""Forward simulation of solid-phase oligonucleotide synthesis.

What this does that a spreadsheet does not
------------------------------------------
The standard industry estimate is

    full-length fraction = (coupling efficiency) ** (n - 1)

which is a scalar. It cannot tell you *which* impurities you have, and so it
cannot tell you whether they will separate. This module instead propagates the
population through each cycle and resolves it by failure pattern, which is what
determines chromatographic behaviour downstream.

The mechanism, per cycle i (adding residue i)
---------------------------------------------
1. Detritylation exposes the 5'-OH of every active chain.
2. Coupling: a fraction c_i of active chains extends. The rest retain a free
   5'-OH.
3. Capping: a fraction k_i of those failures is acetylated and leaves the
   growing population permanently, becoming a TRUNCATION of length i-1-|D|.
4. The failures that escape capping stay active with a free 5'-OH and couple in
   a *later* cycle. They become DELETION sequences: full-length-minus-one but
   missing an internal residue.

That distinction is the reason capping exists. A truncation is short and elutes
far from product. A deletion sequence differs from product by one internal
residue and co-elutes closely, which is why (1 - capping efficiency) drives
crude purity harder than its magnitude suggests.

Sulfurization shortfall -- the chance a phosphorothioate (PS) linkage ends up
as a plain phosphodiester (PO) -- is tracked jointly with the deletion state:
each species carries a mismatch count alongside its deleted-position set, so
mass and identity reflect the actual backbone composition. Only the *count*
is tracked, not which linkages mismatched, since count is what determines
charge and hydrophobicity (and therefore chromatographic behaviour), while
tracking positions would multiply the state space by 2^(n_ps) for no
resolvable benefit.

One consequence: `full_length_fraction == prod(c_i)` exactly (the invariant
tested in test_synthesis.py) only for PO oligos, or for PS oligos where
`max_mismatches` never binds. A chain that stays on the zero-deletion branch
but racks up more than `max_mismatches` PS/PO mismatches is shunted into
`unresolved_fraction`, which strips its "full-length" identity along with its
mismatch state -- there's no way to know it was deletion-free once it's in
that bucket. In practice this is a small effect (see `correct_product_fraction`
vs. `naive_full_length_fraction` on a real PS oligo), but it means
`full_length_fraction` for a PS oligo is a lower bound on `prod(c_i)`, not an
exact equality.

Approximations in v0.1 (each is a named v0.2+ item)
---------------------------------------------------
* Depurination is not modelled. It is acid-catalysed, dA-dominated, and
  accumulates with cycle count.
* n+1 insertions from premature detritylation are not modelled.
* Cleavage and deprotection losses are not modelled.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from math import comb

from .chemistry import Linkage, Oligo, PS_DELTA
from .conditions import ProcessConditions


@dataclass(frozen=True)
class Species:
    """One resolved chemical species in the crude mixture.

    Attributes
    ----------
    deleted
        1-based synthesis positions missing from this species.
    mismatches
        Count of PS linkages that formed as PO instead (sulfurization
        shortfall). Not attributed to a position -- see module docstring.
    length
        Residue count.
    fraction
        Mole fraction of the total support-bound population.
    truncated
        True if the chain was capped and stopped growing; False if it ran to
        the end of the synthesis.
    mass
        Monoisotopic mass, accounting for `mismatches`.
    """

    deleted: frozenset[int]
    mismatches: int
    length: int
    fraction: float
    truncated: bool
    mass: float

    @property
    def label(self) -> str:
        """Short identifier, e.g. 'FLP', 'trunc:12mer', 'del:n-1@7'."""
        if not self.deleted and not self.truncated:
            base = "FLP"
        elif self.truncated:
            base = f"trunc:{self.length}mer"
            if self.deleted:
                base += "+d" + ",".join(str(p) for p in sorted(self.deleted))
        else:
            deficit = len(self.deleted)
            positions = ",".join(str(p) for p in sorted(self.deleted))
            base = f"del:n-{deficit}@{positions}"
        return f"{base}+mm{self.mismatches}" if self.mismatches else base


@dataclass
class SynthesisResult:
    """Outcome of one simulated synthesis."""

    oligo: Oligo
    conditions: ProcessConditions
    species: list[Species]
    unresolved_fraction: float
    max_deletions: int
    max_mismatches: int

    # -- headline numbers ---------------------------------------------------
    @property
    def correct_product_fraction(self) -> float:
        """Mole fraction of the fully correct molecule: full length, no
        deletions, and every PS linkage correctly sulfurized (mismatches==0).
        This is the narrowest "product" definition -- the actual target
        molecule, not just the right length with an arbitrary backbone.
        """
        return sum(
            s.fraction
            for s in self.species
            if not s.deleted and not s.truncated and s.mismatches == 0
        )

    @property
    def full_length_fraction(self) -> float:
        """Mole fraction of the right length (no deletions), any backbone.

        Includes PS/PO mismatches -- a chain with the correct sequence but a
        sulfurization defect still counts here, unlike
        `correct_product_fraction`. Equal to `prod(c_i)` for PO oligos, or
        for PS oligos where `max_mismatches` doesn't bind (see module
        docstring); otherwise a lower bound on it.
        """
        return sum(s.fraction for s in self.species if not s.deleted and not s.truncated)

    @property
    def deletion_fraction(self) -> float:
        """Full-synthesis chains missing one or more internal residues."""
        return sum(s.fraction for s in self.species if s.deleted and not s.truncated)

    @property
    def truncation_fraction(self) -> float:
        return sum(s.fraction for s in self.species if s.truncated)

    @property
    def naive_full_length_fraction(self) -> float:
        """The spreadsheet estimate, for comparison: prod(c_i) over cycles."""
        value = 1.0
        for pos in range(2, self.oligo.n + 1):
            value *= self.conditions.coupling_at(pos, self.oligo)
        return value

    @property
    def mass_balance(self) -> float:
        """Total accounted mole fraction. Should be 1.0 within float error."""
        return sum(s.fraction for s in self.species) + self.unresolved_fraction

    # -- sulfurization (derived from the joint state) ------------------------
    @property
    def fully_sulfurized_fraction(self) -> float:
        """Among *full-length* species (no deletions, not truncated), the
        fraction with zero PO-for-PS mismatches.

        Deliberately restricted to full-length species: a truncation has
        fewer PS linkages by construction, so it is trivially more likely to
        be mismatch-free, and folding it in would inflate this number
        relative to what "fully sulfurized" means for the actual product.

        Approaches `sulfurization_efficiency ** n_ps_linkages(oligo)` as
        `max_mismatches` becomes non-binding, since conditional on a chain
        being full-length, its mismatch count is an independent
        Binomial(n_ps, 1 - sulfurization_efficiency) regardless of overall
        coupling yield.
        """
        denom = self.full_length_fraction
        if denom == 0.0:
            return 0.0
        return self.correct_product_fraction / denom

    def po_mismatch_distribution(self, max_mismatches: int | None = None) -> dict[int, float]:
        """Mole fraction grouped by PO-for-PS mismatch count.

        Derived from the tracked species, not an independent closed form, so
        it can only resolve counts up to `self.max_mismatches` (the cap
        `simulate()` was run with). Population that overflowed that cap --
        or the deletion cap -- lives in `unresolved_fraction` and isn't
        attributable to a specific count, so the returned distribution can
        sum to slightly less than 1.0 when overflow occurred. `mass_balance`,
        not this method, is the quantity guaranteed to equal 1.0.
        """
        cap = self.max_mismatches if max_mismatches is None else min(max_mismatches, self.max_mismatches)
        out: dict[int, float] = defaultdict(float)
        for s in self.species:
            if s.mismatches <= cap:
                out[s.mismatches] += s.fraction
        return dict(sorted(out.items()))

    # -- views --------------------------------------------------------------
    def top_species(self, k: int = 10) -> list[Species]:
        return sorted(self.species, key=lambda s: -s.fraction)[:k]

    def deletions_by_position(self) -> dict[int, float]:
        """Total single-deletion fraction attributable to each position.

        Answers 'which cycle is costing me purity', which is the question a
        process chemist actually asks when a batch misses spec.
        """
        out: dict[int, float] = defaultdict(float)
        for s in self.species:
            if s.truncated or len(s.deleted) != 1:
                continue
            out[next(iter(s.deleted))] += s.fraction
        return dict(sorted(out.items()))

    def summary(self) -> str:
        # Labels vary in length ("Correct length (any backbone)" is the
        # longest), so pad to a fixed column instead of hand-aligning.
        def row(label: str, value: str) -> str:
            return f"{label:<30}: {value}"

        lines = [
            row("Oligo", f"{self.oligo.sequence_5to3}  (n={self.oligo.n})"),
            row("Couplings", str(self.oligo.n_couplings)),
            "",
            row("Correct product (FLP)", f"{self.correct_product_fraction:.4%}"),
            row("Correct length (any backbone)", f"{self.full_length_fraction:.4%}"),
            row("  naive c^(n-1)", f"{self.naive_full_length_fraction:.4%}"),
            row("Deletions", f"{self.deletion_fraction:.4%}"),
            row("Truncations", f"{self.truncation_fraction:.4%}"),
            row(
                "Unresolved",
                f"{self.unresolved_fraction:.4%}  (>{self.max_deletions} deletions)",
            ),
            row("Mass balance", f"{self.mass_balance:.10f}"),
        ]
        n_ps = self.conditions.n_ps_linkages(self.oligo)
        if n_ps:
            dist = self.po_mismatch_distribution()
            expected_mismatches = sum(k * v for k, v in dist.items())
            lines += [
                "",
                row("PS linkages", str(n_ps)),
                row("Fully sulfurized (of full-length)", f"{self.fully_sulfurized_fraction:.4%}"),
                row(
                    "E[PO mismatches]",
                    f"{expected_mismatches:.4f}  (resolved to {self.max_mismatches})",
                ),
            ]
        lines += ["", "Top species:"]
        for s in self.top_species(8):
            lines.append(f"  {s.label:<22} {s.fraction:>9.4%}  {s.mass:>10.3f} Da")
        return "\n".join(lines)


def simulate(
    oligo: Oligo,
    conditions: ProcessConditions | None = None,
    max_deletions: int = 3,
    max_mismatches: int = 3,
) -> SynthesisResult:
    """Propagate the support-bound population through the synthesis.

    Parameters
    ----------
    max_deletions
        Deletion states are tracked exactly up to this many missing residues.
        Species beyond it are lumped into `unresolved_fraction` so that mass
        balance is preserved and the truncation is visible rather than silent.
        Three is ample in practice: at 99% coupling on a 20-mer, four-deletion
        species are well below any analytical limit of quantitation.
    max_mismatches
        PO-for-PS sulfurization shortfall is tracked exactly up to this many
        mismatched linkages, by the same lumping-into-`unresolved_fraction`
        mechanism as `max_deletions`. Only the count is tracked, not which
        linkages mismatched -- see module docstring.

    Returns
    -------
    SynthesisResult
    """
    if conditions is None:
        conditions = ProcessConditions()
    if max_deletions < 0:
        raise ValueError("max_deletions must be >= 0")
    if max_mismatches < 0:
        raise ValueError("max_mismatches must be >= 0")

    # Active population, keyed by (positions already skipped, PO-for-PS
    # mismatch count so far). Position 1 is preloaded on the support, so
    # everything starts intact and fully sulfurized.
    active: dict[tuple[frozenset[int], int], float] = {(frozenset(), 0): 1.0}
    truncations: list[Species] = []
    unresolved = 0.0

    for position in range(2, oligo.n + 1):
        c = conditions.coupling_at(position, oligo)
        k = conditions.capping_efficiency
        s = conditions.sulfurization_efficiency
        is_ps = oligo.residue_at(position).linkage is Linkage.PS
        nxt: dict[tuple[frozenset[int], int], float] = defaultdict(float)

        for (deleted, mm), amount in active.items():
            if amount == 0.0:
                continue

            # 1. Coupled: extends, deletion set unchanged. If this linkage is
            # a PS position, split further by sulfurization shortfall.
            coupled = amount * c
            if coupled > 0.0:
                if is_ps:
                    nxt[(deleted, mm)] += coupled * s
                    mismatched = coupled * (1.0 - s)
                    if mismatched > 0.0:
                        new_mm = mm + 1
                        if new_mm > max_mismatches:
                            unresolved += mismatched
                        else:
                            nxt[(deleted, new_mm)] += mismatched
                else:
                    nxt[(deleted, mm)] += coupled

            failed = amount * (1.0 - c)
            if failed == 0.0:
                continue

            # 2. Capped: leaves the population as a truncation. No new
            # linkage formed this cycle, so mismatch count carries over.
            capped = failed * k
            if capped > 0.0:
                length = (position - 1) - len(deleted)
                truncated_positions = deleted | frozenset(range(position, oligo.n + 1))
                truncations.append(
                    Species(
                        deleted=deleted,
                        mismatches=mm,
                        length=length,
                        fraction=capped,
                        truncated=True,
                        # A truncation is missing everything from `position` up.
                        mass=oligo.mass(truncated_positions) - mm * PS_DELTA,
                    )
                )

            # 3. Escaped capping: stays active, gains a deletion here.
            escaped = failed * (1.0 - k)
            if escaped > 0.0:
                new_deleted = deleted | {position}
                if len(new_deleted) > max_deletions:
                    unresolved += escaped
                else:
                    nxt[(frozenset(new_deleted), mm)] += escaped

        active = dict(nxt)

    full_synthesis = [
        Species(
            deleted=deleted,
            mismatches=mm,
            length=oligo.n - len(deleted),
            fraction=amount,
            truncated=False,
            mass=oligo.mass(deleted) - mm * PS_DELTA,
        )
        for (deleted, mm), amount in active.items()
        if amount > 0.0
    ]

    return SynthesisResult(
        oligo=oligo,
        conditions=conditions,
        species=full_synthesis + truncations,
        unresolved_fraction=unresolved,
        max_deletions=max_deletions,
        max_mismatches=max_mismatches,
    )


def state_space_size(n: int, max_deletions: int) -> int:
    """Number of exactly-tracked deletion states for an n-mer.

    Useful for sanity-checking cost before running a long sequence.
    """
    return sum(comb(n - 1, k) for k in range(0, max_deletions + 1))


__all__ = [
    "Species",
    "SynthesisResult",
    "simulate",
    "state_space_size",
    "combinations",
]
