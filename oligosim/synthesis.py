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

Approximations in v0.1 (each is a named v0.2+ item)
---------------------------------------------------
* Depurination is not modelled. It is acid-catalysed, dA-dominated, and
  accumulates with cycle count.
* Sulfurization shortfall is computed as an independent marginal over formed
  linkages rather than jointly with the deletion state. Reported separately.
* n+1 insertions from premature detritylation are not modelled.
* Cleavage and deprotection losses are not modelled.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from math import comb

from .chemistry import Oligo
from .conditions import ProcessConditions


@dataclass(frozen=True)
class Species:
    """One resolved chemical species in the crude mixture.

    Attributes
    ----------
    deleted
        1-based synthesis positions missing from this species.
    length
        Residue count.
    fraction
        Mole fraction of the total support-bound population.
    truncated
        True if the chain was capped and stopped growing; False if it ran to
        the end of the synthesis.
    mass
        Monoisotopic mass, ignoring sulfurization shortfall.
    """

    deleted: frozenset[int]
    length: int
    fraction: float
    truncated: bool
    mass: float

    @property
    def label(self) -> str:
        """Short identifier, e.g. 'FLP', 'trunc:12mer', 'del:n-1@7'."""
        if not self.deleted and not self.truncated:
            return "FLP"
        if self.truncated:
            tag = f"trunc:{self.length}mer"
            if self.deleted:
                tag += "+d" + ",".join(str(p) for p in sorted(self.deleted))
            return tag
        deficit = len(self.deleted)
        positions = ",".join(str(p) for p in sorted(self.deleted))
        return f"del:n-{deficit}@{positions}"


@dataclass
class SynthesisResult:
    """Outcome of one simulated synthesis."""

    oligo: Oligo
    conditions: ProcessConditions
    species: list[Species]
    unresolved_fraction: float
    max_deletions: int

    # -- headline numbers ---------------------------------------------------
    @property
    def full_length_fraction(self) -> float:
        """Mole fraction of correct full-length product (FLP)."""
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

    # -- sulfurization (independent marginal; see module docstring) ---------
    @property
    def expected_po_mismatches(self) -> float:
        n_ps = self.conditions.n_ps_linkages(self.oligo)
        return n_ps * (1.0 - self.conditions.sulfurization_efficiency)

    @property
    def fully_sulfurized_fraction(self) -> float:
        n_ps = self.conditions.n_ps_linkages(self.oligo)
        return self.conditions.sulfurization_efficiency**n_ps

    def po_mismatch_distribution(self, max_mismatches: int = 3) -> dict[int, float]:
        """Binomial distribution over count of PO-for-PS mismatches."""
        n_ps = self.conditions.n_ps_linkages(self.oligo)
        s = self.conditions.sulfurization_efficiency
        out: dict[int, float] = {}
        for k in range(0, min(max_mismatches, n_ps) + 1):
            out[k] = comb(n_ps, k) * ((1 - s) ** k) * (s ** (n_ps - k))
        return out

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
        lines = [
            f"Oligo            : {self.oligo.sequence_5to3}  (n={self.oligo.n})",
            f"Couplings        : {self.oligo.n_couplings}",
            "",
            f"Full-length (FLP): {self.full_length_fraction:.4%}",
            f"  naive c^(n-1)  : {self.naive_full_length_fraction:.4%}",
            f"Deletions        : {self.deletion_fraction:.4%}",
            f"Truncations      : {self.truncation_fraction:.4%}",
            f"Unresolved       : {self.unresolved_fraction:.4%}"
            f"  (>{self.max_deletions} deletions)",
            f"Mass balance     : {self.mass_balance:.10f}",
        ]
        n_ps = self.conditions.n_ps_linkages(self.oligo)
        if n_ps:
            lines += [
                "",
                f"PS linkages      : {n_ps}",
                f"Fully sulfurized : {self.fully_sulfurized_fraction:.4%}",
                f"E[PO mismatches] : {self.expected_po_mismatches:.4f}",
            ]
        lines += ["", "Top species:"]
        for s in self.top_species(8):
            lines.append(f"  {s.label:<22} {s.fraction:>9.4%}  {s.mass:>10.3f} Da")
        return "\n".join(lines)


def simulate(
    oligo: Oligo,
    conditions: ProcessConditions | None = None,
    max_deletions: int = 3,
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

    Returns
    -------
    SynthesisResult
    """
    if conditions is None:
        conditions = ProcessConditions()
    if max_deletions < 0:
        raise ValueError("max_deletions must be >= 0")

    # Active population, keyed by the set of positions already skipped.
    # Position 1 is preloaded on the support, so everything starts intact.
    active: dict[frozenset[int], float] = {frozenset(): 1.0}
    truncations: list[Species] = []
    unresolved = 0.0

    for position in range(2, oligo.n + 1):
        c = conditions.coupling_at(position, oligo)
        k = conditions.capping_efficiency
        nxt: dict[frozenset[int], float] = defaultdict(float)

        for deleted, amount in active.items():
            if amount == 0.0:
                continue

            # 1. Coupled: extends, deletion set unchanged.
            nxt[deleted] += amount * c

            failed = amount * (1.0 - c)
            if failed == 0.0:
                continue

            # 2. Capped: leaves the population as a truncation.
            capped = failed * k
            if capped > 0.0:
                length = (position - 1) - len(deleted)
                truncations.append(
                    Species(
                        deleted=deleted,
                        length=length,
                        fraction=capped,
                        truncated=True,
                        # A truncation is missing everything from `position` up.
                        mass=oligo.mass(
                            deleted | frozenset(range(position, oligo.n + 1))
                        ),
                    )
                )

            # 3. Escaped capping: stays active, gains a deletion here.
            escaped = failed * (1.0 - k)
            if escaped > 0.0:
                new_deleted = deleted | {position}
                if len(new_deleted) > max_deletions:
                    unresolved += escaped
                else:
                    nxt[frozenset(new_deleted)] += escaped

        active = dict(nxt)

    full_synthesis = [
        Species(
            deleted=deleted,
            length=oligo.n - len(deleted),
            fraction=amount,
            truncated=False,
            mass=oligo.mass(deleted),
        )
        for deleted, amount in active.items()
        if amount > 0.0
    ]

    return SynthesisResult(
        oligo=oligo,
        conditions=conditions,
        species=full_synthesis + truncations,
        unresolved_fraction=unresolved,
        max_deletions=max_deletions,
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
