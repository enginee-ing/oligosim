"""Modification-aware representation of an oligonucleotide.

Design note
-----------
Therapeutic oligonucleotides are almost never unmodified DNA. Nusinersen is
2'-O-MOE throughout with a full phosphorothioate backbone; siRNA duplexes mix
2'-OMe and 2'-F; gapmers alternate wings and core. Any data model that treats
"a sequence" as a string of ACGT has to be rewritten the moment it meets a real
drug substance, so `Residue` carries base, sugar and 3'-linkage from the start.

Masses are monoisotopic and computed from elemental composition, so they are
independently checkable. A "residue mass" here is the mass contributed by one
nucleoside-3'-phosphate unit inside a chain (i.e. condensed, water already lost).
The mass of an assembled oligo is `sum(residue masses) + H2O`, with the 5'-OH
and 3'-OH accounted for by that water.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Iterable, Sequence as TypingSequence

# ---------------------------------------------------------------------------
# Elemental monoisotopic masses (CODATA / IUPAC values)
# ---------------------------------------------------------------------------
_M_H = 1.00782503207
_M_C = 12.0
_M_N = 14.0030740048
_M_O = 15.9949146196
_M_P = 30.97376163
_M_S = 31.97207100
_M_F = 18.99840322

H2O = 2 * _M_H + _M_O

# Metaphosphate group. A "residue" carries a 3'-phosphate; the 3'-terminal
# residue does not, because after cleavage from the support it has a free
# 3'-OH. Subtracting this converts a residue mass to a nucleoside mass.
HPO3 = _M_H + _M_P + 3 * _M_O

# Mass difference for replacing one linkage oxygen with sulfur (PO -> PS).
PS_DELTA = _M_S - _M_O

# Deamination: C -> U and 5-methyl-C -> T are both a net loss of NH and
# gain of O (both bases' elemental formulas above give the identical
# delta). Monoisotopic: O - N - H = 0.98402 Da. Sources quoting +0.9848 Da
# are using AVERAGE atomic masses (~15.9994 O, ~14.0067 N, ~1.00794 H),
# not monoisotopic -- don't conflate the two against literature/vendor
# software.
DEAMINATION_DELTA = _M_O - _M_N - _M_H


class Base(str, Enum):
    """Nucleobase. `T` is thymine; `U` is uracil; `mC` is 5-methylcytosine."""

    A = "A"
    C = "C"
    G = "G"
    T = "T"
    U = "U"
    mC = "mC"

    @property
    def canonical(self) -> str:
        """Plain Watson-Crick letter, ignoring base modification.

        Published sequences write nusinersen as TCACTTTCATAATGCTGG even though
        every C is 5-methylated, so base-level comparison has to normalise.
        """
        return "C" if self is Base.mC else self.value


class Sugar(str, Enum):
    """2' sugar modification."""

    DNA = "DNA"  # 2'-H
    RNA = "RNA"  # 2'-OH
    OME = "2'-OMe"
    F = "2'-F"
    MOE = "2'-MOE"
    LNA = "LNA"  # 2'-O,4'-C-methylene bridge


class Linkage(str, Enum):
    """Internucleotide linkage on the 3' side of this residue."""

    PO = "PO"  # phosphodiester
    PS = "PS"  # phosphorothioate


# Base-specific mass of the 2'-deoxynucleoside-3'-phosphate residue (condensed).
# e.g. dA = C10H12N5O5P
_DNA_RESIDUE_MASS = {
    Base.A: 10 * _M_C + 12 * _M_H + 5 * _M_N + 5 * _M_O + _M_P,
    Base.C: 9 * _M_C + 12 * _M_H + 3 * _M_N + 6 * _M_O + _M_P,
    Base.G: 10 * _M_C + 12 * _M_H + 5 * _M_N + 6 * _M_O + _M_P,
    Base.T: 10 * _M_C + 13 * _M_H + 2 * _M_N + 7 * _M_O + _M_P,
    Base.U: 9 * _M_C + 11 * _M_H + 2 * _M_N + 7 * _M_O + _M_P,
    Base.mC: 10 * _M_C + 14 * _M_H + 3 * _M_N + 6 * _M_O + _M_P,
}

# Sugar modification expressed as a delta from the 2'-deoxy residue.
#   RNA   : 2'-H -> 2'-OH                      +O
#   2'-F  : 2'-H -> 2'-F                       +F -H
#   2'-OMe: 2'-H -> 2'-OCH3                    +O +CH2
#   2'-MOE: 2'-H -> 2'-OCH2CH2OCH3             +C3H6O2
#   LNA   : 2'-O,4'-C-methylene bridge         +O +CH2 -2H  (relative to DNA)
_SUGAR_DELTA = {
    Sugar.DNA: 0.0,
    Sugar.RNA: _M_O,
    Sugar.F: _M_F - _M_H,
    Sugar.OME: _M_O + _M_C + 2 * _M_H,
    Sugar.MOE: 3 * _M_C + 6 * _M_H + 2 * _M_O,
    Sugar.LNA: _M_O + _M_C + 2 * _M_H - 2 * _M_H,
}

# Relative coupling reactivity, used as a default when the caller does not
# supply explicit per-cycle efficiencies. These are ORDINAL PLACEHOLDERS
# encoding the well-documented qualitative ordering (bulky 2' substituents and
# bridged sugars couple more slowly than DNA); they are NOT calibrated rate
# constants and must not be read as such. Replace with literature values before
# using this model for anything quantitative. See README "Parameter status".
_SUGAR_RELATIVE_REACTIVITY = {
    Sugar.DNA: 1.00,
    Sugar.RNA: 0.90,
    Sugar.F: 0.97,
    Sugar.OME: 0.95,
    Sugar.MOE: 0.92,
    Sugar.LNA: 0.88,
}


@dataclass(frozen=True)
class Residue:
    """One monomer in the chain.

    `linkage` describes the bond on the 3' side of this residue, i.e. the bond
    formed when this residue was coupled onto the growing chain. The 3'-terminal
    residue sits on the solid support and its `linkage` is the support linkage,
    which is not a phosphorothioate candidate.
    """

    base: Base
    sugar: Sugar = Sugar.DNA
    linkage: Linkage = Linkage.PO

    @property
    def mass(self) -> float:
        """Monoisotopic residue mass, including the 3' linkage phosphate."""
        m = _DNA_RESIDUE_MASS[self.base] + _SUGAR_DELTA[self.sugar]
        if self.linkage is Linkage.PS:
            m += PS_DELTA
        return m

    @property
    def nucleoside_mass(self) -> float:
        """Mass without the 3' phosphate, for the 3'-terminal residue.

        The residue at synthesis position 1 is attached to the support through
        a succinate linker. On cleavage it is left with a free 3'-OH, so it
        contributes a nucleoside rather than a nucleoside-3'-phosphate. Getting
        this wrong adds a spurious HPO3 (or HPO2S) to every predicted mass.
        """
        return _DNA_RESIDUE_MASS[self.base] + _SUGAR_DELTA[self.sugar] - HPO3

    @property
    def relative_reactivity(self) -> float:
        """Placeholder ordinal reactivity. See module note and README."""
        return _SUGAR_RELATIVE_REACTIVITY[self.sugar]

    def __str__(self) -> str:
        prefix = "" if self.sugar is Sugar.DNA else f"{self.sugar.value}-"
        return f"{prefix}{self.base.value}"


class Oligo:
    """An ordered chain of residues, indexed 5' -> 3' in the usual convention.

    Internally, positions are numbered 1..n from the 3' end, because that is
    the order in which solid-phase synthesis builds the chain: position 1 is
    loaded on the support and positions 2..n are added by n-1 coupling cycles.
    `Oligo.residues` is stored in synthesis order (3' -> 5') to keep the engine
    readable; `Oligo.sequence_5to3` gives the human-facing view.
    """

    def __init__(self, residues_3to5: TypingSequence[Residue]):
        if len(residues_3to5) < 2:
            raise ValueError("Oligo needs at least 2 residues")
        self._residues = tuple(residues_3to5)

    # -- constructors -------------------------------------------------------
    @classmethod
    def from_string(
        cls,
        sequence_5to3: str,
        sugar: Sugar | Iterable[Sugar] = Sugar.DNA,
        linkage: Linkage | Iterable[Linkage] = Linkage.PO,
    ) -> "Oligo":
        """Build from a 5'->3' base string with uniform or per-position mods.

        `sugar` and `linkage` accept either a single value (applied uniformly)
        or an iterable given in 5'->3' order, one entry per base.
        """
        bases = cls._parse_bases(sequence_5to3)
        n = len(bases)

        sugars = [sugar] * n if isinstance(sugar, Sugar) else list(sugar)
        linkages = [linkage] * n if isinstance(linkage, Linkage) else list(linkage)
        if len(sugars) != n or len(linkages) != n:
            raise ValueError("sugar/linkage iterables must match sequence length")

        residues_5to3 = [
            Residue(b, s, l) for b, s, l in zip(bases, sugars, linkages)
        ]
        return cls(tuple(reversed(residues_5to3)))

    @staticmethod
    def _parse_bases(seq: str) -> list[Base]:
        bases: list[Base] = []
        i = 0
        text = seq.strip().replace(" ", "")
        while i < len(text):
            if text[i : i + 2] == "mC":
                bases.append(Base.mC)
                i += 2
                continue
            try:
                bases.append(Base(text[i].upper()))
            except ValueError as exc:
                raise ValueError(f"Unrecognised base {text[i]!r} at index {i}") from exc
            i += 1
        return bases

    # -- views --------------------------------------------------------------
    @property
    def residues(self) -> tuple[Residue, ...]:
        """Residues in synthesis order (3' -> 5'), position 1 first."""
        return self._residues

    @property
    def n(self) -> int:
        return len(self._residues)

    @property
    def n_couplings(self) -> int:
        """Coupling cycles required: the 3'-terminal residue is preloaded."""
        return self.n - 1

    @property
    def bases_5to3(self) -> str:
        """Canonical base string, 5'->3', with modifications normalised."""
        return "".join(r.base.canonical for r in reversed(self._residues))

    @property
    def is_uniformly_modified(self) -> bool:
        first = self._residues[0]
        return all(
            r.sugar is first.sugar and r.linkage is first.linkage
            for r in self._residues
        )

    @property
    def sequence_5to3(self) -> str:
        """Human-facing sequence.

        Uniformly modified oligos render compactly, since annotating all 18
        residues of a fully 2'-MOE/PS drug substance is noise. Mixed-chemistry
        oligos (gapmers, siRNA strands) annotate per residue, where it matters.
        """
        if self.is_uniformly_modified:
            r = self._residues[0]
            tags = []
            if r.sugar is not Sugar.DNA:
                tags.append(r.sugar.value)
            if r.linkage is not Linkage.PO:
                tags.append(r.linkage.value)
            suffix = f" [{', '.join(tags)}]" if tags else ""
            return f"{self.bases_5to3}{suffix}"
        return "".join(str(r) for r in reversed(self._residues))

    def residue_at(self, position: int) -> Residue:
        """Residue at 1-based position counted from the 3' end."""
        return self._residues[position - 1]

    def mass(self, deleted: frozenset[int] = frozenset()) -> float:
        """Monoisotopic mass of the assembled chain, minus any deletions.

        `deleted` holds 1-based synthesis positions absent from the species.
        A deletion removes that residue's mass, including its 3' phosphate.

        The surviving 3'-most residue contributes a nucleoside mass (free
        3'-OH); every residue above it contributes a full residue mass. So an
        n-mer carries n-1 phosphates, as it must.
        """
        present = [i for i in range(1, self.n + 1) if i not in deleted]
        if not present:
            return 0.0
        terminal, rest = present[0], present[1:]
        total = self.residue_at(terminal).nucleoside_mass
        total += sum(self.residue_at(i).mass for i in rest)
        return total + H2O

    def __len__(self) -> int:
        return self.n

    def __repr__(self) -> str:
        return f"Oligo({self.sequence_5to3!r}, n={self.n})"
