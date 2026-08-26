"""Process parameters for one synthesis run.

In v0.1 these are supplied directly. In a later version the `amidite` module
will *derive* the per-cycle coupling efficiency from raw-material quality
attributes (water content, 31P purity, free acid, related substances) plus
activator choice, excess and coupling time. Keeping the per-cycle efficiency
behind `ProcessConditions.coupling_at()` means that substitution will not
require changing the synthesis engine.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping

from .chemistry import Linkage, Oligo


def _check_fraction(name: str, value: float) -> None:
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must lie in [0, 1], got {value}")


@dataclass
class ProcessConditions:
    """Cycle parameters, with optional per-position overrides.

    Parameters
    ----------
    coupling_efficiency
        Fraction of available 5'-OH that couples in a cycle. The single most
        important parameter: a 20-mer at 0.99 gives a very different crude
        profile than the same sequence at 0.995.
    capping_efficiency
        Fraction of *failed* chains acetylated and removed from the growing
        population. Capping does not improve yield; it converts would-be
        deletion sequences into truncations, which are far easier to separate.
        This is why the parameter matters more than its size suggests.
    sulfurization_efficiency
        Fraction of phosphite triesters converted to phosphorothioate at PS
        positions. Shortfall produces PO-for-PS mismatches, which are among
        the hardest impurities to resolve chromatographically.
    apply_sugar_reactivity
        If True, scale coupling efficiency by the incoming residue's ordinal
        sugar reactivity. Off by default because those factors are
        placeholders, not calibrated values.
    coupling_overrides
        Explicit per-position coupling efficiency, keyed by 1-based synthesis
        position (the position being *added*, so keys run 2..n).
    """

    coupling_efficiency: float = 0.992
    capping_efficiency: float = 0.95
    sulfurization_efficiency: float = 0.995
    apply_sugar_reactivity: bool = False
    coupling_overrides: Mapping[int, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        _check_fraction("coupling_efficiency", self.coupling_efficiency)
        _check_fraction("capping_efficiency", self.capping_efficiency)
        _check_fraction("sulfurization_efficiency", self.sulfurization_efficiency)
        for pos, val in self.coupling_overrides.items():
            _check_fraction(f"coupling_overrides[{pos}]", val)

    def coupling_at(self, position: int, oligo: Oligo) -> float:
        """Effective coupling efficiency for the cycle that adds `position`."""
        if position in self.coupling_overrides:
            return self.coupling_overrides[position]
        eff = self.coupling_efficiency
        if self.apply_sugar_reactivity:
            eff *= oligo.residue_at(position).relative_reactivity
        return min(eff, 1.0)

    def n_ps_linkages(self, oligo: Oligo) -> int:
        """Count PS linkages formed during synthesis.

        Position 1 sits on the support, so its `linkage` field describes the
        support attachment and is not a sulfurization target.
        """
        return sum(
            1
            for pos in range(2, oligo.n + 1)
            if oligo.residue_at(pos).linkage is Linkage.PS
        )
