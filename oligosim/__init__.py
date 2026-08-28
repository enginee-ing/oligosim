"""oligosim - a forward process model for solid-phase oligonucleotide synthesis.

Quick start
-----------
>>> from oligosim import Oligo, ProcessConditions, Sugar, Linkage, simulate
>>> oligo = Oligo.from_string("TCACTTTCATAATGCTGG", sugar=Sugar.MOE,
...                           linkage=Linkage.PS)
>>> result = simulate(oligo, ProcessConditions(coupling_efficiency=0.992))
>>> print(result.summary())
"""

from .chemistry import Base, Linkage, Oligo, Residue, Sugar
from .conditions import ProcessConditions
from .deprotection import (
    DeprotectedSpecies,
    DeprotectionConditions,
    DeprotectionResult,
    deprotect,
)
from .synthesis import Species, SynthesisResult, simulate, state_space_size

__version__ = "0.2.0.dev0"

__all__ = [
    "Base",
    "DeprotectedSpecies",
    "DeprotectionConditions",
    "DeprotectionResult",
    "Linkage",
    "Oligo",
    "ProcessConditions",
    "Residue",
    "Species",
    "Sugar",
    "SynthesisResult",
    "deprotect",
    "simulate",
    "state_space_size",
    "__version__",
]
