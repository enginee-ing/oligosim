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
from .synthesis import Species, SynthesisResult, simulate, state_space_size

__version__ = "0.1.0"

__all__ = [
    "Base",
    "Linkage",
    "Oligo",
    "ProcessConditions",
    "Residue",
    "Species",
    "Sugar",
    "SynthesisResult",
    "simulate",
    "state_space_size",
    "__version__",
]
