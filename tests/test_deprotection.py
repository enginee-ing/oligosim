"""Correctness tests for the deamination stage of oligosim/deprotection.py.

Deamination is modelled as a single post-synthesis expansion of each
`SynthesisResult` species by deamination count, binomial over that
species' own count of C/5-methyl-C residues. These tests pin: exact
reproduction at zero probability, the binomial shape, the mass delta,
mass-balance closure with overflow, and the per-species (not per-oligo)
eligible-residue count via nusinersen's four 5-methyl-C positions.
"""

import math

import pytest

from oligosim import (
    DeprotectionConditions,
    Oligo,
    ProcessConditions,
    Sugar,
    Linkage,
    deprotect,
    simulate,
)
from oligosim.chemistry import DEAMINATION_DELTA

TOL = 1e-12


def test_zero_probability_reproduces_input_exactly():
    """Default (and explicit-zero) deamination probability is a no-op."""
    oligo = Oligo.from_string("ACGTACGTACGT")
    synth = simulate(
        oligo, ProcessConditions(coupling_efficiency=0.95, capping_efficiency=0.8)
    )
    result = deprotect(synth)

    assert len(result.species) == len(synth.species)
    for expanded, parent in zip(result.species, synth.species):
        assert expanded.parent is parent
        assert expanded.deaminations == 0
        assert expanded.fraction == parent.fraction
        assert expanded.mass == parent.mass
    assert result.unresolved_fraction == synth.unresolved_fraction
    assert result.mass_balance == pytest.approx(synth.mass_balance, abs=TOL)


def test_deamination_count_is_binomial_over_eligible_residues():
    """A species with exactly 3 eligible C residues expands as
    Binomial(3, p) over deamination count, when the cap can't bind."""
    oligo = Oligo.from_string("ACCCG")  # 3 C's by construction
    synth = simulate(oligo, ProcessConditions(coupling_efficiency=1.0))
    assert len(synth.species) == 1  # perfect coupling: one species

    p = 0.1
    result = deprotect(
        synth,
        DeprotectionConditions(deamination_probability_per_residue=p),
        max_deaminations=3,
    )

    n_eligible = 3
    fractions_by_j = {s.deaminations: s.fraction for s in result.species}
    assert set(fractions_by_j) == set(range(n_eligible + 1))
    assert sum(fractions_by_j.values()) == pytest.approx(1.0, abs=TOL)
    for j in range(n_eligible + 1):
        expected = math.comb(n_eligible, j) * (p**j) * ((1 - p) ** (n_eligible - j))
        assert fractions_by_j[j] == pytest.approx(expected, abs=TOL)


def test_mass_shift_is_deamination_count_times_delta():
    """DEAMINATION_DELTA is +0.98402 Da monoisotopic, and every expanded
    species' mass reflects exactly n_deaminations * that delta."""
    assert DEAMINATION_DELTA == pytest.approx(0.98402, abs=1e-5)

    oligo = Oligo.from_string("ACCCG")
    synth = simulate(oligo, ProcessConditions(coupling_efficiency=1.0))
    parent = synth.species[0]
    result = deprotect(
        synth,
        DeprotectionConditions(deamination_probability_per_residue=0.1),
        max_deaminations=3,
    )
    for s in result.species:
        assert s.mass == pytest.approx(
            parent.mass + s.deaminations * DEAMINATION_DELTA, abs=TOL
        )


def test_mass_balance_closes_with_overflow():
    """Overflow past max_deaminations is reported, not dropped."""
    oligo = Oligo.from_string("ACGTACGTACGTACGT")  # 4 C's
    synth = simulate(
        oligo, ProcessConditions(coupling_efficiency=0.95, capping_efficiency=0.8)
    )
    result = deprotect(
        synth,
        DeprotectionConditions(deamination_probability_per_residue=0.5),
        max_deaminations=1,
    )
    assert result.unresolved_fraction > synth.unresolved_fraction
    assert result.mass_balance == pytest.approx(1.0, abs=TOL)


def test_correct_product_fraction_requires_zero_deaminations():
    """A deaminated FLP is correct length and correct backbone but has the
    wrong base -- a sequence error, not product. DeprotectionResult's
    correct_product_fraction must exclude it, on top of the synthesis-level
    deleted/truncated/mismatches checks: with every cap set wide enough
    that nothing overflows, it equals the full closed form
    prod(d_i * c_i) * s^n_ps * (1-p)^n_eligible."""
    n, c, d, s, p = 6, 0.99, 0.995, 0.995, 0.01
    oligo = Oligo.from_string("C" * n, linkage=Linkage.PS)  # every residue eligible
    n_ps = n - 1
    n_eligible = n
    synth = simulate(
        oligo,
        ProcessConditions(
            coupling_efficiency=c,
            detritylation_efficiency=d,
            sulfurization_efficiency=s,
        ),
        max_deletions=n - 1,
        max_mismatches=n_ps,
    )
    result = deprotect(
        synth,
        DeprotectionConditions(deamination_probability_per_residue=p),
        max_deaminations=n_eligible,
    )
    expected = (d * c) ** (n - 1) * s**n_ps * (1 - p) ** n_eligible
    assert result.correct_product_fraction == pytest.approx(expected, abs=TOL)
    # ...strictly less than the synthesis-stage figure, since deamination
    # can only remove mass from "correct product", never add to it.
    assert result.correct_product_fraction < synth.correct_product_fraction


def test_nusinersen_has_four_mc_positions_and_expands_over_full_range():
    """Eligible count is computed per species, not once for the parent
    oligo: nusinersen's correct-product species has exactly its 4
    5-methyl-C positions available, and expands over 0..4 deaminations."""
    oligo = Oligo.from_string(
        "TmCAmCTTTmCATAATGmCTGG", sugar=Sugar.MOE, linkage=Linkage.PS
    )
    synth = simulate(oligo)
    correct_product = next(
        s
        for s in synth.species
        if not s.deleted and not s.truncated and s.mismatches == 0
    )

    result = deprotect(
        synth,
        DeprotectionConditions(deamination_probability_per_residue=0.02),
        max_deaminations=4,
    )
    expansions = [s for s in result.species if s.parent is correct_product]
    assert {s.deaminations for s in expansions} == {0, 1, 2, 3, 4}
