"""Correctness tests for the v0.1 synthesis engine.

The important ones are the analytic invariants. If the engine is right, the
full-length fraction must equal prod(c_i) exactly and must not depend on
capping efficiency at all -- capping decides whether a *failure* becomes a
truncation or a deletion, but a chain is only full-length if it coupled on
every single cycle. Any bug in the population bookkeeping breaks one of those.
"""

import math

import pytest

from oligosim import (
    Base,
    Linkage,
    Oligo,
    ProcessConditions,
    Residue,
    Sugar,
    simulate,
    state_space_size,
)
from oligosim.chemistry import H2O

TOL = 1e-12


# ---------------------------------------------------------------------------
# Chemistry
# ---------------------------------------------------------------------------
def test_dna_residue_masses_match_published_values():
    """Standard monoisotopic DNA residue masses, to 4 dp."""
    expected = {Base.A: 313.0576, Base.C: 289.0464, Base.G: 329.0525, Base.T: 304.0460}
    for base, mass in expected.items():
        assert Residue(base).mass == pytest.approx(mass, abs=1e-4)


def test_ps_linkage_adds_sulfur_for_oxygen():
    po = Residue(Base.T, Sugar.DNA, Linkage.PO)
    ps = Residue(Base.T, Sugar.DNA, Linkage.PS)
    assert ps.mass - po.mass == pytest.approx(15.9772, abs=1e-4)


def test_moe_delta_is_consistent_via_two_routes():
    """2'-MOE from DNA must equal 2'-MOE from RNA plus the extra C2H4O."""
    dna = Residue(Base.A, Sugar.DNA).mass
    rna = Residue(Base.A, Sugar.RNA).mass
    moe = Residue(Base.A, Sugar.MOE).mass
    # RNA -> MOE replaces the 2'-OH hydrogen with CH2CH2OCH3 (C3H7O - H)
    assert moe - rna == pytest.approx(58.0419, abs=1e-4)
    assert moe - dna == pytest.approx(74.0368, abs=1e-4)


def test_oligo_mass_is_residues_plus_water_minus_terminal_phosphate():
    """An n-mer carries n-1 phosphates, not n."""
    from oligosim.chemistry import HPO3

    oligo = Oligo.from_string("ACGT")
    expected = sum(r.mass for r in oligo.residues) + H2O - HPO3
    assert oligo.mass() == pytest.approx(expected, abs=TOL)


def test_deletion_removes_exactly_one_residue_mass():
    oligo = Oligo.from_string("ACGT")
    # Position 2 counted from the 3' end.
    delta = oligo.mass() - oligo.mass(frozenset({2}))
    assert delta == pytest.approx(oligo.residue_at(2).mass, abs=TOL)


def test_sequence_roundtrip_and_orientation():
    oligo = Oligo.from_string("ACGT")
    assert oligo.sequence_5to3 == "ACGT"
    # Synthesis order is 3'->5', so position 1 is the 3'-terminal T.
    assert oligo.residue_at(1).base is Base.T
    assert oligo.residue_at(4).base is Base.A


def test_per_position_modifications():
    """Gapmer pattern: MOE wings, DNA core."""
    sugars = [Sugar.MOE] * 5 + [Sugar.DNA] * 10 + [Sugar.MOE] * 5
    oligo = Oligo.from_string("A" * 20, sugar=sugars, linkage=Linkage.PS)
    assert oligo.residue_at(20).sugar is Sugar.MOE  # 5' end
    assert oligo.residue_at(10).sugar is Sugar.DNA  # core
    assert oligo.residue_at(1).sugar is Sugar.MOE  # 3' end


def test_mixed_case_and_methylcytosine_parsing():
    oligo = Oligo.from_string("acmCg")
    assert [r.base for r in reversed(oligo.residues)] == [
        Base.A,
        Base.C,
        Base.mC,
        Base.G,
    ]


# ---------------------------------------------------------------------------
# Engine invariants
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("capping", [0.0, 0.5, 0.95, 1.0])
def test_full_length_equals_product_of_coupling_regardless_of_capping(capping):
    """The central invariant. FLP requires success at every cycle."""
    oligo = Oligo.from_string("ACGTACGTACGTACGTACGT")  # 20-mer, 19 couplings
    c = 0.99
    result = simulate(
        oligo,
        ProcessConditions(coupling_efficiency=c, capping_efficiency=capping),
        max_deletions=19,
    )
    assert result.full_length_fraction == pytest.approx(c**19, abs=1e-12)
    assert result.full_length_fraction == pytest.approx(
        result.naive_full_length_fraction, abs=1e-12
    )


@pytest.mark.parametrize("capping", [0.0, 0.7, 1.0])
@pytest.mark.parametrize("coupling", [0.95, 0.99, 0.999])
def test_mass_balance_is_conserved(capping, coupling):
    oligo = Oligo.from_string("ACGTACGTACGTACGT")
    result = simulate(
        oligo,
        ProcessConditions(
            coupling_efficiency=coupling, capping_efficiency=capping
        ),
        max_deletions=3,
    )
    assert result.mass_balance == pytest.approx(1.0, abs=1e-12)


def test_perfect_capping_produces_no_deletions():
    oligo = Oligo.from_string("ACGTACGTACGT")
    result = simulate(
        oligo, ProcessConditions(coupling_efficiency=0.98, capping_efficiency=1.0)
    )
    assert result.deletion_fraction == pytest.approx(0.0, abs=TOL)
    assert result.unresolved_fraction == pytest.approx(0.0, abs=TOL)


def test_no_capping_produces_no_truncations():
    oligo = Oligo.from_string("ACGTACGTACGT")
    result = simulate(
        oligo,
        ProcessConditions(coupling_efficiency=0.98, capping_efficiency=0.0),
        max_deletions=11,
    )
    assert result.truncation_fraction == pytest.approx(0.0, abs=TOL)


def test_perfect_coupling_gives_pure_product():
    oligo = Oligo.from_string("ACGTACGT")
    result = simulate(oligo, ProcessConditions(coupling_efficiency=1.0))
    assert result.full_length_fraction == pytest.approx(1.0, abs=TOL)
    assert len(result.species) == 1


def test_capping_shifts_failures_from_deletions_to_truncations():
    """The chemistry that justifies the capping step existing at all."""
    oligo = Oligo.from_string("ACGTACGTACGTACGTACGT")
    cond = lambda k: ProcessConditions(  # noqa: E731
        coupling_efficiency=0.99, capping_efficiency=k
    )
    low = simulate(oligo, cond(0.50), max_deletions=4)
    high = simulate(oligo, cond(0.99), max_deletions=4)

    assert high.deletion_fraction < low.deletion_fraction
    assert high.truncation_fraction > low.truncation_fraction
    # ...while leaving FLP untouched.
    assert high.full_length_fraction == pytest.approx(
        low.full_length_fraction, abs=1e-12
    )


def test_single_deletion_fraction_matches_closed_form():
    """One deletion at position p: fail there, escape capping, couple elsewhere."""
    n, c, k = 10, 0.98, 0.9
    oligo = Oligo.from_string("A" * n)
    result = simulate(
        oligo, ProcessConditions(coupling_efficiency=c, capping_efficiency=k)
    )
    by_pos = result.deletions_by_position()
    # n-1 couplings at positions 2..n; failing exactly at p means coupling
    # succeeded at the other n-2 cycles.
    expected = (1 - c) * (1 - k) * c ** (n - 2)
    for position in range(2, n + 1):
        assert by_pos[position] == pytest.approx(expected, abs=1e-12)


def test_truncation_lengths_are_physical():
    oligo = Oligo.from_string("ACGTACGTAC")
    result = simulate(
        oligo, ProcessConditions(coupling_efficiency=0.95, capping_efficiency=0.9)
    )
    for s in result.species:
        assert 1 <= s.length <= oligo.n
        if s.truncated:
            assert s.length < oligo.n


def test_species_masses_agree_with_lengths():
    oligo = Oligo.from_string("ACGTACGT", sugar=Sugar.OME, linkage=Linkage.PS)
    result = simulate(oligo)
    flp = next(
        s
        for s in result.species
        if not s.deleted and not s.truncated and s.mismatches == 0
    )
    assert flp.mass == pytest.approx(oligo.mass(), abs=TOL)
    for s in result.species:
        assert s.mass < oligo.mass() + TOL


# ---------------------------------------------------------------------------
# Truncation of the state space
# ---------------------------------------------------------------------------
def test_unresolved_fraction_is_reported_not_dropped():
    oligo = Oligo.from_string("A" * 30)
    result = simulate(
        oligo,
        ProcessConditions(coupling_efficiency=0.90, capping_efficiency=0.10),
        max_deletions=1,
    )
    assert result.unresolved_fraction > 0.0
    assert result.mass_balance == pytest.approx(1.0, abs=1e-12)


def test_raising_max_deletions_shrinks_unresolved():
    oligo = Oligo.from_string("A" * 25)
    cond = ProcessConditions(coupling_efficiency=0.95, capping_efficiency=0.2)
    coarse = simulate(oligo, cond, max_deletions=1)
    fine = simulate(oligo, cond, max_deletions=4)
    assert fine.unresolved_fraction < coarse.unresolved_fraction
    assert fine.full_length_fraction == pytest.approx(
        coarse.full_length_fraction, abs=1e-12
    )


def test_state_space_size_helper():
    assert state_space_size(20, 0) == 1
    assert state_space_size(20, 1) == 1 + 19
    assert state_space_size(20, 3) == 1 + 19 + math.comb(19, 2) + math.comb(19, 3)


# ---------------------------------------------------------------------------
# Sulfurization marginal
# ---------------------------------------------------------------------------
def test_ps_linkage_count_excludes_support_position():
    oligo = Oligo.from_string("ACGTACGT", linkage=Linkage.PS)
    cond = ProcessConditions()
    assert cond.n_ps_linkages(oligo) == oligo.n - 1


def test_po_mismatch_distribution_is_a_binomial():
    """With coupling certain and few enough PS linkages to stay within the
    mismatch cap, the joint-derived distribution is an exact binomial."""
    oligo = Oligo.from_string("ACGT", linkage=Linkage.PS)  # n_ps = 3, at the cap
    s = 0.99
    result = simulate(
        oligo, ProcessConditions(coupling_efficiency=1.0, sulfurization_efficiency=s)
    )
    dist = result.po_mismatch_distribution()
    n_ps = 3
    assert sum(dist.values()) == pytest.approx(1.0, abs=1e-12)
    for k in range(n_ps + 1):
        expected = math.comb(n_ps, k) * ((1 - s) ** k) * (s ** (n_ps - k))
        assert dist[k] == pytest.approx(expected, abs=1e-12)


def test_no_ps_linkages_when_all_po():
    oligo = Oligo.from_string("ACGT", linkage=Linkage.PO)
    result = simulate(oligo)
    dist = result.po_mismatch_distribution()
    assert set(dist) == {0}
    assert dist[0] == pytest.approx(1.0, abs=1e-12)


def test_full_length_fraction_ignores_mismatches():
    """A chain is full-length by virtue of length, regardless of PS/PO
    backbone composition -- mismatches split the FLP population but don't
    remove it. Uses n_ps == the mismatch cap so overflow can't muddy it.
    """
    oligo = Oligo.from_string("ACGT", linkage=Linkage.PS)  # n_ps = 3, at the cap
    result = simulate(
        oligo, ProcessConditions(coupling_efficiency=1.0, sulfurization_efficiency=0.9)
    )
    assert result.full_length_fraction == pytest.approx(1.0, abs=1e-12)
    # ...but is actually split across several mismatch counts.
    mismatch_counts = {s.mismatches for s in result.species if not s.truncated}
    assert len(mismatch_counts) > 1


def test_mismatch_overflow_is_reported_not_dropped():
    oligo = Oligo.from_string("A" * 12, linkage=Linkage.PS)  # 11 PS linkages
    result = simulate(
        oligo,
        ProcessConditions(coupling_efficiency=1.0, sulfurization_efficiency=0.9),
        max_mismatches=3,
    )
    assert result.unresolved_fraction > 0.0
    assert result.mass_balance == pytest.approx(1.0, abs=1e-12)


def test_correct_product_fraction_matches_closed_form_when_cap_not_binding():
    """With max_mismatches set to n_ps, the cap can never bind (the most
    mismatches any chain can carry equals n_ps), so the fully-correct-
    molecule fraction is exactly the closed-form product of the per-cycle
    coupling probability and the per-linkage sulfurization probability."""
    n, c, s = 10, 0.99, 0.995
    oligo = Oligo.from_string("A" * n, linkage=Linkage.PS)
    n_ps = n - 1
    result = simulate(
        oligo,
        ProcessConditions(coupling_efficiency=c, sulfurization_efficiency=s),
        max_mismatches=n_ps,
    )
    expected = c ** (n - 1) * s**n_ps
    assert result.correct_product_fraction == pytest.approx(expected, abs=1e-12)
    assert result.fully_sulfurized_fraction == pytest.approx(s**n_ps, abs=1e-12)


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------
def test_rejects_out_of_range_efficiencies():
    with pytest.raises(ValueError):
        ProcessConditions(coupling_efficiency=1.2)
    with pytest.raises(ValueError):
        ProcessConditions(capping_efficiency=-0.1)


def test_rejects_unparseable_base():
    with pytest.raises(ValueError):
        Oligo.from_string("ACXT")


def test_rejects_too_short_oligo():
    with pytest.raises(ValueError):
        Oligo.from_string("A")


def test_rejects_mismatched_modification_list():
    with pytest.raises(ValueError):
        Oligo.from_string("ACGT", sugar=[Sugar.MOE, Sugar.DNA])


def test_coupling_override_applies():
    oligo = Oligo.from_string("ACGTACGT")
    cond = ProcessConditions(coupling_efficiency=0.99, coupling_overrides={5: 0.80})
    assert cond.coupling_at(5, oligo) == 0.80
    assert cond.coupling_at(4, oligo) == 0.99
    result = simulate(oligo, cond)
    expected = 0.99**6 * 0.80
    assert result.full_length_fraction == pytest.approx(expected, abs=1e-12)


# ---------------------------------------------------------------------------
# Validation against a real drug substance
# ---------------------------------------------------------------------------
def test_terminal_residue_carries_no_phosphate():
    """An n-mer has n-1 phosphates. Position 1 has a free 3'-OH."""
    oligo = Oligo.from_string("ACGT")
    residue_sum = sum(r.mass for r in oligo.residues) + H2O
    # The model must be exactly one HPO3 lighter than the naive residue sum.
    from oligosim.chemistry import HPO3

    assert oligo.mass() == pytest.approx(residue_sum - HPO3, abs=TOL)


def test_nusinersen_monoisotopic_mass():
    """Nusinersen: 18-mer, uniform 2'-O-MOE, full PS, 5-methyl-C.

    Reference formula C234H340N61O128P17S17 derived from nucleoside
    composition: 18 x (2'-deoxynucleoside + C3H6O2) condensed into 17
    phosphorothioate diesters (each linkage +H3PO3S -2H2O). Sequence is
    public. This is an end-to-end check on sugar deltas, the PS delta, the
    5-methyl-C base, and the terminal-phosphate correction simultaneously.
    """
    M = {
        "H": 1.00782503207,
        "C": 12.0,
        "N": 14.0030740048,
        "O": 15.9949146196,
        "P": 30.97376163,
        "S": 31.97207100,
    }
    formula = {"C": 234, "H": 340, "N": 61, "O": 128, "P": 17, "S": 17}
    reference = sum(M[e] * n for e, n in formula.items())

    nusinersen = Oligo.from_string(
        "TmCAmCTTTmCATAATGmCTGG", sugar=Sugar.MOE, linkage=Linkage.PS
    )
    assert nusinersen.n == 18
    assert nusinersen.bases_5to3 == "TCACTTTCATAATGCTGG"
    assert nusinersen.mass() == pytest.approx(reference, abs=1e-6)


def test_uniform_modification_renders_compactly():
    oligo = Oligo.from_string("ACGT", sugar=Sugar.MOE, linkage=Linkage.PS)
    assert oligo.sequence_5to3 == "ACGT [2'-MOE, PS]"
    assert Oligo.from_string("ACGT").sequence_5to3 == "ACGT"


def test_mixed_modification_renders_per_residue():
    sugars = [Sugar.MOE, Sugar.MOE, Sugar.DNA, Sugar.DNA]
    oligo = Oligo.from_string("ACGT", sugar=sugars)
    assert not oligo.is_uniformly_modified
    assert "2'-MOE-A" in oligo.sequence_5to3


def test_species_labels_are_unambiguous():
    oligo = Oligo.from_string("ACGTACGT")
    result = simulate(
        oligo, ProcessConditions(coupling_efficiency=0.95, capping_efficiency=0.8)
    )
    labels = {s.label for s in result.species}
    assert "FLP" in labels
    assert any(l.startswith("trunc:") and l.endswith("mer") for l in labels)
    assert any(l.startswith("del:n-1@") for l in labels)
