"""Fast baseline tests of the two distinct applications."""

import ast
from pathlib import Path

import pytest

from stereoapps.core import (
    dependency_svg,
    detect_stereo_elements,
    describe_sites,
    exact_dependency_analysis,
    heuristic_edges,
    molecule_svg,
    parse_smiles,
    rdkit_stereoisomer_count,
    run_enumeration,
)
from stereoapps.examples import EXAMPLES


def test_both_apps_are_independent_python_entrypoints():
    for filename in ("stereocenter_detection.py", "stereoisomer_generation.py"):
        script = Path("apps") / filename
        assert script.is_file()
        ast.parse(script.read_text(encoding="utf-8"))


def test_simple_tetrahedral_detection_and_enumeration():
    bundle = parse_smiles("C[C@H](O)N", matched_rdkit_sites=True)
    assert len(bundle.retained_tetrahedral_atoms) == 1
    assert not bundle.missing_rdkit_candidates
    analysis = detect_stereo_elements(bundle.graph)
    rows = describe_sites(analysis)
    assert len(rows) == 1
    assert rows[0]["local_states"] == 2
    result = run_enumeration(bundle.graph)
    assert result["formal_assignments"] == 2
    assert result["unique_stereo_graphs"] == 2
    assert rdkit_stereoisomer_count(bundle.molecule) == 2


def test_no_descriptor_sites_are_inserted():
    bundle = parse_smiles("C[C@H](O)N", matched_rdkit_sites=True)
    original = set(bundle.graph.atom_stereo)
    detect_stereo_elements(bundle.graph)
    run_enumeration(bundle.graph)
    assert set(bundle.graph.atom_stereo) == original
    assert all(
        stereo.parity is None for stereo in bundle.graph.atom_stereo.values()
    )


def test_excluded_isotopes_and_multiple_bonds():
    with pytest.raises(ValueError, match="Isotope"):
        parse_smiles("[2H][C@@H](Br)Cl")
    with pytest.raises(ValueError, match="Only single-bond"):
        parse_smiles("C/C=C/C")


def test_inactive_tetrahedral_candidate():
    bundle = parse_smiles("CC(C)(F)Cl")
    if not bundle.retained_tetrahedral_atoms:
        pytest.skip("SMG import did not nominate an inactive tetrahedral center")
    result = run_enumeration(bundle.graph)
    assert result["unique_stereo_graphs"] == 1
    exact = exact_dependency_analysis(bundle.graph)
    assert all(
        x == "inactive in all tested contexts" for x in exact.activity.values()
    )


def test_complex_example_dependency_graph_rendering():
    bundle = parse_smiles(
        EXAMPLES["VS294 — bridged oxygen rings"],
        matched_rdkit_sites=True,
    )
    analysis = detect_stereo_elements(bundle.graph)
    assert len(analysis.sites) == 5
    vertices = set(bundle.retained_tetrahedral_atoms)
    for source, target in heuristic_edges(analysis):
        assert source in vertices
        assert target in vertices
        assert source != target
    svg = dependency_svg(analysis)
    assert "<svg " in svg or "<svg\n" in svg
    mol_svg = molecule_svg(bundle)
    assert "<svg" in mol_svg


def test_two_ordinary_sites_have_no_exact_dependency():
    molecule = "F[C@H](Br)[C@H](I)Cl"
    bundle = parse_smiles(molecule, matched_rdkit_sites=True)
    analysis = detect_stereo_elements(bundle.graph)
    assert len(analysis.sites) == 2
    exact = exact_dependency_analysis(bundle.graph, analysis)
    assert set(exact.activity.values()) == {"always active"}
    assert not exact.edges
