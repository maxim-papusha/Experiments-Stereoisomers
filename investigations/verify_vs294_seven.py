"""Independent investigation of the seven VS294/VS295 tetrahedral assignments.

Uses RDKit + NetworkX for a group-theoretical calculation independent of
StereoMolGraph's equality implementation, then cross-checks StereoMolGraph's
seven generated representatives. Also probes both C3 configurations using
ETKDG and MMFF without RDKit's normal stereochemical cleanup.

Run:
  python investigations/verify_vs294_seven.py --conformers 10 \
    --output investigations/results/vs294_analysis.json
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from itertools import product
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import networkx as nx
import numpy as np
from networkx.algorithms.isomorphism import GraphMatcher
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
    StereoEnumerationOptions,
)

from stereoapps.core import (
    detect_stereo_elements,
    enumerate_stereo_assignments,
    parse_smiles,
)
from stereoapps.examples import EXAMPLES
from stereomolgraph.algorithms.isomorphism import vf2pp_all_isomorphisms


CENTERS = (1, 3, 5, 8, 9)
SMILES = EXAMPLES["VS294 — bridged oxygen rings"]


def parity(permutation: tuple[int, ...] | list[int]) -> int:
    """Permutation sign, +1 for even and -1 for odd."""
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def reference_graph():
    molecule = Chem.MolFromSmiles(SMILES)
    if molecule is None:
        raise ValueError("Invalid example SMILES")
    Chem.RemoveStereochemistry(molecule)
    explicit = Chem.AddHs(molecule)
    graph = nx.Graph()
    for atom in explicit.GetAtoms():
        graph.add_node(atom.GetIdx(), element=atom.GetAtomicNum())
    for bond in explicit.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            raise ValueError("Example must contain only single bonds")
        graph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    atom_match = lambda a, b: a["element"] == b["element"]
    automorphisms = tuple(
        GraphMatcher(graph, graph, node_match=atom_match).isomorphisms_iter()
    )
    assert all(len(tuple(graph.neighbors(c))) == 4 for c in CENTERS)
    return molecule, explicit, graph, automorphisms


def action(automorphism, signs, slots):
    """Induced automorphism action on tetrahedral signs.

    slots[c] defines an ordered reference tuple of four neighbors of c.
    Mapping and parity transport preserve handedness under the atom
    permutation. This implementation does NOT call StereoMolGraph.
    """
    index = {center: i for i, center in enumerate(CENTERS)}
    out = [0] * len(CENTERS)
    for i, center in enumerate(CENTERS):
        mapped_center = automorphism[center]
        target = slots[mapped_center]
        mapped_slots = [automorphism[n] for n in slots[center]]
        permutation = [target.index(atom) for atom in mapped_slots]
        out[index[mapped_center]] = signs[i] * parity(permutation)
    return tuple(out)


def partition(automorphisms, slots):
    """All 2^5 local assignments, with exact symmetry orbit partition."""
    assignments = tuple(product((-1, 1), repeat=len(CENTERS)))
    orbit_of = {}
    orbits = []
    for assignment in assignments:
        if assignment in orbit_of:
            continue
        orbit = frozenset(
            action(auto, assignment, slots) for auto in automorphisms
        )
        idx = len(orbits)
        for other in orbit:
            if other in orbit_of:
                raise AssertionError("Orbits overlap!")
            orbit_of[other] = idx
        orbits.append(orbit)
    assert len(orbit_of) == len(assignments)
    fixed_counts = [
        sum(action(auto, x, slots) == x for x in assignments)
        for auto in automorphisms
    ]
    assert sum(fixed_counts) % len(automorphisms) == 0
    assert len(orbits) == sum(fixed_counts) // len(automorphisms)
    return assignments, orbits, orbit_of, fixed_counts


def rdkit_stereo_cleanup(molecule, rdkit_slots, automorphisms):
    """Group all fully tagged configurations by RDKit-cleaned SMILES."""
    assignments, orbits, orbit_of, fixes = partition(
        automorphisms, rdkit_slots
    )
    groups = defaultdict(list)
    center3_tags = Counter()
    before_cleanup = {}
    for signs in assignments:
        candidate = Chem.Mol(molecule)
        for center, sign in zip(CENTERS, signs, strict=True):
            candidate.GetAtomWithIdx(center).SetChiralTag(
                Chem.ChiralType.CHI_TETRAHEDRAL_CW
                if sign == 1 else Chem.ChiralType.CHI_TETRAHEDRAL_CCW
            )
        before_cleanup[signs] = Chem.MolToSmiles(
            candidate, isomericSmiles=True
        )
        Chem.AssignStereochemistry(candidate, cleanIt=True, force=True)
        center3_tags[str(candidate.GetAtomWithIdx(3).GetChiralTag())] += 1
        key = Chem.MolToSmiles(candidate, isomericSmiles=True)
        groups[key].append(signs)

    # Crucial validation: one full stereochemical graph orbit should never
    # lead to multiple distinct cleaned RDKit structures.
    keys_per_orbit = {
        i: [key for key, members in groups.items()
            if any(orbit_of[x] == i for x in members)]
        for i in range(len(orbits))
    }
    assert all(len(keys) == 1 for keys in keys_per_orbit.values())
    clean_groups = []
    for key, members in sorted(groups.items()):
        clean_groups.append({
            "smiles": key,
            "formal_assignments": len(members),
            "distinct_full_tetrahedral_orbits": sorted(
                {orbit_of[x] for x in members}
            ),
        })
    collapsed = [g for g in clean_groups
                 if len(g["distinct_full_tetrahedral_orbits"]) > 1]

    opts = StereoEnumerationOptions(
        onlyUnassigned=True, unique=True, maxIsomers=0,
        tryEmbedding=False,
    )
    library_count = sum(
        1 for _ in EnumerateStereoisomers(molecule, options=opts)
    )
    assert len(groups) == library_count
    return {
        "explicitly_tagged_orbits": len(orbits),
        "cleaned_rdkit_outputs": library_count,
        "center3_tags_after_cleanup": dict(center3_tags),
        "clean_groups": clean_groups,
        "collapsed_orbit_groups": collapsed,
        "orbits_from_full_configurations": [
            {
                "index": i,
                "representative": list(sorted(orbit)[0]),
                "size": len(orbit),
                "mirror_orbit": orbit_of[
                    tuple(-v for v in sorted(orbit)[0])
                ],
            }
            for i, orbit in enumerate(orbits)
        ],
        "fixed_counts": sorted(fixes),
        "example_raw_smiles_with_c3": list(before_cleanup.values())[0],
    }


def check_stereomolgraph(automorphisms, explicit):
    """Compare graph-theoretical orbits with the actual seven SMG outputs."""
    bundle = parse_smiles(SMILES, matched_rdkit_sites=True)
    graph = bundle.graph
    assert set(graph.atom_stereo) == set(CENTERS)
    assert not graph.bond_stereo
    slots = {
        center: tuple(graph.atom_stereo[center].atoms[1:])
        for center in CENTERS
    }
    assert all(
        len(v) == 4 and set(v) == set(
            a.GetIdx() for a in explicit.GetAtomWithIdx(c).GetNeighbors()
        )
        for c, v in slots.items()
    )
    _, orbits, orbit_of, fixes = partition(automorphisms, slots)
    analysis = detect_stereo_elements(graph)
    raw = tuple(enumerate_stereo_assignments(graph, analysis, unique=False))
    unique = tuple(enumerate_stereo_assignments(graph, analysis, unique=True))
    fingerprints = [
        tuple(int(g.atom_stereo[c].parity) for c in CENTERS)
        for g in unique
    ]
    ids = [orbit_of[bits] for bits in fingerprints]

    # This checks every pair against SMG's VF2 implementation without
    # the hash-bucket shortcut of Python set deduplication.
    direct_matches = []
    for i in range(len(unique)):
        for j in range(i + 1, len(unique)):
            isomorphic = any(
                vf2pp_all_isomorphisms(
                    unique[i], unique[j],
                    stereo=True, stereo_change=False, subgraph=False,
                )
            )
            if isomorphic:
                direct_matches.append([i, j])
    assert len(raw) == 32, len(raw)
    assert len(unique) == 7, len(unique)
    assert len(orbits) == 7
    assert len(ids) == len(set(ids)) == 7
    assert direct_matches == []
    assert len(automorphisms) == 8
    return {
        "raw_assignments": len(raw),
        "unique_smg_outputs": len(unique),
        "independent_orbit_count": len(orbits),
        "burnside_sum_fixed": sum(fixes),
        "graph_automorphisms": len(automorphisms),
        "unique_representative_orbit_ids": ids,
        "unique_representative_parities": [list(x) for x in fingerprints],
        "pairwise_vf2_matches_among_seven": direct_matches,
        "independent_orbit_sizes": sorted(len(orbit) for orbit in orbits),
    }


def signed_volume(mol, conformer_id, atom_index):
    """Oriented volume of four ligand positions in RDKit neighbor order."""
    positions = mol.GetConformer(conformer_id).GetPositions()
    atom = mol.GetAtomWithIdx(atom_index)
    neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
    x = positions[neighbors]
    return float(np.linalg.det(np.array([
        x[0] - x[3], x[1] - x[3], x[2] - x[3]
    ])))


def embedding_probe(molecule, n_conformers):
    """Compare two opposite C3 states with other four sites unchanged.

    A force-field-optimized conformer supports 3D local realizability,
    not isolation, optical activity or resistance to interconversion.
    """
    if n_conformers < 1:
        return {"skipped": True}
    reports = []
    for c3 in (-1, 1):
        signs = (-1, c3, -1, -1, -1)
        mol = Chem.Mol(molecule)
        for c, sign in zip(CENTERS, signs, strict=True):
            mol.GetAtomWithIdx(c).SetChiralTag(
                Chem.ChiralType.CHI_TETRAHEDRAL_CW
                if sign > 0 else Chem.ChiralType.CHI_TETRAHEDRAL_CCW
            )
        # Preserve the center's tag even if RDKit's cleanup considers
        # it nonstereogenic.
        Chem.AssignStereochemistry(mol, cleanIt=False, force=True)
        assert (mol.GetAtomWithIdx(3).GetChiralTag()
                != Chem.ChiralType.CHI_UNSPECIFIED)
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 21067
        params.numThreads = 1
        ids = list(AllChem.EmbedMultipleConfs(
            mol, numConfs=n_conformers, params=params
        ))
        props = AllChem.MMFFGetMoleculeProperties(mol)
        energies = []
        volumes = []
        for conformer_id in ids:
            ff = AllChem.MMFFGetMoleculeForceField(
                mol, props, confId=conformer_id
            )
            if ff is None:
                raise RuntimeError("MMFF unsupported for the molecule")
            ff.Minimize(maxIts=500)
            e = float(ff.CalcEnergy())
            v = [signed_volume(mol, conformer_id, c) for c in CENTERS]
            if all(abs(value) > 0.1 for value in v):
                energies.append(e)
                volumes.append(v)

        reports.append({
            "central_tag": c3,
            "requested_tags": list(signs),
            "embedded": len(ids),
            "nondegenerate_optimized": len(volumes),
            "min_mmff_energy_kcal_mol": min(energies) if energies else None,
            "max_mmff_energy_kcal_mol": max(energies) if energies else None,
            "oriented_volume_signatures": [
                [int(math.copysign(1, v)) for v in vs]
                for vs in volumes
            ],
            "first_volume_values": volumes[0] if volumes else None,
        })
    assert all(x["nondegenerate_optimized"] == n_conformers
               for x in reports)
    # All unchanged centers retain their chirality; C3 flips.
    for other in reports[1:]:
        assert all(
            left[0] == right[0]
            and left[2:] == right[2:]
            and left[1] == -right[1]
            for left, right in zip(
                reports[0]["oriented_volume_signatures"],
                other["oriented_volume_signatures"],
            )
        )
    return {"conformers_per_configuration": n_conformers,
            "two_c3_configurations": reports}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--conformers", type=int, default=10)
    parser.add_argument(
        "--output", type=Path,
        default=Path("investigations/results/vs294_analysis.json")
    )
    args = parser.parse_args()

    rd_mol, explicit, _, automorphisms = reference_graph()
    rd_slots = {
        c: tuple(n.GetIdx() for n in explicit.GetAtomWithIdx(c).GetNeighbors())
        for c in CENTERS
    }
    result = {
        "molecule": "VS294 (same constitution as VS295)",
        "smiles": SMILES,
        "rdkit_version": rdBase.rdkitVersion,
        "centers": list(CENTERS),
        "independent_group": {
            "automorphism_count": len(automorphisms)
        },
        "stereomolgraph": check_stereomolgraph(
            automorphisms, explicit
        ),
        "rdkit": rdkit_stereo_cleanup(
            rd_mol, rd_slots, automorphisms
        ),
        "3d_feasibility_probe": embedding_probe(
            rd_mol, args.conformers
        ),
    }
    assert result["stereomolgraph"]["unique_smg_outputs"] == 7
    assert result["rdkit"]["cleaned_rdkit_outputs"] == 6
    assert sum(len(x["distinct_full_tetrahedral_orbits"]) > 1
               for x in result["rdkit"]["clean_groups"]) == 1
    print(json.dumps(result, indent=2))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print("RESULT_JSON=" + str(args.output))


if __name__ == "__main__":
    main()
