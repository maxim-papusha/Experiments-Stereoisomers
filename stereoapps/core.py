"""Shared, importable computations for the two independent applications.

The *enumeration* functions NEVER introduce stereo descriptor sites.
All candidate elements come from the supplied StereoMolGraph. The optional
SMILES adapter is an upstream convenience for preparing such an input.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from math import prod
from pathlib import Path
from typing import Literal

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
    StereoEnumerationOptions,
)

from stereomolgraph import StereoMolGraph
from stereomolgraph.experimental.stereoisomer_draft import (
    StereoAnalysis,
    detect_stereo_elements,
    enumerate_stereo_assignments,
    locally_stereogenic_at_assignment,
)
from stereomolgraph.stereodescriptors import Tetrahedral


@dataclass(frozen=True)
class InputBundle:
    graph: StereoMolGraph
    molecule: Chem.Mol
    smiles: str
    rdkit_candidate_atoms: frozenset[int]
    imported_tetrahedral_atoms: frozenset[int]
    retained_tetrahedral_atoms: frozenset[int]
    missing_rdkit_candidates: frozenset[int]


@dataclass(frozen=True)
class ExactDependencies:
    """Exact with respect to the SMG isomorphism model, not physical stability.

    The directed edge t -> s means there is a completed assignment of other
    stereo descriptors for which changing t changes the local activity of s.
    """

    edges: frozenset[tuple[int, int]]
    activity: dict[int, str]
    evaluations: int
    witnesses: dict[str, dict[str, object]]


def parse_smiles(
    smiles: str,
    *,
    unassign: bool = True,
    matched_rdkit_sites: bool = False,
    restrict_single_bonds: bool = True,
) -> InputBundle:
    """Create a tetrahedral-only SMG from a SMILES input.

    The importer nominates geometry sites first. Here non-tetrahedral
    sites are removed because this application is intentionally restricted
    to tetrahedral stereochemistry. No descriptor is ever added by the
    analysis or enumeration engine.

    matched_rdkit_sites=True is for *controlled RDKit comparisons*, not
    for independent detection: it filters existing SMG sites to the RDKit
    potential tetrahedral candidate set, without adding missing sites.
    """
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("RDKit cannot parse the input SMILES")
    if any(atom.GetIsotope() for atom in molecule.GetAtoms()):
        raise ValueError("Isotope-labelled atoms are outside the present scope")
    if restrict_single_bonds and any(
        bond.GetBondType() != Chem.BondType.SINGLE
        for bond in molecule.GetBonds()
    ):
        raise ValueError(
            "Only single-bond molecules are supported for this comparison "
            "(double, triple and aromatic bonds are excluded)."
        )

    if unassign:
        Chem.RemoveStereochemistry(molecule)

    potential = Chem.FindPotentialStereo(molecule)
    rdkit_sites = frozenset(
        int(stereo.centeredOn)
        for stereo in potential
        if stereo.type == Chem.StereoType.Atom_Tetrahedral
    )

    graph = StereoMolGraph.from_rdmol(
        Chem.AddHs(Chem.Mol(molecule)),
        stereo_complete=False,
        resonance=False,
    )

    # This is an explicit *application adapter*, not part of the two
    # analysis algorithms. Never add an atom/bond descriptor.
    for bond in tuple(graph.bond_stereo):
        graph.delete_bond_stereo(bond)
    for atom, stereo in tuple(graph.atom_stereo.items()):
        if not isinstance(stereo, Tetrahedral):
            graph.delete_atom_stereo(atom)

    imported = frozenset(int(atom) for atom in graph.atom_stereo)
    if matched_rdkit_sites:
        for atom in tuple(graph.atom_stereo):
            if atom not in rdkit_sites:
                graph.delete_atom_stereo(atom)

    retained = frozenset(int(atom) for atom in graph.atom_stereo)
    return InputBundle(
        graph=graph,
        molecule=molecule,
        smiles=smiles,
        rdkit_candidate_atoms=rdkit_sites,
        imported_tetrahedral_atoms=imported,
        retained_tetrahedral_atoms=retained,
        missing_rdkit_candidates=rdkit_sites - imported,
    )


def describe_sites(analysis: StereoAnalysis) -> list[dict[str, object]]:
    """Stable JSON-friendly table, restricted to tetrahedral descriptors."""
    rows = []
    for site in analysis.sites:
        if not isinstance(site.template, Tetrahedral):
            continue
        atom_id = int(site.key[1])
        rows.append(
            {
                "atom": atom_id,
                "configuration": site.template.parity,
                "ligands": list(site.template.atoms[1:]),
                "classification": site.classification,
                "local_states": len(site.states),
                "candidate_influences": sorted(
                    int(source[1]) for source in analysis.dependencies[site.key]
                    if source[0] == "atom"
                ),
            }
        )
    return sorted(rows, key=lambda row: int(row["atom"]))


def heuristic_edges(analysis: StereoAnalysis) -> frozenset[tuple[int, int]]:
    """Return directed atom-site edges (source, target). Not proof."""
    return frozenset(
        (int(source[1]), int(target[1]))
        for target, influences in analysis.dependencies.items()
        if target[0] == "atom"
        for source in influences
        if source[0] == "atom"
    )


def exact_dependency_analysis(
    graph: StereoMolGraph,
    analysis: StereoAnalysis | None = None,
    *,
    max_sites: int = 5,
    max_assignments: int = 32,
) -> ExactDependencies:
    """Enumerate every represented context and test site-anchored activity.

    This is an exact (and potentially costly) dependency test *relative to
    the current SMG stereo-isomorphism implementation*. No approximate
    dependence edge is used to discard configurations.

    Every site must be a supplied tetrahedral descriptor. The method does
    not add, remove or infer stereo sites.
    """
    if analysis is None:
        analysis = detect_stereo_elements(graph)
    sites = tuple(analysis.sites)
    if any(
        site.key[0] != "atom" or not isinstance(site.template, Tetrahedral)
        for site in sites
    ):
        raise ValueError("Only existing tetrahedral atom sites are supported")
    if len(sites) > max_sites:
        raise ValueError(
            f"{len(sites)} sites exceed the exact analysis limit {max_sites}"
        )
    domain_sizes = tuple(len(site.states) for site in sites)
    n_assignments = prod(domain_sizes)
    if n_assignments > max_assignments:
        raise ValueError(
            f"{n_assignments} assignments exceed the limit {max_assignments}"
        )

    configs = tuple(product(*(range(n) for n in domain_sizes)))
    activity: dict[tuple[int, tuple[int, ...]], bool] = {}

    for config in configs:
        for target_index, site in enumerate(sites):
            # A_s depends only on the *other* sites. Cache that context.
            context = config[:target_index] + config[target_index + 1:]
            key = (target_index, context)
            if key in activity:
                continue
            completed = graph.copy()
            for assigned_site, state_index in zip(sites, config, strict=True):
                if assigned_site.key[1] not in completed.atom_stereo:
                    raise ValueError("Input site unexpectedly disappeared")
                completed.set_atom_stereo(assigned_site.states[state_index])
            activity[key] = locally_stereogenic_at_assignment(
                completed, site.key
            )

    def is_active(target_index: int, config: tuple[int, ...]) -> bool:
        return activity[
            (target_index, config[:target_index] + config[target_index + 1:])
        ]

    edges: set[tuple[int, int]] = set()
    witnesses: dict[str, dict[str, object]] = {}
    status: dict[int, str] = {}

    for target_index, target in enumerate(sites):
        atom = int(target.key[1])
        outcomes = [
            value for (i, _), value in activity.items()
            if i == target_index
        ]
        if all(outcomes):
            status[atom] = "always active"
        elif any(outcomes):
            status[atom] = "conditionally active"
        else:
            status[atom] = "inactive in all tested contexts"

        for source_index, source in enumerate(sites):
            if source_index == target_index or domain_sizes[source_index] < 2:
                continue
            for config in configs:
                for new_state in range(domain_sizes[source_index]):
                    if new_state == config[source_index]:
                        continue
                    altered = list(config)
                    altered[source_index] = new_state
                    altered = tuple(altered)
                    before = is_active(target_index, config)
                    after = is_active(target_index, altered)
                    if before != after:
                        t = int(source.key[1])
                        edges.add((t, atom))
                        witnesses[f"{t}->{atom}"] = {
                            "before": list(config),
                            "after": list(altered),
                            "active_before": before,
                            "active_after": after,
                        }
                        break
                if (int(source.key[1]), atom) in edges:
                    break

    return ExactDependencies(
        edges=frozenset(edges),
        activity=status,
        evaluations=len(activity),
        witnesses=witnesses,
    )


def run_enumeration(
    graph: StereoMolGraph,
    *,
    max_sites: int = 8,
    max_assignments: int = 256,
    max_output: int = 128,
) -> dict[str, object]:
    """Enumerate states of existing tetrahedral descriptors, no new sites."""
    analysis = detect_stereo_elements(graph)
    if any(
        s.key[0] != "atom" or not isinstance(s.template, Tetrahedral)
        for s in analysis.sites
    ):
        raise ValueError("Only tetrahedral atom descriptors are allowed")
    if len(analysis.sites) > max_sites:
        raise ValueError("Too many candidate sites for this baseline")
    theoretical = prod(len(site.states) for site in analysis.sites)
    if theoretical > max_assignments:
        raise ValueError(
            f"Formal space of {theoretical} assignments exceeds "
            f"the limit {max_assignments}"
        )

    entries = []
    count = 0
    truncated = False
    # Count all unique outputs if the formal product fits within the cap.
    for result in enumerate_stereo_assignments(graph, analysis, unique=True):
        count += 1
        if count <= max_output:
            entries.append({
                "index": count,
                "configuration": {
                    str(a): stereo.parity
                    for a, stereo in sorted(result.atom_stereo.items())
                },
            })
        else:
            truncated = True

    return {
        "candidate_sites": len(analysis.sites),
        "formal_assignments": theoretical,
        "unique_stereo_graphs": count,
        "output_truncated": truncated,
        "configurations": entries,
        "site_rows": describe_sites(analysis),
    }


def rdkit_stereoisomer_count(
    molecule: Chem.Mol,
    *,
    max_sites: int = 8,
) -> int | None:
    """Reference count, NOT a ground-truth correctness certificate."""
    candidates = Chem.FindPotentialStereo(molecule)
    if len(candidates) > max_sites:
        return None
    options = StereoEnumerationOptions(
        onlyUnassigned=True,
        unique=True,
        tryEmbedding=False,
        maxIsomers=0,
    )
    return sum(1 for _ in EnumerateStereoisomers(molecule, options=options))


def molecule_svg(bundle: InputBundle, *, width: int = 750, height: int = 380) -> str:
    """RDKit depiction of the molecular skeleton with candidate sites marked."""
    mol = Chem.Mol(bundle.molecule)
    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.drawOptions().addAtomIndices = True
    drawer.DrawMolecule(
        mol, highlightAtoms=sorted(
            a for a in bundle.retained_tetrahedral_atoms
            if a < mol.GetNumAtoms()
        )
    )
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def dependency_svg(
    analysis: StereoAnalysis,
    *,
    exact: ExactDependencies | None = None,
) -> str:
    """Render a directed dependency graph as a standalone SVG.

    Dashed gray = heuristic *positional influence*, not proved dependence.
    Solid = exact conditional stereogenicity influence within the tested
    complete local configuration space.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import networkx as nx
    from io import StringIO

    sites = [site for site in analysis.sites if isinstance(site.template, Tetrahedral)]
    graph = nx.DiGraph()
    for site in sites:
        graph.add_node(int(site.key[1]))
    estimated = set(heuristic_edges(analysis))
    confirmed = set(exact.edges) if exact is not None else set()
    graph.add_edges_from(estimated | confirmed)

    fig, ax = plt.subplots(figsize=(9, 6))
    if graph.number_of_nodes() > 0:
        # Deterministic layout across runs and atom order.
        layout = nx.spring_layout(graph, seed=17)
        colors = []
        for node in graph.nodes:
            state = exact.activity.get(node) if exact is not None else None
            if state == "conditionally active":
                colors.append("#F5B050")
            elif state == "always active":
                colors.append("#7BBDA3")
            elif state == "inactive in all tested contexts":
                colors.append("#CED3DA")
            else:
                colors.append("#98BCE5")
        nx.draw_networkx_nodes(
            graph, layout, node_color=colors, node_size=1550,
            edgecolors="#35404b", linewidths=1.2, ax=ax,
        )
        nx.draw_networkx_labels(
            graph, layout, labels={n: f"C{n}" for n in graph.nodes},
            font_size=11, font_weight="bold", ax=ax,
        )
        possible = sorted(estimated - confirmed)
        if possible:
            nx.draw_networkx_edges(
                graph, layout, edgelist=possible, ax=ax,
                edge_color="#9099a4", style="dashed", arrows=True,
                arrowsize=22, node_size=1550,
                connectionstyle="arc3,rad=0.12", width=1.8,
            )
        if confirmed:
            nx.draw_networkx_edges(
                graph, layout, edgelist=sorted(confirmed), ax=ax,
                edge_color="#287d71", arrows=True, arrowsize=24,
                node_size=1550, connectionstyle="arc3,rad=0.12",
                width=2.7,
            )
    else:
        ax.text(0.5, 0.5, "No candidate tetrahedral sites",
                ha="center", va="center", transform=ax.transAxes)

    ax.set_title(
        "Tetrahedral-site dependency network\n"
        "dashed: possible positional influence; "
        "solid: exact activity change",
        fontsize=12,
    )
    ax.axis("off")
    fig.tight_layout()
    buffer = StringIO()
    fig.savefig(buffer, format="svg", bbox_inches="tight")
    plt.close(fig)
    return buffer.getvalue()


def write_detection_examples(
    output_dir: Path,
    *,
    exact: bool = False,
    matched_rdkit_sites: bool = True,
) -> list[dict[str, object]]:
    """Generate reproducible example molecule depictions and graphs."""
    from .examples import EXAMPLES

    output_dir.mkdir(parents=True, exist_ok=True)
    reports = []
    for name in (
        "VS006 — bridged cage",
        "VS294 — bridged oxygen rings",
        "VS295 — related oxygen rings",
    ):
        sample_id = name.split(" ")[0]
        bundle = parse_smiles(
            EXAMPLES[name], matched_rdkit_sites=matched_rdkit_sites
        )
        analysis = detect_stereo_elements(bundle.graph)
        exact_result = None
        if exact:
            exact_result = exact_dependency_analysis(bundle.graph, analysis)
        (output_dir / f"{sample_id}_molecule.svg").write_text(
            molecule_svg(bundle), encoding="utf-8"
        )
        (output_dir / f"{sample_id}_dependencies.svg").write_text(
            dependency_svg(analysis, exact=exact_result), encoding="utf-8"
        )
        reports.append({
            "name": name,
            "smiles": EXAMPLES[name],
            "candidate_sites": list(sorted(bundle.retained_tetrahedral_atoms)),
            "classification": describe_sites(analysis),
            "heuristic_edges": [list(e) for e in sorted(heuristic_edges(analysis))],
            "exact_edges": (
                [list(e) for e in sorted(exact_result.edges)]
                if exact_result is not None else None
            ),
            "exact_site_activity": exact_result.activity if exact_result else None,
            "exact_witnesses": exact_result.witnesses if exact_result else None,
        })
    import json
    (output_dir / "examples.json").write_text(
        json.dumps(reports, indent=2), encoding="utf-8"
    )
    return reports
