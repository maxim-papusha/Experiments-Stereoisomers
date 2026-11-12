"""Application 1: tetrahedral stereocenter detection and dependencies.

Run with: streamlit run apps/stereocenter_detection.py
"""

import json
import sys
from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from stereoapps.core import (
    dependency_svg,
    describe_sites,
    detect_stereo_elements,
    exact_dependency_analysis,
    heuristic_edges,
    molecule_svg,
    parse_smiles,
)
from stereoapps.examples import EXAMPLES


st.set_page_config(page_title="Stereocenter detection", layout="wide")
st.title("Application 1 — Tetrahedral stereocenter detection")
st.caption(
    "Investigate potential stereocenters and stereochemical dependencies. "
    "No stereodescriptors are introduced during detection."
)

with st.form("detection"):
    example = st.selectbox(
        "Example molecule", list(EXAMPLES),
        index=list(EXAMPLES).index("VS294 — bridged oxygen rings"),
    )
    smiles = st.text_area(
        "SMILES (editable)", EXAMPLES[example],
        key="smiles_" + example, height=100,
    )
    c1, c2, c3 = st.columns(3)
    with c1:
        unassign = st.checkbox(
            "Remove specified stereo for discovery", value=True
        )
    with c2:
        matched = st.checkbox(
            "Use RDKit's tetrahedral candidate-site list",
            value=False,
            help=(
                "For controlled comparisons only. Filters existing "
                "StereoMolGraph descriptors, never adds new ones."
            ),
        )
    with c3:
        exact = st.checkbox(
            "Compute exact contextual dependencies (may be slow)",
            value=False,
            help=(
                "Enumerates all local assignments and tests site-anchored "
                "graph isomorphism. Restricted to 5 sites/32 assignments."
            ),
        )
    submitted = st.form_submit_button(
        "Analyze tetrahedral stereochemistry", type="primary"
    )

if submitted:
    try:
        with st.spinner("Analyzing input stereochemistry"):
            bundle = parse_smiles(
                smiles, unassign=unassign, matched_rdkit_sites=matched
            )
            analysis = detect_stereo_elements(bundle.graph)
            exact_result = exact_dependency_analysis(
                bundle.graph, analysis
            ) if exact else None
            rows = describe_sites(analysis)
            possible = heuristic_edges(analysis)
            network = dependency_svg(analysis, exact=exact_result)
            molecule = molecule_svg(bundle)

        st.success("StereoMolGraph analysis completed")
        m1, m2, m3 = st.columns(3)
        m1.metric("Supplied tetrahedral sites", len(rows))
        m2.metric("Potential influence edges", len(possible))
        m3.metric(
            "Confirmed dependency edges",
            len(exact_result.edges) if exact_result is not None else "not tested",
        )

        if bundle.missing_rdkit_candidates:
            st.warning(
                "RDKit sites absent from the input StereoMolGraph: "
                + str(sorted(bundle.missing_rdkit_candidates))
                + ". Detection never adds missing descriptors."
            )

        c1, c2 = st.columns([1, 1])
        with c1:
            st.subheader("Molecular skeleton")
            components.html(molecule, height=420, scrolling=False)
        with c2:
            st.subheader("Stereochemical dependency graph")
            components.html(network, height=620, scrolling=False)

        st.caption(
            "An arrow t → s means a change at t may influence the "
            "stereogenicity of s. Gray dashed arrows are heuristic "
            "positional influences, NOT confirmed chemical dependencies. "
            "Green solid arrows are confirmed within the represented "
            "StereoMolGraph configuration model if exact analysis was enabled."
        )
        st.subheader("Candidate site classifications")
        st.dataframe(rows, use_container_width=True, hide_index=True)

        if exact_result is not None:
            st.subheader("Contextual stereogenicity")
            st.dataframe(
                [
                    {"atom": atom, "stereogenicity": state}
                    for atom, state in sorted(exact_result.activity.items())
                ],
                hide_index=True,
                use_container_width=True,
            )
            with st.expander("Exact dependency witnesses"):
                st.json(exact_result.witnesses)

        report = {
            "input_smiles": smiles,
            "example": example,
            "supplied_sites": rows,
            "heuristic_directed_edges": sorted(possible),
            "exact_directed_edges": (
                sorted(exact_result.edges) if exact_result else None
            ),
            "exact_activity": exact_result.activity if exact_result else None,
            "site_source": "RDKit-filtered" if matched else "StereoMolGraph",
        }
        st.download_button(
            "Download detection report (JSON)",
            data=json.dumps(report, indent=2),
            file_name="stereocenter_detection.json",
            mime="application/json",
        )
        st.download_button(
            "Download dependency graph (SVG)",
            data=network,
            file_name="stereochemical_dependencies.svg",
            mime="image/svg+xml",
        )
    except Exception as exc:
        st.error(f"{type(exc).__name__}: {exc}")

st.info(
    "A geometry-qualified candidate is not necessarily a stereogenic "
    "center. A site is conditionally stereogenic when its activity "
    "depends on assignments at other sites. Physical ring feasibility "
    "and stereochemical stability are outside the present model."
)
