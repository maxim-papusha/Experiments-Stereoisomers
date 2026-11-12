"""Application 2: tetrahedral stereoisomer assignment generation.

Run with: streamlit run apps/stereoisomer_generation.py
"""

import json
import sys
from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from stereoapps.core import (
    molecule_svg,
    parse_smiles,
    rdkit_stereoisomer_count,
    run_enumeration,
)
from stereoapps.examples import EXAMPLES


st.set_page_config(page_title="Stereoisomer generation", layout="wide")
st.title("Application 2 — Tetrahedral stereoisomer generation")
st.caption(
    "Exhaustive assignment of existing local tetrahedral states, "
    "followed by StereoMolGraph stereochemical graph-isomorphism deduplication."
)

with st.form("generation"):
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
        unassign = st.checkbox("Enumerate from unassigned stereo", value=True)
    with c2:
        matched = st.checkbox(
            "Use RDKit-matched tetrahedral sites",
            value=True,
            help=(
                "Controls the starting site set for a fair count comparison. "
                "Only deletes extra imported descriptors; never adds a site."
            ),
        )
    with c3:
        max_output = st.number_input(
            "Maximum configurations displayed", value=64,
            min_value=1, max_value=256, step=1,
        )
    submitted = st.form_submit_button("Generate configurations", type="primary")

if submitted:
    try:
        bundle = parse_smiles(
            smiles, unassign=unassign, matched_rdkit_sites=matched
        )
        if bundle.missing_rdkit_candidates:
            st.warning(
                "RDKit potential tetrahedral sites absent in the supplied "
                "StereoMolGraph descriptors: "
                + str(sorted(bundle.missing_rdkit_candidates))
                + ". No new descriptors are created."
            )

        with st.spinner("Enumerating stereochemical assignments"):
            result = run_enumeration(
                bundle.graph, max_sites=8, max_assignments=256,
                max_output=int(max_output),
            )
            reference = rdkit_stereoisomer_count(bundle.molecule)

        st.success("Enumeration completed")
        m1, m2, m3 = st.columns(3)
        m1.metric("Formal local assignments", result["formal_assignments"])
        m2.metric("Unique StereoMolGraph outputs", result["unique_stereo_graphs"])
        m3.metric(
            "RDKit outputs (reference only)",
            reference if reference is not None else "limit exceeded",
        )

        if reference is not None and reference != result["unique_stereo_graphs"]:
            st.warning(
                "The engines disagree. This does not by itself identify "
                "which result is chemically correct; differences can arise "
                "from cyclic geometry and inactive local descriptors."
            )

        left, right = st.columns([1, 1])
        with left:
            st.subheader("Input molecule")
            components.html(molecule_svg(bundle), height=420)
        with right:
            st.subheader("Stereoisomer assignment table")
            rows = [
                {"#": x["index"], **x["configuration"]}
                for x in result["configurations"]
            ]
            st.dataframe(rows, hide_index=True, use_container_width=True)
            if result["output_truncated"]:
                st.caption(
                    "The displayed table is truncated, but all unique "
                    "outputs within the formal assignment cap were counted."
                )

        with st.expander("Candidate tetrahedral site details"):
            st.dataframe(
                result["site_rows"], use_container_width=True, hide_index=True
            )

        report = {
            "input_smiles": smiles,
            "example": example,
            "rdkit_comparison": reference,
            "matched_rdkit_sites": matched,
            **result,
        }
        st.download_button(
            "Download enumeration results (JSON)",
            data=json.dumps(report, indent=2),
            file_name="stereoisomers.json",
            mime="application/json",
        )
    except Exception as exc:
        st.error(f"{type(exc).__name__}: {exc}")

st.info(
    "The output is a set of distinct StereoMolGraph configurations "
    "under the current graph representation, not a guarantee of "
    "physically stable stereoisomers. Automorphism pruning and "
    "inactive-site normalization are not implemented."
)
