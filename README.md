# Experiments — Stereoisomers

Two **separate, locally runnable Streamlit applications** built on
[the StereoMolGraph stereoisomer branch](https://github.com/maxim-papusha/StereoMolGraph/tree/stereoisomer):

## Application 1: Stereocenter detection and dependency graph

Run:

    streamlit run apps/stereocenter_detection.py

- Read tetrahedral candidates from an already prepared StereoMolGraph.
- Optionally build a StereoMolGraph from SMILES using the existing RDKit import adapter.
- Classify existing candidate sites using fast color refinement.
- Show the molecular skeleton with candidate atoms highlighted.
- Draw an **oriented dependency graph** among candidate tetrahedral sites.
- Optionally perform **exact, context-dependent activity analysis** by
  enumerating other local assignments and testing site-anchored stereo
  isomorphism (at most 5 sites/32 formal assignments).
- Export the report as JSON and the network as SVG.

**Dashed gray arrows** represent heuristic positional influence, not proven
stereochemical dependence. If exact analysis is selected, **solid green
arrows** indicate dependencies certified *under the present StereoMolGraph
graph-isomorphism model*. Neither graph proves physical ring feasibility.

## Application 2: Stereoisomer generation

Run:

    streamlit run apps/stereoisomer_generation.py

- Use *only the stereo descriptor sites already present*.
- Enumerate the local tetrahedral states supplied by
  StereoMolGraph descriptors; specified parities stay fixed.
- Remove duplicate completed graphs using the existing stereo-aware
  equality/isomorphism implementation.
- Compare the output count against RDKit as a reference; no claim of
  independent validation or equality of complete isomer **sets**.
- Export the configuration table and counts as JSON.

This is a complete baseline for **represented local descriptor assignments**,
not a complete stereochemical feasibility or stereocenter normalization
algorithm. No McKay-style automorphism pruning is used yet.

## Installation

Python 3.11+ and Git are required.

    git clone https://github.com/maxim-papusha/Experiments-Stereoisomers.git
    cd Experiments-Stereoisomers
    pip install -e ".[test]"

The installer retrieves StereoMolGraph from its experimental
\`stereoisomer\` branch. For reproducible publications pin the exact
StereoMolGraph commit and RDKit version.

## Scope

- Current experiments focus on **tetrahedral atom configurations**.
- Isotope-labelled inputs and inputs containing aromatic or multiple
  bonds are excluded (a stricter restriction than E/Z removal).
- Geometry-nomination is performed *before* the algorithm; neither
  application creates any stereo descriptor site.
- If a SMILES is used, the upstream StereoMolGraph importer may nominate
  many geometry-eligible but non-stereogenic candidates. For direct
  RDKit comparisons an optional checkbox filters those **existing**
  candidates to the RDKit potential tetrahedral site list. It never
  inserts missing sites.
- Molecular models and stereo-activity tests are graph-theoretical.
  3D ring realizability, interconversion barriers, and inactive-site
  normalization are not yet covered.

## Generated complex-molecule dependency examples

The following structures originate from the
[canonicalization experiment's 300 CIP Validation Suite SMILES](https://github.com/maxim-papusha/Experiments-CanonicalStereoMolGraph/blob/main/data/smiles.txt):

- **VS006:** bridged tetrahedral cage.
- **VS294:** cyclic oxygen-rich molecule with five candidate centers.
- **VS295:** related cyclic molecule with five candidate centers.

Generate both molecular drawings and dependency networks:

    python examples/generate_dependency_figures.py --output examples/figures

For an optional expensive exact dependency analysis:

    python examples/generate_dependency_figures.py --output examples/figures --exact

Expected files include
\`VS294_molecule.svg\`, \`VS294_dependencies.svg\`, and
\`examples.json\`. CI uploads these as artifacts, and the repository
keeps a static set of generated example graphs for inspection.

## Tests

    pytest -q

The GitHub Actions workflow independently tests the applications and
regenerates complex-molecule dependency figures. The applications can
be run separately and share only the underlying science/model layer.
