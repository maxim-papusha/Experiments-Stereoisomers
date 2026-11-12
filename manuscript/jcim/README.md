# JCIM manuscript — working outline

This folder contains an **initial research-article outline**, not a submission-ready manuscript.

- `main.tex`: ACS `achemso` class (`journal=jcisd8`, `manuscript=article`), proposed title, preliminary abstract, equations, drafted methodology, existing pilot findings, figure placeholders, and explicitly marked open work.
- `references.bib`: verified starting references, including the published 2026 StereoMolGraph paper, Rieder *et al.* on *enu*, McKay's isomorph-free generation, and RDKit documentation.

Build with a TeX Live installation containing `achemso` and BibTeX:

```bash
cd manuscript/jcim
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

## Scope and planned claims

The **input** is a StereoMolGraph with all potential tetrahedral sites already represented as descriptors. Unassigned configurations have `parity=None`; detection and enumeration **do not add descriptor sites**. Initially ignore isotope labels, multiple bonds, and non-tetrahedral stereochemistry.

The paper will investigate (1) exact site-anchored stereogenicity, (2) configuration-dependent stereo activation networks, (3) exhaustive and symmetry-orbit-validated stereoisomer generation, and (4) eventual group-theoretical pruning.

The detailed VS294/VS295 case has preliminary evidence for seven configured graph classes versus six default RDKit outputs; the manuscript carefully separates this combinatorial result from experimental isolability.

## Before submission

Confirm authors, affiliation, and correspondence details. Replace all **To develop** markers with validated material. Run a larger curated benchmark, compare complete isomer sets, add final chemically checked molecular figures and a graphical TOC entry, provide Supporting Information as a separate file, pin computational versions, and verify the latest JCIM author guidelines.

The `achemso` class is provided as requested. The package's current documentation describes its class as legacy-supported; confirm the publication-time ACS template requirements before submitting.
