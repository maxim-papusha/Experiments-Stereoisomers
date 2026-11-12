# Scientific and editorial review — JCIM outline

**Review date:** 13 November 2026  
**Manuscript:** [achemso main.tex](jcim/main.tex)

## Changes made

1. Rewrote the abstract and Introduction to start with accessible tetrahedral examples rather than group theory.
2. Added a worked conceptual illustration using **2,3-dichlorobutane**: four local configurations but three molecular stereoisomers because of a meso equivalence.
3. Defined chemical and mathematical terms in the order in which they are used:
   local tetrahedral parity, allowed local states, complete assignments,
   site-anchored stereo graph isomorphism, activation function,
   directional dependency witnesses, atom orbits, and global stereo graph classes.
4. Distinguished graph-theoretical **activity** from standard IUPAC
   stereogenicity and from experimentally isolable stereoisomers.
5. Explicitly separated the local test fixing the selected atom from the
   global whole-molecule equivalence used to eliminate duplicated outputs.
6. Clarified that a dependency arrow only asserts the **existence** of
   an assignment pair that changes the target's activity. A directed
   graph alone cannot encode higher-order joint-configuration effects.
7. Explained automorphisms and Burnside's counting relation in words,
   with the parity-transport caveat. Described McKay pruning as
   future work rather than an implemented result.
8. Tightened the VS294 text and stated its independently validated
   pilot count (seven tetrahedral graph classes vs six default
   RDKit outputs) without asserting an experimentally established
   chemical stereoisomer count.
9. Rewrote the Methods to explain chemical scope, controls,
   independent validation, exact-vs-heuristic distinctions and
   physical feasibility tests.

## Literature review and bibliography

The BibTeX file now includes **15 references**, spanning:

- **IUPAC definitions:** Moss (1996), Favre and Powell (2014).
- **CIP software and validation suite:** Hanson et al. (2018).
- **Prior stereo representations and enumeration:** Sobez and
  Reiher (2020), Rieder et al. (2023), Papusha and Leonhard (2026).
- **Graph algorithms:** Morgan (1965), Cordella et al. (2004),
  McKay (1998), NetworkX/Hagberg et al. (2008).
- **RDKit:** documentation and version-specific software release (2026.09.2).
- **Geometry:** Riniker and Landrum (2015), Wang et al. (2020),
  Halgren (1996).

Publication metadata were checked with publisher and other
authoritative records where available. All cite keys used in
the LaTeX manuscript resolve in references.bib. DOI coverage
is broad but not every software/proceedings record has one.
These are substantive starting references, not a guarantee
that every specialized prior-art study has been identified.

## Scientific points to retain in later revisions

- The input is already a **StereoMolGraph** with existing
  tetrahedral descriptors. Undefined parity is not a wildcard
  for equality and does not authorize adding candidate sites.
- The two parity signs are relative to local ligand ordering,
  **not** CIP \(R/S\) labels.
- The local activity function fixes the target atom during
  stereo graph isomorphism; completed-molecule enumeration does not.
- A positional color-refinement edge is a heuristic,
  not a certified configuration dependency.
- Global symmetry orbit counts require properly transporting
  tetrahedral parity under graph automorphisms.
- Graph distinctness does not alone establish geometric
  realizability, configurational stability or stereogenicity
  according to all IUPAC chemical conventions.
- The VS294/VS295 pilot is especially strong because
  independent group theory and MMFF 3D probes complement
  StereoMolGraph graph comparisons, but its scope remains limited.

## Important work still required before JCIM submission

- Add verified, carefully drawn molecular examples
  for **orbit-preserving** versus **orbit-splitting**
  stereochemical decoration. Do not re-use an unchecked
  diagram with unreliable stereochemical configurations.
- Add a molecular figure and activation truth table
  for VS294 with exact assignment witnesses.
- Curate a larger comparison dataset with independent
  equivalence certificates, not only RDKit output counts.
- Run complete atom-relabeling and ligand-slot permutation
  invariance checks, and characterize errors on cyclic molecules.
- Report reproducible performance data; do not claim
  automorphism-pruned speedups before implementing them.
- Check and describe all physical constraints for the
  proposed structures, where needed.
- Add a validated author affiliation and correspondence
  before journal submission, as well as a TOC graphic,
  graphical figures, separate Supporting Information,
  archival software versions and a data release.

The author is listed as **Maxim Papusha**, without
affiliation, at the user's request.
