# VS294 and VS295: seven stereoisomer classes or six?

**Conclusion (under the molecular model tested): seven distinct tetrahedral
stereochemical graph classes.** RDKit 2026.09.2's default stereochemical
cleanup merges two graph-theoretically distinct classes by discarding the
configuration at atom **C3**. Direct 3D embedding and MMFF minimization
also produce valid, stereochemically distinct local geometries for this
pair, although this does not establish experimental isolability.

- SMILES of VS294:
  \`O[C@H]1O[C@@]2(O[C@H]1O)O[C@H]([C@H](O2)O)O\`
- VS295 has the same connectivity, with a different initial set of
  specified tetrahedral tags. Full de novo enumeration has the same
  candidate assignment space.
- [Reproducible independent analysis](verify_vs294_seven.py)
- [Output JSON from GitHub Actions](results/vs294_analysis.json)
- [Validated test run](https://github.com/maxim-papusha/Experiments-Stereoisomers/actions/runs/43322221338)

## 1. What was independently tested?

Three separate checks were applied:

1. **SMG stereo-isomorphism:** enumerate all 32 assignments from
   pre-existing tetrahedral descriptors in StereoMolGraph. Using
   \`unique=True\`, obtain exactly **7** graphs. Compare all
   \(\binom{7}{2}=21\) pairs again using the direct
   \`vf2pp_all_isomorphisms(stereo=True)\` implementation, independently
   of Python set/hash-based deduplication: **zero** isomorphic pairs.
2. **Independent group theory:** construct an explicit-hydrogen graph
   using RDKit (atomic numbers and connectivity), calculate *all*
   constitutional graph automorphisms using NetworkX, and implement
   parity transport on the local tetrahedral ligand orderings without
   calling the StereoMolGraph equality algorithm.
3. **RDKit cleanup and physical-geometry probe:** force every local
   tetrahedral tag, inspect how \`Chem.AssignStereochemistry(cleanIt=True)\`
   modifies them, and embed/minimize the configurations differing only
   at central site C3 using ETKDGv3 and MMFF.

The analysis was run against **RDKit 2026.09.2** on GitHub Actions,
and the independent RDKit cleanup test was also reproduced locally
with RDKit **2025.09.4**.

## 2. Independent proof of seven formal classes

The potential tetrahedral centers in the original RDKit atom indexing are

$$
S=\{1,3,5,8,9\}.
$$

Each has a local state space \(\{-1,+1\}\), so the unrestricted
configuration space consists of

$$
|\Omega|=2^5=32
$$

assignments.

The hydrogen-explicit constitutional graph has **8 automorphisms**.
Each automorphism \(\phi\) induces a permutation of the five tetrahedral
sites, together with a handedness transformation at each site.

Let \(L_s=(\ell_{s,1},\ldots,\ell_{s,4})\) be a fixed ordered
reference list of neighbors of site \(s\).
The mapped neighbor list
\(\phi(L_s)\) is a permutation \(\pi_{\phi,s}\) of
\(L_{\phi(s)}\). The correct transported parity is

$$
(\phi\cdot \sigma)_{\phi(s)}
=
\operatorname{sgn}(\pi_{\phi,s})\,\sigma_s.
\tag{1}
$$

The two assignments are symmetry equivalent if and only if they
are in the same orbit under this action.

Applying all 8 automorphisms gives the following **7 orbit sizes**:

$$
2,\ 2,\ 4,\ 4,\ 4,\ 8,\ 8.
\tag{2}
$$

They add up to 32 and are mutually disjoint. Burnside's lemma provides
an additional check:

$$
|\Omega/\operatorname{Aut}(G)|
=
\frac{1}{8}\sum_{\phi\in\operatorname{Aut}(G)}
|\operatorname{Fix}(\phi)|
=
\frac{32+8+8+8}{8}
=
\boxed{7}.
\tag{3}
$$

The independent algorithm returns **seven**, and the seven actual
StereoMolGraph output representatives belong to **all seven different
independent orbits**: one representative from each orbit, no duplicates.

This checks the graph-isomorphism result without trusting either
the SMG graph hash or its specialized VF2++ stereochemical feasibility
predicate.

## 3. The precise RDKit discrepancy

RDKit's \`FindPotentialStereo\` reports *all five* tetrahedral
candidates, including C3. RDKit's default enumeration then gives
**six** unique canonical isomeric SMILES.

When we explicitly set all \(2^5=32\) stereotags and call

~~~python
Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
~~~

the C3 tag becomes \`CHI_UNSPECIFIED\` **in all 32 assignments**,
even when the four remote configurations make C3 stereogenic.

Consequently, two distinct full stereo orbits (independent orbit
indices **0** and **3**) are assigned the **same** cleaned canonical SMILES:

~~~text
O[C@@H]1OC2(O[C@@H](O)[C@@H](O)O2)O[C@@H]1O
~~~

Those two orbits each contain 4 of the 32 formal assignments.
All other five orbits map one-to-one to their RDKit-cleaned
SMILES classes. The specific difference is therefore **7 → 6**,
not an accidental discrepancy somewhere else in the enumeration.

Furthermore, inversion of *all* local tetrahedral configurations
maps orbit 0 onto orbit 3: they are an **enantiomeric pair** in the
graph-stereochemical model. RDKit's cleanup thus loses this
pairwise distinction.

### Minimal explicit stereo-pair example

Holding the configurations at C1, C5, C8 and C9 fixed, opposite C3
stereotags produce the following strings before RDKit cleanup:

~~~text
O[C@@H]1O[C@@]2(O[C@@H](O)[C@@H](O)O2)O[C@@H]1O
O[C@@H]1O[C@]2(O[C@@H](O)[C@@H](O)O2)O[C@@H]1O
~~~

After \`AssignStereochemistry(cleanIt=True)\`, both lose the
\`[C@]\`/\`[C@@]\` tag at C3 and become the *same* isomeric SMILES.

The problem is **not** a limitation of detecting the site:
RDKit's \`FindPotentialStereo\` already identifies C3.
It is a discrepancy in the later cleanup/serialization and
enumeration of the conditionally stereogenic central atom.

The official [RDKit Book](https://www.rdkit.org/docs/RDKit_Book.html)
describes both the legacy \`AssignStereochemistry\` route and newer
\`FindPotentialStereo\` algorithm, and notes improved handling of
parastereochemistry in the latter. The experiment demonstrates
that the two stages do not yield the same distinction here.

## 4. Are the extra configurations geometrically possible?

To test whether the pair is merely an impossible formal ring assignment:

- Preserve the existing C1, C5, C8, C9 tetrahedral tags.
- Set C3 to one or the opposite tetrahedral tag.
- Call \`Chem.AssignStereochemistry(cleanIt=False)\` to keep it.
- Add explicit H atoms.
- Embed **10** conformers of each with RDKit ETKDGv3.
- Minimize all of them using MMFF.
- Compute the signed volume of the four neighboring atom
  positions at each center:

$$
V_s=
\det\left(
\mathbf r_1-\mathbf r_4,\,
\mathbf r_2-\mathbf r_4,\,
\mathbf r_3-\mathbf r_4
\right).
\tag{4}
$$

**Results:** all 20 embeddings and optimizations succeeded.
Every structure retained nonzero signed volumes at the five centers.
The four unchanged centers retained the same orientation in both
configurations, and the signed C3 volume changed sign. The optimization
also yielded finite MMFF local minima on both sides.

The numerical minimum energies from this run were approximately
\(-84.48\) and \(-84.60\) kcal/mol, respectively. These are force-field
energies for the explored conformers, not thermodynamic free energies,
and should not be interpreted as a quantitative stereoisomer
stability comparison.

**Interpretation:** the pair is realizable as distinct tetrahedral
3D assignments within the force-field model; the extra configuration
is not eliminated by ordinary ring-closure geometry constraints.

This is strong evidence for **RDKit undercounting the configured
stereochemical states** for this molecule.

## 5. Limits of the conclusion

The evidence establishes:

- seven pairwise non-isomorphic **stereochemically labelled graphs**;
- seven classes from independent graph automorphism action;
- one pair of these classes erased during RDKit's cleanup;
- successful 3D geometries for both opposite C3 configurations.

It does **not** establish an experimentally measured number of
isolable compounds, a kinetic barrier to interconversion, or the
correct IUPAC stereodescriptor names. To make that stronger claim,
we would need to assess whether the tetrahedral configurations
can interconvert by an allowable chemical process, perform suitable
energy/barrier calculations if relevant, and/or consult experimental
structure/optical activity data.

For software development, however, **seven is the correct reference
count under the stated tetrahedral stereo-graph model**.
The unpruned StereoMolGraph enumerator should retain this case as
a regression benchmark, and an automated cleanup should **not**
eliminate C3 solely because its ligands are equivalent in
the unassigned constitutional graph.

## 6. Reproduce

From the root of
\`maxim-papusha/Experiments-Stereoisomers\`:

~~~bash
python investigations/verify_vs294_seven.py \
    --conformers 10 \
    --output investigations/results/vs294_analysis.json
~~~

The independent script asserts the mathematical orbit count,
all 21 direct pairwise VF2 comparisons, RDKit's 6-count cleanup
and the 3D signed-volume checks. The GitHub Actions workflow
executes the analysis and commits the resulting JSON report.
