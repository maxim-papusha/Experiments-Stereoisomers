"""Selected isotope-free, single-bond examples from the CIP Validation Suite.

Source: maxim-papusha/Experiments-CanonicalStereoMolGraph/data/smiles.txt
The examples are structural test inputs, not independently validated
enumeration reference counts.
"""

EXAMPLES = {
    "Simple acyclic": "C[C@H](O)N",
    "VS006 — bridged cage": "C1[C@H]2C[C@H]3C[C@@H]1C[C@@H](C2)C3",
    "VS280 — symmetric cage": "C1C[C@]23CCC[C@@]24CCC[C@@]45CCC[C@@]56CCC[C@]36C1",
    "VS294 — bridged oxygen rings": "O[C@H]1O[C@@]2(O[C@H]1O)O[C@H]([C@H](O2)O)O",
    "VS295 — related oxygen rings": "O[C@H]1O[C@@]2(O[C@H]1O)O[C@@H]([C@@H](O2)O)O",
    "VS298 — bicyclic dichloride": "Cl[C@H]1C[C@]2(C1)C[C@H](C2)Cl",
}
