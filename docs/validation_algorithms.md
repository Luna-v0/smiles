# Simple Molecular Validation Algorithms

The SMILES Validator uses two fundamental chemistry rules to validate molecular structures: the **Octet Rule** for valency checking and **Huckel's Rule** for aromaticity validation. These algorithms operate on the molecular graph built during SMILES parsing.

## Octet Rule (Valency Validation)

The octet rule states that atoms tend to form bonds until they are surrounded by **8 valence electrons**, achieving the electron configuration of a noble gas. Lighter elements (hydrogen, helium) follow the **duet rule**, requiring only **2 electrons**.

### How It Works

The validator checks valency for bracket atoms (`[...]` notation in SMILES) that are **not aromatic** and **not part of a ring**. The process is:

1. **Parse the electron configuration** from the periodic table data (e.g., Carbon: `1s2 2s2 2p2` = 4 valence electrons).
2. **Adjust for charge and hydrogens** specified in the bracket atom:
    - Negative charge adds electrons (e.g., `[O-]` adds 1 electron).
    - Positive charge removes electrons (e.g., `[Li+]` removes 1 electron).
    - Each explicit hydrogen adds 1 bonding electron (e.g., `[NH3]` adds 3 electrons).
3. **Count bond electrons from the molecular graph** — each bond to a neighbor contributes shared electrons (single = 1, double = 2, triple = 3).
4. **Check the effective electron count** against the target:
    - **H, He**: exactly 2 electrons (duet rule).
    - **Li, Be**: must reach 2 electrons through internal configuration only (charge), not covalent bonds.
    - **All other elements**: at least 8 electrons (octet rule, with expanded octets allowed).

### Validation Scope

| Atom Type | Valency Check |
|---|---|
| Regular atoms (not bracket) | Skipped (assumed valid) |
| Aromatic bracket atoms | Skipped (delocalized bonding) |
| Ring bracket atoms | Skipped (ring bonds contribute) |
| Acyclic bracket atoms | Full octet/duet check |

### Examples

| SMILES | Atom | Valence e- | Adjustment | Graph Bonds | Effective | Target | Valid |
|---|---|---|---|---|---|---|---|
| `[H]c1ccccc1` | `[H]` | 1 | none | 1 (to c) | 2 | 2 | Yes |
| `[NH3]` | `[NH3]` | 5 | +3 (H) | 0 | 8 | 8 | Yes |
| `[Li+]` | `[Li+]` | 1 | -1 (charge) | 0 | 2 | 2 | Yes |
| `CC[Li]` | `[Li]` | 1 | none | 1 (to C) | 2 (internal=1) | 2 | No (Li requires charge) |
| `[O-]c1ccccc1` | `[O-]` | 6 | +1 (charge) | 1 (to c) | 8 | 8 | Yes |

## Huckel's Rule (Aromaticity Validation)

Huckel's rule states that a planar, cyclic, fully conjugated molecule is aromatic if it contains **4n + 2 pi electrons**, where n is a non-negative integer (n = 0, 1, 2, ...). This gives valid pi electron counts of **2, 6, 10, 14, 18, ...** for aromatic systems.

### How It Works

After parsing the SMILES and building the molecular graph, the validator:

1. **Identifies all cycles** in the molecular graph.
2. **Groups cycles into fused ring systems** — two cycles are fused if they share 2 or more atoms.
3. **For isolated aromatic cycles**, validates each independently:
    - Count pi electrons based on atom type.
    - Check if the count satisfies 4n + 2.
4. **For fused aromatic systems**, uses relaxed validation:
    - If at least one individual cycle satisfies Huckel's rule, the system is valid.
    - For complex systems (3+ rings where all atoms are aromatic), the SMILES notation is trusted.

### Pi Electron Contributions

Each aromatic atom contributes pi electrons based on its type and bonding:

| Atom | Condition | Pi Electrons | Example |
|---|---|---|---|
| Carbon (c) | Aromatic | 1 | Benzene ring carbon |
| Nitrogen (n) | No explicit H (pyridine-like) | 1 | Pyridine nitrogen |
| Nitrogen ([nH]) | Has explicit H (pyrrole-like) | 2 | Pyrrole nitrogen |
| Oxygen (o) | Aromatic | 2 | Furan oxygen |
| Sulfur (s) | Aromatic | 2 | Thiophene sulfur |
| Boron (b) | Aromatic | 1 | Borole boron |
| Phosphorus (p) | Aromatic | 1 | Phosphorine phosphorus |

### Examples

**Benzene** (`c1ccccc1`): 6 aromatic carbons, each contributing 1 pi electron = **6 pi electrons**. Check: 4(1) + 2 = 6. Valid (n=1).

**Pyrrole** (`c1cc[nH]c1`): 4 aromatic carbons (4 pi e-) + 1 pyrrole nitrogen with H (2 pi e-) = **6 pi electrons**. Valid (n=1).

**Cyclobutadiene** (hypothetical `c1ccc1`): 4 aromatic carbons = **4 pi electrons**. Check: 4n + 2 = 4 has no integer solution. Invalid (antiaromatic).

**Naphthalene** (`c1ccc2ccccc2c1`): Fused system with 2 rings. Each 6-membered ring individually has 6 pi electrons. At least one satisfies 4n + 2, so the fused system is valid.

### Fused Ring System Handling

Fused aromatic systems (e.g., naphthalene, anthracene) require special treatment because pi electrons are delocalized across the entire system. The validator uses a tiered approach:

1. **Simple fused systems (2 rings)**: Check if at least one ring satisfies Huckel's rule individually.
2. **Complex fused systems (3+ rings, all aromatic)**: Trust the SMILES notation as valid. This handles molecules like perylene where individual cycle detection may not perfectly decompose the ring system.
3. **Mixed systems**: If a fused system contains both aromatic and non-aromatic atoms, it is considered valid (the non-aromatic parts are aliphatic substituents).

### Validation Flow

```
Cycles detected?
  |
  No --> Valid (no aromaticity to check)
  |
  Yes --> Group into fused systems
            |
            For each system:
              |
              Single cycle? --> Count pi electrons --> Check 4n+2
              |
              Fused system? --> Check individual cycles
                                  |
                                  Any cycle satisfies 4n+2? --> Valid
                                  |
                                  No? --> Check total pi electrons
```
