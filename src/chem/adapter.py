"""Adapter from our ``MolecularGraph`` to a backend chemistry engine.

This module is the integration seam described in
``docs/chemistry_backend_study.md``.  It converts the molecule graph that our
SLY parser builds into a *backend-neutral* descriptor — a list of atom
dictionaries plus a list of ``(i, j, order)`` bond tuples — and then drives the
chemistry engine of an external library to validate valence and aromaticity.

The whole point is that the backend never sees the SMILES string: we hand it a
graph that *we* constructed, and it performs only chemistry perception.  The
backend's own SMILES parser (``pysmiles.read_smiles``) is never called.
"""

from typing import List, Optional, Tuple

from chem.atomic import Atom, BracketAtom
from chem.structure import MolecularGraph

# Map our edge ``bond_type`` strings onto bond orders understood by the backend.
# Aromatic ``:`` -> 1.5, stereo ``/`` ``\`` collapse to single (they do not
# affect valence or aromaticity), ``.`` is no bond.
_BOND_ORDER = {"-": 1, "=": 2, "#": 3, "$": 4, ":": 1.5, "/": 1, "\\": 1, ".": 0}

AtomDescriptor = dict
BondDescriptor = Tuple[int, int, float]


def to_descriptor(
    graph: MolecularGraph,
) -> Tuple[List[AtomDescriptor], List[BondDescriptor]]:
    """Translate a ``MolecularGraph`` into a backend-neutral descriptor.

    Args:
        graph: The molecule graph produced by the parser/``GraphBuilder``.

    Returns:
        A tuple ``(atoms, bonds)`` where ``atoms`` is a list of dicts (one per
        atom, indexed by position) carrying ``element``, ``charge``,
        ``aromatic`` and — only for bracket atoms — an explicit ``hcount`` and
        optional ``isotope``; and ``bonds`` is a list of ``(i, j, order)``
        tuples with ``order`` in ``{1, 2, 3, 4, 1.5}``.

    Notes:
        Hydrogen handling is deliberately asymmetric and mirrors SMILES
        semantics: a **bracket** atom fully specifies its hydrogens, so we
        always emit ``hcount`` (``None`` means zero explicit H); an **organic**
        (non-bracket) atom has *implicit* hydrogens, so we omit ``hcount`` and
        let the backend fill the valence.  Disconnected components (the ``.``
        bond) simply produce no edge.
    """
    index = {atom: i for i, atom in enumerate(graph.adjacency_list)}

    atoms: List[AtomDescriptor] = []
    for atom in graph.adjacency_list:
        symbol = atom.symbol
        element = "*" if symbol == "*" else symbol.capitalize()
        descriptor: AtomDescriptor = {
            "element": element,
            "charge": getattr(atom, "charge", None) or 0,
            "aromatic": bool(getattr(atom, "aromatic", False)),
        }
        if isinstance(atom, BracketAtom):
            # Brackets fully specify hydrogens: None means zero explicit H.
            descriptor["hcount"] = atom.hcount if atom.hcount is not None else 0
            if atom.isotope is not None:
                descriptor["isotope"] = atom.isotope
        atoms.append(descriptor)

    bonds: List[BondDescriptor] = []
    seen = set()
    for atom, neighbours in graph.adjacency_list.items():
        for neighbour, bond_type in neighbours:
            if neighbour is atom:
                continue  # drop any self-loop artifact
            key = frozenset((index[atom], index[neighbour]))
            if key in seen:
                continue  # undirected: emit each edge once
            seen.add(key)
            order = _BOND_ORDER.get(bond_type, 1)
            if order == 0:
                continue  # dot disconnection -> no bond
            bonds.append((index[atom], index[neighbour], order))

    return atoms, bonds


def valid_pysmiles(
    atoms: List[AtomDescriptor], bonds: List[BondDescriptor]
) -> Tuple[bool, Optional[str]]:
    """Validate a descriptor using the pysmiles chemistry engine.

    Builds a ``networkx`` graph with the node/edge attributes pysmiles expects
    and calls its valence and aromaticity helpers directly, **bypassing**
    ``pysmiles.read_smiles``.

    Args:
        atoms: Atom descriptors from :func:`to_descriptor`.
        bonds: Bond descriptors from :func:`to_descriptor`.

    Returns:
        ``(is_valid, reason)``.  ``reason`` is ``None`` when valid, otherwise a
        short human-readable message.
    """
    # Imported lazily so importing this module does not hard-depend on the
    # backend until chemistry validation is actually requested.
    import networkx as nx
    from pysmiles import PTE
    from pysmiles.smiles_helper import (
        bonds_missing,
        correct_aromatic_rings,
        fill_valence,
    )

    mol = nx.Graph()
    for i, atom in enumerate(atoms):
        node = {
            "element": atom["element"],
            "charge": atom["charge"],
            "aromatic": atom["aromatic"],
        }
        if "hcount" in atom:
            node["hcount"] = atom["hcount"]
        if "isotope" in atom:
            node["isotope"] = atom["isotope"]
        mol.add_node(i, **node)
    for i, j, order in bonds:
        mol.add_edge(i, j, order=order)

    try:
        # Aromaticity perception + kekulization (raises if a region marked
        # aromatic cannot be kekulized).
        correct_aromatic_rings(mol, strict=True)
        # Fill implicit hydrogens from the per-element valence model.
        fill_valence(mol)
    except Exception as exc:  # noqa: BLE001 - backend signals invalidity by raising
        return False, f"aromaticity/kekulization failed: {exc}"

    for node in mol:
        element = mol.nodes[node].get("element", "*")
        if element == "*":
            continue
        if element not in PTE:
            return False, f"unknown element {element}"
        # A disconnected, charged atom is a free ion (e.g. [Na+], [Li+], [Cl-]):
        # it has no bonds to satisfy a valence, so the charge alone makes it
        # well-formed.  pysmiles' valence model is undefined for several such
        # noble-gas-configuration ions, so we accept them directly (matching
        # RDKit).
        if mol.degree(node) == 0 and mol.nodes[node].get("charge", 0) != 0:
            continue
        try:
            missing = bonds_missing(mol, node)
        except (ValueError, KeyError):
            # Valence undeterminable (e.g. some transition metals): no
            # constraint to apply, treat as acceptable.
            continue
        # ``bonds_missing`` is positive when the atom is *under* valence (a
        # radical / electron-deficient species) and negative when it is *over*
        # valence (beyond the element's maximum modelled valence, e.g.
        # pentavalent carbon).  This validator is deliberately permissive: it
        # accepts radicals and expanded-octet/hypervalent species (which RDKit
        # rejects) and only rejects the physically impossible over-valence.
        if missing < 0:
            return False, f"impossible valence on {element} (atom {node})"

    return True, None


# Map descriptor bond orders onto RDKit bond types (filled lazily).
_RDKIT_BOND_TYPE = None


def valid_rdkit(
    atoms: List[AtomDescriptor], bonds: List[BondDescriptor]
) -> Tuple[bool, Optional[str]]:
    """Validate a descriptor using RDKit's chemistry engine.

    Builds a ``Chem.RWMol`` atom-by-atom and bond-by-bond and runs
    ``Chem.SanitizeMol`` (valence + kekulization + aromaticity perception),
    **bypassing** ``Chem.MolFromSmiles``.

    Note:
        RDKit's sanitizer uses strict default valences, so this driver rejects
        the hypervalent / expanded-octet species that :func:`valid_pysmiles`
        accepts under the permissive policy.  It is therefore the *strict*
        (RDKit-equivalent) backend; choose it when maximum agreement with RDKit
        is wanted rather than permissiveness.

    Args:
        atoms: Atom descriptors from :func:`to_descriptor`.
        bonds: Bond descriptors from :func:`to_descriptor`.

    Returns:
        ``(is_valid, reason)``.  ``reason`` is ``None`` when valid.
    """
    global _RDKIT_BOND_TYPE
    from rdkit import Chem

    if _RDKIT_BOND_TYPE is None:
        _RDKIT_BOND_TYPE = {
            1: Chem.BondType.SINGLE,
            2: Chem.BondType.DOUBLE,
            3: Chem.BondType.TRIPLE,
            4: Chem.BondType.QUADRUPLE,
            1.5: Chem.BondType.AROMATIC,
        }

    rwmol = Chem.RWMol()
    idx = []
    for atom in atoms:
        element = atom["element"]
        rd_atom = Chem.Atom(0 if element == "*" else element)
        rd_atom.SetFormalCharge(atom.get("charge", 0))
        if atom.get("aromatic"):
            rd_atom.SetIsAromatic(True)
        if "hcount" in atom:
            rd_atom.SetNumExplicitHs(atom["hcount"])
            rd_atom.SetNoImplicit(True)
        if "isotope" in atom:
            rd_atom.SetIsotope(atom["isotope"])
        idx.append(rwmol.AddAtom(rd_atom))
    for i, j, order in bonds:
        rwmol.AddBond(idx[i], idx[j], _RDKIT_BOND_TYPE.get(order, Chem.BondType.SINGLE))
        if order == 1.5:
            rwmol.GetBondBetweenAtoms(idx[i], idx[j]).SetIsAromatic(True)

    mol = rwmol.GetMol()
    try:
        Chem.SanitizeMol(mol)  # valence + kekulize + set-aromaticity
    except Exception as exc:  # noqa: BLE001 - backend signals invalidity by raising
        return False, f"rdkit sanitize failed: {exc}"
    return True, None
