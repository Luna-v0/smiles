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


# Atomic-number ranges of the d- and f-block (transition metals, lanthanides,
# actinides).  Their valence depends on oxidation state and coordination, which
# no simple electron-count model captures, so no valence limit is applied.
_D_F_BLOCK_RANGES = ((21, 30), (39, 48), (57, 80), (89, 112))

# Second-row elements have no d orbitals and cannot expand their octet: only
# the lowest valence of the isoelectronic element is reachable (N 3, N+ 4,
# O 2, O- 1, ...).  pysmiles' electron-count model would otherwise also allow
# N 5 and O 4/6.
_NO_OCTET_EXPANSION = {"B", "C", "N", "O", "F"}

# Neutral donors whose bond to a d/f-block metal is read as a coordinate
# (dative) bond: it uses the donor's lone pair, so it does not count towards
# the donor's own valence (``[NH3][Pt]`` is cisplatin's ammine ligand).
_DATIVE_DONORS = {"N", "O"}

# Electron counts of the noble gases.  An ion isoelectronic with one of them
# (``[F-]``, ``[Cl-]``, ``[O-2]``) has a closed shell and cannot form bonds.
_NOBLE_GAS_ELECTRONS = {2, 10, 18, 36, 54, 86}


def _is_d_f_block(element: str) -> bool:
    """Return whether ``element`` is a transition metal, lanthanide or actinide."""
    from pysmiles import PTE

    number = PTE[element]["AtomicNumber"]
    return any(lo <= number <= hi for lo, hi in _D_F_BLOCK_RANGES)


def _valence_limit(mol, node, explicit_h: int) -> Tuple[float, float]:
    """Return ``(bond_total, max_valence)`` for one atom of a perceived graph.

    ``bond_total`` sums bond orders plus the explicit (bracket) hydrogens;
    implicit hydrogens of organic-subset atoms are never counted, since the
    backend only adds them up to a valid valence.

    Raises:
        ValueError, KeyError: If the backend cannot determine a valence.
    """
    from pysmiles import PTE
    from pysmiles.smiles_helper import valence

    attrs = mol.nodes[node]
    element = attrs["element"]
    charge = attrs.get("charge", 0)

    # An aromatic bond is one sigma bond; the atom's share of the delocalised
    # pi system adds one more, once (benzene C 2+1, naphthalene ring-fusion C
    # 3+1, pyridine N 2+1).  Summing 1.5 per aromatic bond would overcount.
    total = float(explicit_h)
    if any(mol.edges[node, n].get("order", 1) == 1.5 for n in mol[node]):
        total += 1
    for neighbour in mol[node]:
        neighbour_element = mol.nodes[neighbour].get("element", "*")
        if (
            element in _DATIVE_DONORS
            and charge == 0
            and neighbour_element != "*"
            and _is_d_f_block(neighbour_element)
        ):
            continue  # coordinate bond: the donor's lone pair, not a valence
        order = mol.edges[node, neighbour].get("order", 1)
        total += 1 if order == 1.5 else order

    if PTE[element]["AtomicNumber"] - charge in _NOBLE_GAS_ELECTRONS:
        return total, 0

    valences = valence(attrs)
    if not valences:
        raise ValueError(f"no valence model for {element}")
    if element in _NO_OCTET_EXPANSION:
        limit = valences[0]
        if _is_pentavalent_n_convention(mol, node, element, charge, explicit_h, total):
            limit = 5
        return total, limit
    return total, valences[-1]


def _is_pentavalent_n_convention(mol, node, element, charge, explicit_h, total) -> bool:
    """Neutral pentavalent N written without charge separation.

    ``CN(=O)=O`` (nitro), ``c1ccn(=O)cc1`` (N-oxide), ``C=N#N`` (diazo) and
    ``N=N#N`` (azide) are conventional shorthands for ``[N+]``/``[O-]`` or
    ``[N+]``/``[N-]`` forms; RDKit accepts them by charge-separating, so they
    are accepted here too.  A hydrogen on that N, or a fifth bond that is not
    a multiple bond to O/N, is still impossible.
    """
    if element != "N" or charge != 0 or explicit_h or total != 5:
        return False
    for neighbour in mol[node]:
        order = mol.edges[node, neighbour].get("order", 1)
        if order >= 2 and mol.nodes[neighbour].get("element") in ("O", "N"):
            return True
    return False


def _perceive(atoms: List[AtomDescriptor], bonds: List[BondDescriptor]):
    """Build the pysmiles graph and run aromaticity perception + H filling.

    Raises:
        Exception: Whatever the backend raises when a region marked aromatic
            cannot be kekulized.
    """
    import networkx as nx
    from pysmiles.smiles_helper import correct_aromatic_rings, fill_valence

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

    # Aromaticity perception + kekulization (raises if a region marked
    # aromatic cannot be kekulized).
    correct_aromatic_rings(mol, strict=True)
    # Fill implicit hydrogens from the per-element valence model.
    fill_valence(mol)
    return mol


def valence_violations(
    atoms: List[AtomDescriptor], bonds: List[BondDescriptor]
) -> List[Tuple[int, str]]:
    """List the atoms whose valence is impossible, with a reason for each.

    Raises:
        Exception: If aromaticity perception fails (see :func:`_perceive`).
    """
    from pysmiles import PTE

    mol = _perceive(atoms, bonds)
    violations: List[Tuple[int, str]] = []
    for node in mol:
        element = mol.nodes[node].get("element", "*")
        if element == "*":
            continue
        if element not in PTE:
            violations.append((node, f"unknown element {element}"))
            continue
        # A disconnected, charged atom is a free ion (e.g. [Na+], [Li+], [Cl-]):
        # it has no bonds to satisfy a valence, so the charge alone makes it
        # well-formed.  pysmiles' valence model is undefined for several such
        # noble-gas-configuration ions, so we accept them directly (matching
        # RDKit).
        if mol.degree(node) == 0 and mol.nodes[node].get("charge", 0) != 0:
            continue
        if _is_d_f_block(element):
            continue
        try:
            total, limit = _valence_limit(mol, node, atoms[node].get("hcount", 0))
        except (ValueError, KeyError):
            # Valence undeterminable: no constraint to apply, treat as acceptable.
            continue
        # Only over-valence is rejected.  This validator is deliberately
        # permissive: it accepts radicals (under-valence) and the
        # expanded-octet/hypervalent species of period 3+ (which RDKit rejects).
        if total > limit:
            charge = mol.nodes[node].get("charge", 0)
            charge_text = f"{charge:+d}" if charge else "neutral"
            violations.append(
                (
                    node,
                    f"impossible valence on {element} (atom {node}): "
                    f"{total:g} bond{'s' if total != 1 else ''}, {charge_text} {element} "
                    f"allows at most {limit}",
                )
            )
    return violations


def open_shell_atoms(
    atoms: List[AtomDescriptor], bonds: List[BondDescriptor]
) -> List[Tuple[int, int]]:
    """List atoms left with an odd number of unpaired ground-state electrons.

    This is the opt-in **closed-shell check** (see
    ``docs/validation.md``, "Closed-shell check"): it never feeds into
    :func:`valence_violations` and never changes a verdict, it only tells the
    caller which atoms are radicals.

    An atom's bond total (bonds + explicit H, aromaticity-aware, dative bonds
    to d/f-block metals excluded — the same ``total`` :func:`_valence_limit`
    computes) is compared against its *ground-state* valence (``valence()[0]``
    from pysmiles — the lowest, unexcited value, the same one used as the
    ceiling for second-row atoms). A shortfall is a lone pair when it is even
    (the unfilled electrons can pair with each other, e.g. a singlet carbene)
    and a forced radical when it is odd (one electron has no partner to pair
    with). Only the odd case is reported.

    Exempt, matching :func:`valence_violations`: wildcard atoms, elements
    with no entry in the periodic table, free (disconnected) charged ions,
    closed-shell (noble-gas-isoelectronic) ions, and d/f-block metals (no
    simple valence model applies to any of them).  An atom already reported
    by :func:`valence_violations` (over-valent) is never open-shell: when the
    bond total reaches or exceeds the ground-state valence there is no
    shortfall to report, expanded-octet species included.

    Args:
        atoms: Atom descriptors from :func:`to_descriptor`.
        bonds: Bond descriptors from :func:`to_descriptor`.

    Returns:
        A list of ``(atom_index, unpaired_electron_count)`` pairs, in atom
        order. ``unpaired_electron_count`` is always odd and positive.

    Raises:
        Exception: If aromaticity perception fails (see :func:`_perceive`).
    """
    from pysmiles import PTE
    from pysmiles.smiles_helper import valence

    mol = _perceive(atoms, bonds)
    flagged: List[Tuple[int, int]] = []
    for node in mol:
        attrs = mol.nodes[node]
        element = attrs.get("element", "*")
        if element == "*" or element not in PTE:
            continue
        if mol.degree(node) == 0 and attrs.get("charge", 0) != 0:
            continue
        if _is_d_f_block(element):
            continue
        try:
            # Unlike valence_violations (over-valence only, where ignoring an
            # organic atom's implicit H is a safe lower bound), the implicit
            # H an organic atom was filled with here (fill_valence, inside
            # _perceive) must be counted: it is what already closes that
            # atom's shell, and skipping it would flag every ordinary
            # organic-subset atom as a false radical.
            total, _ = _valence_limit(mol, node, mol.nodes[node].get("hcount", 0))
        except (ValueError, KeyError):
            continue
        if PTE[element]["AtomicNumber"] - attrs.get("charge", 0) in _NOBLE_GAS_ELECTRONS:
            continue  # closed-shell ion: no valence to fill, never open-shell
        ground = valence(attrs)
        if not ground:
            continue
        deficit = ground[0] - total
        if deficit > 0 and int(deficit) % 2 == 1:
            flagged.append((node, int(deficit)))
    return flagged


def valid_pysmiles(
    atoms: List[AtomDescriptor], bonds: List[BondDescriptor]
) -> Tuple[bool, Optional[str]]:
    """Validate a descriptor using the pysmiles chemistry engine.

    Builds a ``networkx`` graph with the node/edge attributes pysmiles expects
    and calls its valence and aromaticity helpers directly, **bypassing**
    ``pysmiles.read_smiles``.  pysmiles' electron-count valence model is
    tightened where it is chemically wrong: second-row atoms cannot expand
    their octet, closed-shell ions cannot bond, and d/f-block metals get no
    limit (see :func:`valence_violations`).

    Args:
        atoms: Atom descriptors from :func:`to_descriptor`.
        bonds: Bond descriptors from :func:`to_descriptor`.

    Returns:
        ``(is_valid, reason)``.  ``reason`` is ``None`` when valid, otherwise a
        short human-readable message.
    """
    try:
        violations = valence_violations(atoms, bonds)
    except ImportError:
        raise  # a missing backend is an internal error, not a verdict
    except Exception as exc:  # noqa: BLE001 - backend signals invalidity by raising
        return False, f"aromaticity/kekulization failed: {exc}"
    if violations:
        return False, violations[0][1]
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
