from smiles_checker.chem.structure import Graph
from smiles_checker.chem.atomic import Atom, BracketAtom
from smiles_checker.chem.chemistry import chemistry

def test_add_edge():
    graph = Graph()
    atom1 = Atom("C")
    atom2 = Atom("N")
    graph.add_edge(atom1, atom2, bond_type="-")
    assert atom2 in [a for a, b in graph.adjacency_list[atom1]]
    assert atom1 in [a for a, b in graph.adjacency_list[atom2]]

def test_get_acyclic_subgraphs():
    graph = Graph()
    atom1 = Atom("C")
    atom2 = Atom("N")
    atom3 = Atom("O")
    atom4 = Atom("P")
    graph.add_edge(atom1, atom2, bond_type="-")
    graph.add_edge(atom3, atom4, bond_type="-")
    subgraphs = graph.get_acyclic_subgraphs()
    assert len(subgraphs) == 2
    assert [atom1, atom2] in subgraphs or [atom2, atom1] in subgraphs
    assert [atom3, atom4] in subgraphs or [atom4, atom3] in subgraphs

def test_add_cycle():
    graph = Graph()
    atom1 = Atom("C")
    atom2 = Atom("N")
    atom3 = Atom("O")
    cycle = [atom1, atom2, atom3]
    graph.add_cycle(cycle)
    assert cycle in graph.cycles

def test_check_valency_for_aba():
    graph = Graph()
    atom1 = chemistry.BracketAtom("C", hidrogens=3) # Methane, should be valid
    atom2 = chemistry.BracketAtom("O", hidrogens=0) # Hydroxyl, should be valid
    atom3 = chemistry.BracketAtom("N", hidrogens=1) # Amine, should be valid
    atom4 = chemistry.BracketAtom("F") # Fluorine, should be valid

    graph.add_edge(atom1, atom2, bond_type="-")
    graph.add_edge(atom2, atom3, bond_type="-")
    graph.add_edge(atom3, atom4, bond_type="-")

    assert graph.check_valency_for_aba() == True

    atom5 = chemistry.BracketAtom("C", hidrogens=1) # Invalid carbon
    graph.add_edge(atom4, atom5, bond_type="-")
    assert graph.check_valency_for_aba() == False

def test_huckel():
    graph = Graph()
    # Benzene ring (aromatic)
    c1 = chemistry.BracketAtom("C", aromatic=True)
    c2 = chemistry.BracketAtom("C", aromatic=True)
    c3 = chemistry.BracketAtom("C", aromatic=True)
    c4 = chemistry.BracketAtom("C", aromatic=True)
    c5 = chemistry.BracketAtom("C", aromatic=True)
    c6 = chemistry.BracketAtom("C", aromatic=True)

    graph.add_edge(c1, c2, bond_type=":")
    graph.add_edge(c2, c3, bond_type=":")
    graph.add_edge(c3, c4, bond_type=":")
    graph.add_edge(c4, c5, bond_type=":")
    graph.add_edge(c5, c6, bond_type=":")
    graph.add_edge(c6, c1, bond_type=":")

    graph.add_cycle([c1, c2, c3, c4, c5, c6])

    assert graph.huckel() == True

    # Cyclobutadiene (anti-aromatic)
    c7 = chemistry.BracketAtom("C", aromatic=True)
    c8 = chemistry.BracketAtom("C", aromatic=True)
    c9 = chemistry.BracketAtom("C", aromatic=True)
    c10 = chemistry.BracketAtom("C", aromatic=True)

    graph2 = Graph()
    graph2.add_edge(c7, c8, bond_type=":")
    graph2.add_edge(c8, c9, bond_type=":")
    graph2.add_edge(c9, c10, bond_type=":")
    graph2.add_edge(c10, c7, bond_type=":")

    graph2.add_cycle([c7, c8, c9, c10])

    assert graph2.huckel() == False

def test_get_acyclic_subgraphs_single_node():
    graph = Graph()
    atom1 = Atom("C")
    graph.add_edge(atom1, atom1, bond_type="-") # Self-loop, but still a single node for DFS
    subgraphs = graph.get_acyclic_subgraphs()
    assert len(subgraphs) == 1
    assert subgraphs[0] == [atom1]

def test_get_acyclic_subgraphs_disconnected_nodes():
    graph = Graph()
    atom1 = Atom("C")
    atom2 = Atom("N")
    atom3 = Atom("O")
    graph.add_edge(atom1, atom1, bond_type="-")
    graph.add_edge(atom2, atom2, bond_type="-")
    graph.add_edge(atom3, atom3, bond_type="-")
    subgraphs = graph.get_acyclic_subgraphs()
    assert len(subgraphs) == 3
    assert [atom1] in subgraphs
    assert [atom2] in subgraphs
    assert [atom3] in subgraphs

def test_get_acyclic_subgraphs_connected_nodes():
    graph = Graph()
    atom1 = Atom("C")
    atom2 = Atom("N")
    atom3 = Atom("O")
    graph.add_edge(atom1, atom2, bond_type="-")
    graph.add_edge(atom2, atom3, bond_type="-")
    subgraphs = graph.get_acyclic_subgraphs()
    assert len(subgraphs) == 1
    assert len(subgraphs[0]) == 3
    assert atom1 in subgraphs[0]
    assert atom2 in subgraphs[0]
    assert atom3 in subgraphs[0]

def test_get_acyclic_subgraphs_with_cycle():
    graph = Graph()
    atom1 = Atom("C")
    atom2 = Atom("N")
    atom3 = Atom("O")
    graph.add_edge(atom1, atom2, bond_type="-")
    graph.add_edge(atom2, atom3, bond_type="-")
    graph.add_edge(atom3, atom1, bond_type="-") # Create a cycle
    subgraphs = graph.get_acyclic_subgraphs()
    assert len(subgraphs) == 1
    assert len(subgraphs[0]) == 3
    assert atom1 in subgraphs[0]
    assert atom2 in subgraphs[0]
    assert atom3 in subgraphs[0]

def test_check_valency_for_aba_non_bracket_atom():
    graph = Graph()
    atom1 = Atom("C")
    graph.add_edge(atom1, atom1, bond_type="-")
    assert graph.check_valency_for_aba() == True

def test_huckel_non_aromatic_cycle():
    graph = Graph()
    atom1 = Atom("C")
    atom2 = Atom("C")
    atom3 = Atom("C")
    atom4 = Atom("C")
    graph.add_edge(atom1, atom2, bond_type=":")
    graph.add_edge(atom2, atom3, bond_type=":")
    graph.add_edge(atom3, atom4, bond_type=":")
    graph.add_edge(atom4, atom1, bond_type=":")
    graph.add_cycle([atom1, atom2, atom3, atom4])
    assert graph.huckel() == False

def test_get_acyclic_subgraphs_complex():
    graph = Graph()
    atom1 = Atom("C")
    atom2 = Atom("N")
    atom3 = Atom("O")
    atom4 = Atom("P")
    atom5 = Atom("S")

    graph.add_edge(atom1, atom2, bond_type="-")
    graph.add_edge(atom2, atom3, bond_type="-")
    graph.add_edge(atom3, atom1, bond_type="-") # Cycle 1

    graph.add_edge(atom4, atom5, bond_type="-") # Disconnected component

    subgraphs = graph.get_acyclic_subgraphs()
    assert len(subgraphs) == 2
    # Check if the subgraphs contain the correct atoms
    assert (len(subgraphs[0]) == 3 and atom1 in subgraphs[0] and atom2 in subgraphs[0] and atom3 in subgraphs[0]) or \
           (len(subgraphs[1]) == 3 and atom1 in subgraphs[1] and atom2 in subgraphs[1] and atom3 in subgraphs[1])
    assert (len(subgraphs[0]) == 2 and atom4 in subgraphs[0] and atom5 in subgraphs[0]) or \
           (len(subgraphs[1]) == 2 and atom4 in subgraphs[1] and atom5 in subgraphs[1])

def test_check_valency_for_aba_empty_graph():
    graph = Graph()
    assert graph.check_valency_for_aba() == True

def test_huckel_no_cycles():
    graph = Graph()
    atom1 = Atom("C")
    atom2 = Atom("N")
    graph.add_edge(atom1, atom2, bond_type="-")
    assert graph.huckel() == True

def test_check_valency_for_aba_mixed_atoms():
    graph = Graph()
    atom1 = Atom("C")
    atom2 = chemistry.BracketAtom("O", hidrogens=0)
    atom3 = Atom("N")
    atom4 = chemistry.BracketAtom("F")

    graph.add_edge(atom1, atom2, bond_type="-")
    graph.add_edge(atom2, atom3, bond_type="-")
    graph.add_edge(atom3, atom4, bond_type="-")

    assert graph.check_valency_for_aba() == True

def test_huckel_non_huckel_rule():
    graph = Graph()
    c1 = chemistry.BracketAtom("C", aromatic=True)
    c2 = chemistry.BracketAtom("C", aromatic=True)
    c3 = chemistry.BracketAtom("C", aromatic=True)
    c4 = chemistry.BracketAtom("C", aromatic=True)
    c5 = chemistry.BracketAtom("C", aromatic=True)

    graph.add_edge(c1, c2, bond_type=":")
    graph.add_edge(c2, c3, bond_type=":")
    graph.add_edge(c3, c4, bond_type=":")
    graph.add_edge(c4, c5, bond_type=":")
    graph.add_edge(c5, c1, bond_type=":")

    graph.add_cycle([c1, c2, c3, c4, c5])

    assert graph.huckel() == False

def test_get_acyclic_subgraphs_linear_chain():
    graph = Graph()
    atom1 = Atom("A")
    atom2 = Atom("B")
    atom3 = Atom("C")
    atom4 = Atom("D")
    graph.add_edge(atom1, atom2, bond_type="-")
    graph.add_edge(atom2, atom3, bond_type="-")
    graph.add_edge(atom3, atom4, bond_type="-")
    subgraphs = graph.get_acyclic_subgraphs()
    assert len(subgraphs) == 1
    assert len(subgraphs[0]) == 4
    assert atom1 in subgraphs[0]
    assert atom2 in subgraphs[0]
    assert atom3 in subgraphs[0]
    assert atom4 in subgraphs[0]