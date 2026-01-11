from itertools import combinations
from sly import Lexer, Parser
import re


# =============================================================================
# LEXER
# =============================================================================

PT_SYMBOLS = [
    'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si',
    'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',
    'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',
    'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb',
    'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu',
    'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re',
    'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr',
    'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es',
    'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg',
    'Cn', 'Fl', 'Lv'
]

AROMATIC_SYMBOLS = ['b', 'c', 'n', 'o', 'p', 's', 'se', 'as']


def generate_regex_from_list(elem_list):
    return "|".join(re.escape(e) for e in elem_list)


def generate_all_symbols(elem_list):
    result = []
    for elem in elem_list:
        result.append(elem)
        result.append(elem.lower())
    return sorted(set(result), key=len, reverse=True)


# Pre-compute the regex pattern
_ALL_SYMBOLS = generate_all_symbols(PT_SYMBOLS) + AROMATIC_SYMBOLS
_SEMI_SYMBOL_PATTERN = generate_regex_from_list(_ALL_SYMBOLS)


class SmilesLex(Lexer):
    """Tokenizer for SMILES strings."""
    
    literals = {'.', '@', '-', '+', ':', '%', 'H', ')', '(', ']', '['}
    tokens = {'semi_bond', 'digit', 'semi_symbol'}
    
    semi_symbol = _SEMI_SYMBOL_PATTERN
    semi_bond = r'[=#$/\\]'
    
    @_(r'\d')
    def digit(self, t):
        t.value = int(t.value)
        return t


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def generate_combinations(rule: str) -> list[str]:
    """Generate all combinations for optional elements in a rule."""
    parts = rule.split()
    required = [p.rstrip("?") for p in parts if not p.endswith("?")]
    optional = [p.rstrip("?") for p in parts if p.endswith("?")]

    all_combinations = []
    for i in range(len(optional) + 1):
        for combo in combinations(optional, i):
            ordered_combo = [
                p.rstrip("?") for p in parts
                if p.rstrip("?") in combo or p.rstrip("?") in required
            ]
            all_combinations.append(" ".join(ordered_combo))
    return all_combinations


# =============================================================================
# AST NODES
# =============================================================================

class ASTNode:
    def __repr__(self):
        attrs = {k: v for k, v in self.__dict__.items() if v is not None}
        return f"{self.__class__.__name__}({attrs})"


class Chain(ASTNode):
    def __init__(self, elements=None):
        self.elements = elements or []
    
    def append(self, elem):
        self.elements.append(elem)
        return self


class ChainElem(ASTNode):
    def __init__(self, kind, value, bond=None):
        self.kind = kind
        self.value = value
        self.bond = bond


class Atom(ASTNode):
    def __init__(self, symbol=None, bracket_atom=None):
        self.symbol = symbol
        self.bracket_atom = bracket_atom


class BracketAtom(ASTNode):
    def __init__(self, isotope=None, symbol=None, chiral=None, 
                 hcount=None, charge=None, mol_map=None):
        self.isotope = isotope
        self.symbol = symbol
        self.chiral = chiral
        self.hcount = hcount
        self.charge = charge
        self.mol_map = mol_map


class Branch(ASTNode):
    def __init__(self, bond=None, chain=None):
        self.bond = bond
        self.chain = chain


class RingNum(ASTNode):
    def __init__(self, value):
        self.value = value


class Bond(ASTNode):
    def __init__(self, symbol):
        self.symbol = symbol


# =============================================================================
# ZERO-CONFLICT PARSER
# =============================================================================

class SmilesParser(Parser):
    """
    SMILES Parser with ZERO shift/reduce conflicts.
    
    Grammar:
    --------
    
    smiles      -> chain
    
    chain       -> atom                     (base case)
                 | chain chain_elem         (LEFT recursion)
    
    chain_elem  -> rnum                     (ring number)
                 | branch                   (branch)
                 | bond atom                (bonded atom)
                 | bond rnum                (bonded ring)
                 | DOT atom                 (disconnection)
                 | atom                     (implicit bond)
    
    The key is that:
    1. chain -> atom can ONLY reduce when no chain_elem can follow
    2. chain -> chain chain_elem shifts the chain_elem and then reduces
    3. No ambiguity about whether to stop or continue
    """

    debugfile = "parser_zero_conflict.out"
    tokens = SmilesLex.tokens
    
    # Precedence to resolve any remaining potential ambiguity
    # Lower in list = higher precedence
    # This tells the parser: when in doubt, shift these tokens
    precedence = (
        ('left', '.', '-', 'semi_bond', '[', 'H', 'semi_symbol', 'digit', '%', '('),
    )
    
    def error(self, t):
        if t:
            raise SyntaxError(f"Syntax error at '{t.value}' (type: {t.type}) position {t.index}")
        else:
            raise SyntaxError("Syntax error at end of input")

    # =========================================================================
    # CHAIN (top-level)
    # =========================================================================
    
    @_("chain chain_elem")
    def chain(self, p):
        """Extend the chain with another element (LEFT recursive)."""
        p.chain.append(p.chain_elem)
        return p.chain
    
    @_("atom")
    def chain(self, p):
        """Base case: a single atom starts a chain."""
        return Chain([ChainElem('atom', p.atom)])

    # =========================================================================
    # CHAIN_ELEM - elements that extend a chain
    # =========================================================================
    
    @_("rnum")
    def chain_elem(self, p):
        """Ring closure."""
        return ChainElem('rnum', p.rnum)
    
    @_("branch")
    def chain_elem(self, p):
        """Branch."""
        return ChainElem('branch', p.branch)
    
    @_("bond atom")
    def chain_elem(self, p):
        """Explicit bond to atom."""
        return ChainElem('atom', p.atom, bond=p.bond)
    
    @_("bond rnum")
    def chain_elem(self, p):
        """Explicit bond in ring closure."""
        return ChainElem('rnum', p.rnum, bond=p.bond)
    
    @_('"." atom')
    def chain_elem(self, p):
        """Disconnection (dot)."""
        return ChainElem('dot', p.atom)
    
    @_("atom")
    def chain_elem(self, p):
        """Adjacent atom (implicit single bond)."""
        return ChainElem('atom', p.atom)

    # =========================================================================
    # ATOM
    # =========================================================================
    
    @_("symbol")
    def atom(self, p):
        return Atom(symbol=p.symbol)

    @_("bracket_atom")
    def atom(self, p):
        return Atom(bracket_atom=p.bracket_atom)

    # =========================================================================
    # BRACKET_ATOM
    # =========================================================================
    
    @_('"[" internal_bracket "]"')
    def bracket_atom(self, p):
        return p.internal_bracket

    @_(*generate_combinations("isotope? symbol chiral? hcount? charge? mol_map?"))
    def internal_bracket(self, p):
        return BracketAtom(
            isotope=getattr(p, 'isotope', None),
            symbol=getattr(p, 'symbol', None),
            chiral=getattr(p, 'chiral', None),
            hcount=getattr(p, 'hcount', None),
            charge=getattr(p, 'charge', None),
            mol_map=getattr(p, 'mol_map', None)
        )

    # =========================================================================
    # SYMBOL
    # =========================================================================
    
    @_("semi_symbol")
    def symbol(self, p):
        return p.semi_symbol
    
    @_('"H"')
    def symbol(self, p):
        return 'H'

    # =========================================================================
    # BRANCH
    # Branches can start with an optional bond/dot
    # =========================================================================
    
    @_('"(" chain ")"')
    def branch(self, p):
        """Branch without explicit bond."""
        return Branch(chain=p.chain)
    
    @_('"(" bond chain ")"')
    def branch(self, p):
        """Branch with explicit bond."""
        return Branch(bond=p.bond, chain=p.chain)
    
    @_('"(" "." chain ")"')
    def branch(self, p):
        """Branch with disconnection."""
        return Branch(bond='.', chain=p.chain)

    # =========================================================================
    # BOND
    # =========================================================================
    
    @_("semi_bond")
    def bond(self, p):
        return Bond(p.semi_bond)
    
    @_('"-"')
    def bond(self, p):
        return Bond('-')

    # =========================================================================
    # RING NUMBER
    # =========================================================================
    
    @_("digit")
    def rnum(self, p):
        return RingNum(p.digit)

    @_('"%" digit digit')
    def rnum(self, p):
        return RingNum(10 * p[1] + p[2])

    # =========================================================================
    # ISOTOPE
    # =========================================================================
    
    @_("digit")
    def isotope(self, p):
        return p.digit

    @_("digit digit")
    def isotope(self, p):
        return 10 * p[0] + p[1]

    @_("digit digit digit")
    def isotope(self, p):
        return 100 * p[0] + 10 * p[1] + p[2]

    # =========================================================================
    # HCOUNT
    # =========================================================================
    
    @_('"H"')
    def hcount(self, p):
        return 1

    @_('"H" digit')
    def hcount(self, p):
        return p.digit

    # =========================================================================
    # CHARGE
    # =========================================================================
    
    @_('"+"')
    def charge(self, p):
        return 1

    @_('"-"')
    def charge(self, p):
        return -1

    @_('"+" fifteen')
    def charge(self, p):
        return p.fifteen

    @_('"-" fifteen')
    def charge(self, p):
        return -p.fifteen

    @_('"+" "+"')
    def charge(self, p):
        return 2

    @_('"-" "-"')
    def charge(self, p):
        return -2

    # =========================================================================
    # MOL_MAP
    # =========================================================================
    
    @_('":" digit')
    def mol_map(self, p):
        return p.digit

    @_('":" digit digit')
    def mol_map(self, p):
        return 10 * p[1] + p[2]

    @_('":" digit digit digit')
    def mol_map(self, p):
        return 100 * p[1] + 10 * p[2] + p[3]

    # =========================================================================
    # CHIRAL
    # =========================================================================
    
    @_('"@"')
    def chiral(self, p):
        return '@'

    @_('"@" "@"')
    def chiral(self, p):
        return '@@'

    # =========================================================================
    # FIFTEEN
    # =========================================================================
    
    @_("digit")
    def fifteen(self, p):
        return p.digit

    @_("digit digit")
    def fifteen(self, p):
        return 10 * p[0] + p[1]


# =============================================================================
# TESTING
# =============================================================================

def validate_smiles(mol: str) -> tuple:
    lexer = SmilesLex()
    parser = SmilesParser()
    
    try:
        result = parser.parse(lexer.tokenize(mol))
        return True, result
    except Exception as e:
        return False, e


def test_parser():
    """Test the parser with various SMILES strings."""
    test_cases = [
        # Simple atoms
        ("C", True),
        ("N", True),
        ("O", True),
        ("H", True),
        
        # Chains
        ("CC", True),
        ("CCC", True),
        ("CCCC", True),
        
        # Bonds
        ("C=C", True),
        ("C#C", True),
        ("C=CC", True),
        ("C-C", True),
        
        # Ring closures
        ("C1CC1", True),
        ("C1CCCCC1", True),
        ("C%12CC%12", True),
        ("C=1CC=1", True),
        
        # Branches
        ("CC(C)C", True),
        ("CC(C)(C)C", True),
        ("CC(=O)O", True),
        ("C(C)(C)C", True),
        
        # Bracket atoms
        ("[C]", True),
        ("[CH4]", True),
        ("[13C]", True),
        ("[C@H]", True),
        ("[C@@H]", True),
        ("[Cu+2]", True),
        ("[O-]", True),
        ("[NH4+]", True),
        ("[C:1]", True),
        
        # Complex examples
        ("CC(C)CC", True),
        ("c1ccccc1", True),
        ("C1=CC=CC=C1", True),
        ("CCO", True),
        ("CC(=O)O", True),
        ("CC(=O)OC", True),
        ("c1ccc2ccccc2c1", True),
        
        # Disconnected structures
        ("C.C", True),
        ("[Na+].[Cl-]", True),
        
        # Edge cases
        ("C(C)C", True),
        ("C(=O)(O)C", True),
    ]
    
    print("Testing SMILES Parser (Zero Conflicts)")
    print("=" * 60)
    
    passed = 0
    failed = 0
    
    for smiles, expected in test_cases:
        valid, result = validate_smiles(smiles)
        
        if valid == expected:
            status = "✓"
            passed += 1
        else:
            status = "✗"
            failed += 1
            
        print(f"{status} {smiles}")
        if not valid and expected:
            print(f"  Error: {result}")
    
    print("=" * 60)
    print(f"Results: {passed} passed, {failed} failed")
    
    return failed == 0


if __name__ == "__main__":
    test_parser()

print(validate_smiles("[H]c1c([H])c([H])c([H])c([H])c1[H]"))

