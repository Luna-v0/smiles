from sly import Parser
from validator.lex import SmilesLex
from src.chem.chem import chem
from itertools import combinations
from validator.parser_manager import parser_manager

def generate_combinations(rule: str) -> list[str]:
    """
    Generate all combinations of a rule with optional elements.
    Example: 
        rule = "X? Y Z?"
        combinations = [ X Y Z, X Y, Y Z, Y ]
    
    Args:
        rule: A string with the rule to generate combinations from.
    Returns:
        A list of strings with all combinations of the rule.
    """
    parts = rule.split()

    # Separate required and optional elements
    required = [p.rstrip('?') for p in parts if not p.endswith('?')]
    optional = [p.rstrip('?') for p in parts if p.endswith('?')]

    all_combinations = []

    for i in range(len(optional) + 1):
        for combo in combinations(optional, i):
            ordered_combo = [p.rstrip('?') for p in parts if p.rstrip(
                '?') in combo or p.rstrip('?') in required]
            all_combinations.append(" ".join(ordered_combo))

    return all_combinations


def getAttributes(rules, properties):
    """
    Get the attributes of the rules.
    Args:
        rules: The rules to get the attributes from.
        properties: The properties to get the attributes from.
    Returns:
        The attributes of the rules.
    """
    if type(properties) != list:
        return getattr(rules, properties, None)

    values = []

    for prop in properties:
        values.append(getAttributes(rules, prop))

    return values


class SmilesParser(Parser):
    """
    Parser using the SLY library to parse SMILES strings.
    """
    debugfile = 'parser.out'
    tokens = SmilesLex.tokens
    use_only_grammar = False

    def error(self, t):
        raise Exception(f"Error on {str(t)}")

    @_('atom', 'atom chain_branch')  # type: ignore
    def line(self, rules):
        return parser_manager.listify(*rules)        

    @_('chains', 'branch', 'chains chain_branch', 'branch chain_branch')  # type: ignore
    def chain_branch(self, rules):
        return parser_manager.listify(*rules)

    @_('chain', 'chain chains')  # type: ignore
    def chains(self, rules):
        if len(rules) == 1:
            return rules[0]

        if type(rules.chains) == list:
            return [rules[0]] + rules.chains

        return [rules[0], rules.chains]

    @_('"[" internal_bracket "]"')  # type: ignore
    def bracket_atom(self, rules):
        return rules.internal_bracket

    @_(*generate_combinations('isotope? symbol chiral? hcount? charge? map?')) # type: ignore
    def internal_bracket(self, rules):
        
        mol = getAttributes(rules, ['isotope', 'symbol', 'chiral','hcount', 'charge', 'map'])
        
        isotope, symbol, chiral, hcount, charge, mol_map = mol

        if self.use_only_grammar: return

        if not chem.validate_valency_bracket(isotope, symbol, chiral, hcount, charge, mol_map): 
            raise Exception(f"Invalid valency in Bracket [{','.join([str(x) for x in mol if x is not None])}]")
        
        

    @_('dot_proxy', 'bond atom', 'bond rnum', 'atom', 'rnum')  # type: ignore
    def chain(self, rules):
        if len(rules) == 1:
            return rules[0]
        
        if rules.bond == ':' and rules.atom and type(rules.atom) == str and rules.atom[0].isupper():
            raise Exception(f"Aromatic bond cannot be use with Uppercase and collon {rules.atom}")
        
        # TODO: need to check if the atom is not bracketed too
        
        return [rules.atom, chem.number_of_electrons_per_bond(rules.bond)]
    

    @_('"." atom')  # type: ignore
    def dot_proxy(self, rules):
        pass

    @_('semi_symbol', '"H"')  # type: ignore
    def symbol(self, rules):
        return rules[0]

    @_('"(" inner_branch ")"')  # type: ignore
    def branch(self, rules):
        return rules.inner_branch

    @_('bond_dot line', 'line', 'bond_dot line inner_branch', 'line inner_branch')    # type: ignore
    def inner_branch(self, rules):
        pass

    @_('bond', '"."')  # type: ignore
    def bond_dot(self, rules):
        return rules[0]

    @_('semi_bond_rule', '"-"')  # type: ignore
    def bond(self, rules):
        return rules[0]

    @_('semi_bond')  # type: ignore
    def semi_bond_rule(self, rules):
        return rules[0]

    @_('symbol', 'bracket_atom')  # type: ignore
    def atom(self, rules):
        return parser_manager.atom(rules[0])

    @_('digit', '"%" digit digit ')  # type: ignore
    def rnum(self, rules):
        return parser_manager.ring_number(*rules)

    @_('digit digit digit', 'digit digit', 'digit')  # type: ignore
    def isotope(self, rules):
        return parser_manager.int(*rules)

    @_('"H" digit', '"H"')  # type: ignore
    def hcount(self, rules):
        return parser_manager.hcount(*rules)

    @_('"+"', '"+" "+"', '"-"', '"-" "-"', '"-" fifteen', '"+" fifteen')  # type: ignore
    def charge(self, rules):
        return parser_manager.charge(*rules)

    @_('":" digit digit digit', '":" digit digit', '":" digit')  # type: ignore
    def map(self, rules):
        return parser_manager.int(*rules[1:])

    @_('"@"', '"@" "@"')  # type: ignore
    def chiral(self, rules):
        return parser_manager.chiral(*rules)

    @_('digit digit', 'digit')  # type: ignore
    def fifteen(self, rules):
        return parser_manager.fifteen(*rules)


parser = SmilesParser()
lexer = SmilesLex()


def validate_smiles(mol: str, use_only_grammar: bool = False) -> tuple[bool, Exception | None]:
    """
    Function for valdiating a SMILES molecule.
    
    Args:
        mol: Chemical formula as a string.
        use_only_grammar: For valdiating only the Grammar
        
    Returns:
        A Tuple containg in the first element if it is a valid SMILES and the second element a Exception.
    """
    try:
        parser.use_only_grammar = use_only_grammar
        parser.parse(lexer.tokenize(mol))
        return True, None
    except Exception as e:
        return False, e
