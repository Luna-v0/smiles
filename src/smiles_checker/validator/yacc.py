from itertools import combinations

from sly import Parser

from smiles_checker.chem.chemistry import chemistry as chem
from smiles_checker.exceptions import ParserException
from smiles_checker.validator.lex import SmilesLex

from .parser_manager import parser_manager as pm


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
    required = [p.rstrip("?") for p in parts if not p.endswith("?")]
    optional = [p.rstrip("?") for p in parts if p.endswith("?")]

    all_combinations = []

    for i in range(len(optional) + 1):
        for combo in combinations(optional, i):
            ordered_combo = [
                p.rstrip("?")
                for p in parts
                if p.rstrip("?") in combo or p.rstrip("?") in required
            ]
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

    debugfile = "parser.out"
    tokens = SmilesLex.tokens
    use_only_grammar = False

    precedence = (
        ("right", "."),
        ("right", "-"),
        ("right", "%"),
        ("right", "digit"),
        ("right", "semi_bond"),
    )

    def error(self, t):
        raise Exception(f"Error on {str(t)}")

    @_("atom")  # type: ignore
    def line(self, rules):
        return pm.line(atom=rules.atom)

    @_("atom chain_branch")  # type: ignore
    def line(self, rules):
        return pm.line(atom=rules.atom, chain_branch=rules.chain_branch)

    @_("chain_branch chain_branch_item %prec semi_bond")  # type: ignore
    def chain_branch(self, rules):
        return pm.chain_branch(chains=rules.chains, chain_branch=rules.chain_branch)

    @_("chain_branch_item")  # type: ignore
    def chain_branch(self, rules):
        return pm.chain_branch(chain_branch=rules.chain_branch_item)

    @_("chain")  # type: ignore
    def chains(self, rules):
        return pm.chains(chain=rules.chain)

    @_("chains chain")  # type: ignore
    def chains(self, rules):
        return pm.chains(chain=rules.chain, chains=rules.chains)

    @_('"[" internal_bracket "]"')  # type: ignore
    def bracket_atom(self, rules):
        return pm.bracket_atom(internal_bracket=rules.internal_bracket)

    @_(*generate_combinations("isotope? symbol chiral? hcount? charge? mol_map?"))  # type: ignore
    def internal_bracket(self, rules):
        attrs = getAttributes(
            rules, ["isotope", "symbol", "chiral", "hcount", "charge", "mol_map"]
        )
        kwargs = {
            "isotope": attrs[0],
            "symbol": attrs[1],
            "chiral": attrs[2],
            "hidrogens": attrs[3],
            "charge": attrs[4],
            "mol_map": attrs[5],
        }
        # Create a dictionary from the list of attributes and their corresponding keys
        return pm.internal_bracket(**kwargs)

    @_("dot_proxy")  # type: ignore
    def chain(self, rules):
        return pm.chain(dot_proxy=rules.dot_proxy)

    @_("bond atom")  # type: ignore
    def chain(self, rules):
        return pm.chain(bond=rules.bond, atom=rules.atom)

    @_("bond rnum")  # type: ignore
    def chain(self, rules):
        return pm.chain(bond=rules.bond, rnum=rules.rnum)

    @_("rnum")  # type: ignore
    def chain(self, rules):
        return pm.chain(rnum=rules.rnum)

    @_('"." atom')  # type: ignore
    def dot_proxy(self, rules):
        return pm.dot_proxy(atom=rules.atom)

    @_("semi_symbol", '"H"')  # type: ignore
    def symbol(self, rules):
        return pm.symbol(semi_symbol=rules[0])

    @_('"(" inner_branch ")"')  # type: ignore
    def branch(self, rules):
        return pm.branch(inner_branch=rules.inner_branch)

    @_("bond_dot line")  # type: ignore
    def inner_branch(self, rules):
        return pm.inner_branch(bond_dot=rules.bond_dot, line=rules.line)

    @_("line")  # type: ignore
    def inner_branch(self, rules):
        return pm.inner_branch(line=rules.line)

    @_("bond_dot line inner_branch")  # type: ignore
    def inner_branch(self, rules):
        return pm.inner_branch(
            bond_dot=rules.bond_dot, line=rules.line, inner_branch=rules.inner_branch
        )

    @_("line inner_branch")  # type: ignore
    def inner_branch(self, rules):
        return pm.inner_branch(line=rules.line, inner_branch=rules.inner_branch)

    @_("bond", '"."')  # type: ignore
    def bond_dot(self, rules):
        return pm.bond_dot(bond=rules[0])

    @_("semi_bond_rule", '"-"')  # type: ignore
    def bond(self, rules):
        return pm.bond(semi_bond_rule=rules[0])

    @_("semi_bond")  # type: ignore
    def semi_bond_rule(self, rules):
        return pm.semi_bond_rule(semi_bond=rules.semi_bond)

    @_("symbol")  # type: ignore
    def atom(self, rules):
        return pm.atom(symbol=rules.symbol)

    @_("bracket_atom")  # type: ignore
    def atom(self, rules):
        return pm.atom(bracket_atom=rules.bracket_atom)

    @_("digit")  # type: ignore
    def rnum(self, rules):
        return pm.rnum(digit=rules.digit)

    @_('"%" digit digit ')  # type: ignore
    def rnum(self, rules):
        return pm.rnum(digit=(rules[1], rules[2]))

    @_("digit digit digit")  # type: ignore
    def isotope(self, rules):
        return pm.isotope(digit=(rules[0], rules[1], rules[2]))

    @_("digit digit")  # type: ignore
    def isotope(self, rules):
        return pm.isotope(digit=(rules[0], rules[1]))

    @_("digit")  # type: ignore
    def isotope(self, rules):
        return pm.isotope(digit=rules.digit)

    @_('"H" digit')  # type: ignore
    def hcount(self, rules):
        return pm.hcount(digit=rules.digit)

    @_('"H"')  # type: ignore
    def hcount(self, rules):
        return pm.hcount()

    @_('"+"')  # type: ignore
    def charge(self, rules):
        return pm.charge(charge=1)

    @_('"+" "+"')  # type: ignore
    def charge(self, rules):
        return pm.charge(charge=2)

    @_('"-"')  # type: ignore
    def charge(self, rules):
        return pm.charge(charge=-1)

    @_('"-" "-"')  # type: ignore
    def charge(self, rules):
        return pm.charge(charge=-2)

    @_('"-" fifteen')  # type: ignore
    def charge(self, rules):
        return pm.charge(charge=-rules.fifteen)

    @_('"+" fifteen')  # type: ignore
    def charge(self, rules):
        return pm.charge(rules)

    @_('":" digit digit digit')  # type: ignore
    def mol_map(self, rules):
        return pm.mol_map(digit=(rules[1], rules[2], rules[3]))

    @_('":" digit digit')  # type: ignore
    def mol_map(self, rules):
        return pm.mol_map(digit=(rules[1], rules[2]))

    @_('":" digit')  # type: ignore
    def mol_map(self, rules):
        return pm.mol_map(digit=rules.digit)

    @_('"@"')  # type: ignore
    def chiral(self, rules):
        return pm.chiral(chiral="clockwise")

    @_('"@" "@"')  # type: ignore
    def chiral(self, rules):
        return pm.chiral(chiral="counterclockwise")

    @_("digit digit")  # type: ignore
    def fifteen(self, rules):
        return pm.fifteen(digit=(rules[0], rules[1]))

    @_("digit")  # type: ignore
    def fifteen(self, rules):
        return pm.fifteen(digit=rules.digit)

    @_("chains", "branch")
    def chain_branch_item(self, rules):
        return rules[0]


parser = SmilesParser()
lexer = SmilesLex()


def validate_smiles(mol: str) -> tuple[bool, Exception | None]:
    """
    Function for valdiating a SMILES molecule.

    Args:
        mol: Chemical formula as a string.
        use_only_grammar: For valdiating only the Grammar

    Returns:
        A Tuple containg in the first element if it is a valid SMILES and the second element a Exception.
    """
    try:
        parser.parse(lexer.tokenize(mol))
        pm._reset()
        return True, None
    except Exception as e:
        raise e
        return False, e
