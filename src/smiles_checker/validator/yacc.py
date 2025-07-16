from itertools import combinations

from sly import Parser

from smiles_checker.chem.chemistry import chemistry as chem
from smiles_checker.validator.lex import SmilesLex
from smiles_checker.validator.parser_manager import parser_manager as pm


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


def dictify(obj):
    result = {}
    if hasattr(obj, "_namemap") and hasattr(obj, "_slice"):
        print(obj._namemap)
        for name, getter in obj._namemap.items():
            try:
                value = getter(obj._slice)
            except Exception as e:
                value = f"<error: {e}>"
            result[name] = value
    return result


class SmilesParser(Parser):
    """
    Parser using the SLY library to parse SMILES strings.
    """

    debugfile = "parser.out"
    tokens = SmilesLex.tokens
    use_only_grammar = False

    def error(self, t):
        raise Exception(f"Error on {str(t)}")

    @_("atom chain_branch")  # type: ignore
    def line(self, rules):
        return pm.line(atom=rules.atom, chain_branch=rules.chain_branch)

    @_("atom")
    def line(self, rules):
        return pm.line(atom=rules.atom)

    @_("chains chain_branch")  # type: ignore
    def chain_branch(self, rules):
        return pm.chain_branch(chains=rules.chains, branch=rules.chain_branch)

    @_("branch chain_branch")  # type: ignore
    def chain_branch(self, rules):
        return pm.chain_branch(branch=rules.branch, chains=rules.chain_branch)

    @_("branch")
    def chain_branch(self, rules):
        return pm.chain_branch(branch=rules.branch)

    @_("chains")
    def chain_branch(self, rules):
        return pm.chain_branch(chains=rules.chains)

    @_("chain")  # type: ignore
    def chains(self, rules):
        return pm.chains(chain=rules.chain)

    @_("chain chains")  # type: ignore
    def chains(self, rules):
        return pm.chains(chain=rules.chain, chains=rules.chains)

    @_('"[" internal_bracket "]"')  # type: ignore
    def bracket_atom(self, rules):
        return rules.internal_bracket

    @_(*generate_combinations("isotope? symbol chiral? hcount? charge? mol_map?"))  # type: ignore
    def internal_bracket(self, rules):
        return pm.internal_bracket(
            *getAttributes(
                rules, ["isotope", "symbol", "chiral", "hcount", "charge", "mol_map"]
            )
        )

    @_("dot_proxy")
    def chain(self, rules):
        return pm.chain(dot_proxy=rules.dot_proxy)

    @_("bond atom")
    def chain(self, rules):
        return pm.chain(bond=rules.bond, atom=rules.atom)

    @_("bond rnum")
    def chain(self, rules):
        return pm.chain(bond=rules.bond, rnum=rules.rnum)

    @_("atom")
    def chain(self, rules):
        return pm.chain(atom=rules.atom)

    @_("rnum")  # type: ignore
    def chain(self, rules):
        return pm.chain(rnum=rules.rnum)

    @_('"." atom')  # type: ignore
    def dot_proxy(self, rules):
        return pm.dot_proxy(rules.atom)

    @_("semi_symbol", '"H"')  # type: ignore
    def symbol(self, rules):
        return rules[0]

    @_('"(" inner_branch ")"')  # type: ignore
    def branch(self, rules):
        return pm.branch(rules.inner_branch)

    @_("bond_dot line")  # type: ignore
    def inner_branch(self, rules):
        return pm.inner_branch(bond_dot=rules.bond_dot, line=rules.line)

    @_("line")
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
        return rules[0]

    @_("semi_bond_rule", '"-"')  # type: ignore
    def bond(self, rules):
        return rules[0]

    @_("semi_bond")  # type: ignore
    def semi_bond_rule(self, rules):
        return rules[0]

    @_("symbol")  # type: ignore
    def atom(self, rules):
        return pm.atom(symbol=rules.symbol)

    @_("bracket_atom")
    def atom(self, rules):
        return pm.atom(bracket_atom=rules.bracket_atom)

    @_("digit")  # type: ignore
    def rnum(self, rules):
        return pm.rnum(ring_number=rules.digit)

    @_('"%" digit digit')
    def rnum(self, rules):
        return pm.rnum(ring_number=10 * rules[1] + rules[2])

    @_("digit")  # type: ignore
    def isotope(self, rules):
        return pm.isotope(value=rules.digit)

    @_("digit digit")
    def isotope(self, rules):
        return pm.isotope(value=10 * rules[0] + rules[1])

    @_("digit digit digit")  # type: ignore
    def isotope(self, rules):
        return pm.isotope(value=100 * rules[0] + 10 * rules[1] + rules[2])

    @_('"H"')  # type: ignore
    def hcount(self, rules):
        return pm.hcount(hcount=1)

    @_('"H" digit')
    def hcount(self, rules):
        return pm.hcount(hcount=rules.digit)

    @_('"+"')
    def charge(self, rules):
        return pm.charge(charge=1)

    @_('"-"')
    def charge(self, rules):
        return pm.charge(charge=-1)

    @_('"+" fifteen')
    def charge(self, rules):
        return pm.charge(charge=rules.fifteen)

    @_('"-" fifteen')
    def charge(self, rules):
        return pm.charge(charge=-rules.fifteen)

    @_('"+" "+"')
    def charge(self, rules):
        return pm.charge(charge=2)

    @_('"-" "-"')
    def charge(self, rules):
        return pm.charge(charge=-2)

    @_('":" digit digit digit')
    def mol_map(self, rules):
        return pm.mol_map(value=100 * rules[1] + 10 * rules[2] + rules[3])

    @_('":" digit digit')
    def mol_map(self, rules):
        return pm.mol_map(value=10 * rules[1] + rules[2])

    @_('":" digit')
    def mol_map(self, rules):
        return pm.mol_map(value=rules.digit)

    @_('"@"')
    def chiral(self, rules):
        return pm.chiral(rotation="clockwise")

    @_('"@" "@"')
    def chiral(self, rules):
        return pm.chiral(rotation="counterclockwise")

    @_("digit digit")  # type: ignore
    def fifteen(self, rules):
        return pm.fifteen(value=10 * rules[0] + rules[1])

    @_("digit")
    def fifteen(self, rules):
        return pm.fifteen(value=rules.digit)


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
        pm.validate_branch()
        pm._reset()
        return True, None
    except Exception as e:
        raise e
        return False, e
