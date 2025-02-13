from sly import Parser
from lex import SmilesLex
import chem
import re
from itertools import combinations

# https://depth-first.com/articles/2020/04/20/smiles-formal-grammar/


def generate_combinations(rule: str):
    parts = rule.split()

    # Separate required and optional elements
    required = [p.rstrip('?') for p in parts if not p.endswith('?')]
    optional = [p.rstrip('?') for p in parts if p.endswith('?')]

    all_combinations = []

    # Generate all possible ordered combinations of optional elements
    for i in range(len(optional) + 1):
        for combo in combinations(optional, i):
            ordered_combo = [p.rstrip('?') for p in parts if p.rstrip(
                '?') in combo or p.rstrip('?') in required]
            all_combinations.append(" ".join(ordered_combo))

    return all_combinations


class SmilesParser(Parser):
    debugfile = 'parser.out'
    tokens = SmilesLex.tokens

    last_rule = ""

    def error(self, t):
        raise Exception(f"Error on {str(t)}, last rule {self.last_rule}")

    def update_last_rule(self, rule_name):
        self.last_rule = rule_name

    @_('atom', 'atom chain_branch')  # type: ignore
    def line(self, rules):
        self.update_last_rule('line')
        pass

    @_('chains', 'branch', 'chain_branch chains', 'chain_branch branch')  # type: ignore
    def chain_branch(self, rules):
        self.update_last_rule('chain_branch')
        pass

    @_('chain', 'chain chains')  # type: ignore
    def chains(self, rules):
        self.update_last_rule('chains')
        pass

    @_('"[" internal_bracket "]"')  # type: ignore
    def bracket_atom(self, rules):
        self.update_last_rule('bracket_atom')
        pass

    @_(*generate_combinations('isotope? hcount? symbol chiral? charge? map?')) # type: ignore
    def internal_bracket(self, rules):
        chem.validate_valency_mol(rules.isotope, rules.symbol,
                                  rules.chiral, rules.hcount,
                                  rules.charge, rules.map)
        pass

    @_('dot_proxy', 'bond atom', 'bond rnum', 'atom', 'rnum')  # type: ignore
    def chain(self, rules):
        self.update_last_rule('chain')
        pass

    @_('"." atom')  # type: ignore
    def dot_proxy(self, rules):
        pass

    @_('semi_symbol', 'organic_symbol')  # type: ignore
    def symbol(self, rules):
        self.update_last_rule('symbol')
        return rules[0] if len(rules) > 0 else None

    @_('"(" inner_branch ")"')  # type: ignore
    def branch(self, rules):
        self.update_last_rule('branch')
        pass

    @_('bond_dot line', 'line', 'inner_branch bond_dot line', 'inner_branch line') # type: ignore
    def inner_branch(self, rules):
        self.update_last_rule('inner_branch')
        pass

    @_('bond', '"."')  # type: ignore
    def bond_dot(self, rules):
        self.update_last_rule('bond_dot')
        return rules[0] if len(rules) > 0 else None

    @_('semi_bond_rule', '"-"')  # type: ignore
    def bond(self, rules):
        self.update_last_rule('bond')
        return rules[0] if len(rules) > 0 else None

    @_('semi_bond')  # type: ignore
    def semi_bond_rule(self, rules):
        self.update_last_rule('semi_bond_rule')
        pass

    @_('"H"', 'semi_organic_rule')  # type: ignore
    def organic_symbol(self, rules):
        self.update_last_rule('organic_symbol')
        pass

    @_('semi_organic_symbol')  # type: ignore
    def semi_organic_rule(self, rules):
        pass

    @_('organic_symbol', 'bracket_atom')  # type: ignore
    def atom(self, rules):
        self.update_last_rule('atom')
        pass

    @_('digit', '"%" digit digit ')  # type: ignore
    def rnum(self, rules):
        self.update_last_rule('rnum')
        pass

    @_('digit digit digit', 'digit digit', 'digit')  # type: ignore
    def isotope(self, rules):
        self.update_last_rule('isotope')
        pass

    @_('"H" digit', '"H"')  # type: ignore
    def hcount(self, rules):
        self.update_last_rule('hcount')
        pass

    @_('"+"', '"+" "+"', '"-"', '"-" "-"', '"-" fifteen', '"+" fifteen')  # type: ignore
    def charge(self, rules):
        self.update_last_rule('charge')
        pass

    @_('":" digit digit digit', '":" digit digit', '":" digit')  # type: ignore
    def map(self, rules):
        self.update_last_rule('map')
        pass

    @_('"@"', '"@" "@"')  # type: ignore
    def chiral(self, rules):
        self.update_last_rule('chiral')
        pass

    @_('digit digit', 'digit')  # type: ignore
    def fifteen(self, rules):
        self.update_last_rule('fifteen')
        pass


parser = SmilesParser()
lexer = SmilesLex()


def validateSMILES(mol: str) -> bool:
    try:
        parser.parse(lexer.tokenize(mol))
        return True
    except Exception as e:
        print(e)
        return False


print(validateSMILES("NCCNCCNCCN"))
