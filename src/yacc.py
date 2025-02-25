from sly import Parser
from src.lex import SmilesLex
import src.chem as chem
import re
from itertools import combinations
import json  # temporary

# https://depth-first.com/articles/2020/04/20/smiles-formal-grammar/


def generate_combinations(rule: str):
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
    if type(properties) != list:
        return getattr(rules, properties, None)

    values = []

    for prop in properties:
        values.append(getAttributes(rules, prop))

    return values


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
        return rules.internal_bracket


    @_(*generate_combinations('isotope? hcount? symbol chiral? charge? map?')) # type: ignore
    def internal_bracket(self, rules):
        
        mol = getAttributes(rules, ['isotope', 'hcount', 'symbol', 'chiral', 'charge', 'map'])
        
        isotope, hcount, symbol, chiral, charge, mol_map = mol

        print(f"Isotope: {isotope}")
        print(f"Hcount: {hcount}")
        print(f"Symbol: {symbol}")
        print(f"Chiral: {chiral}")
        print(f"Charge: {charge}")
        print(f"Map: {mol_map}")

        chem.validate_valency_mol(isotope, symbol, chiral, hcount, charge, mol_map)
        pass

    @_('dot_proxy', 'bond atom', 'bond rnum', 'atom', 'rnum')  # type: ignore
    def chain(self, rules):

        pass

    @_('"." atom')  # type: ignore
    def dot_proxy(self, rules):
        pass

    @_('semi_symbol', 'organic_symbol')  # type: ignore
    def symbol(self, rules):
        return rules[0]

    @_('"(" inner_branch ")"')  # type: ignore
    def branch(self, rules):
        self.update_last_rule('branch')
        pass


    @_('bond_dot line', 'line', 'bond_dot line inner_branch', 'line inner_branch')    # type: ignore
    def inner_branch(self, rules):
        self.update_last_rule('inner_branch')
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

    @_('"H"', 'semi_organic_rule')  # type: ignore
    def organic_symbol(self, rules):
        return rules[0]

    @_('semi_organic_symbol')  # type: ignore
    def semi_organic_rule(self, rules):
        return rules[0]

    @_('organic_symbol', 'bracket_atom')  # type: ignore
    def atom(self, rules):
        return rules[0]

    @_('digit', '"%" digit digit ')  # type: ignore
    def rnum(self, rules):
        if rules[0] == '%':
            return int(''.join(rules[1:]))

        return int(rules[0])

    @_('digit digit digit', 'digit digit', 'digit')  # type: ignore
    def isotope(self, rules):
        return int(''.join(rules))

    @_('"H" digit', '"H"')  # type: ignore
    def hcount(self, rules):
        if len(rules) == 1:
            return 1

        return int(rules.digit)

    @_('"+"', '"+" "+"', '"-"', '"-" "-"', '"-" fifteen', '"+" fifteen')  # type: ignore
    def charge(self, rules):
        if len(rules) == 1:
            return 1 if rules[0] == "+" else -1

        if rules[1] == '-':
            return -2

        if rules[1] == '+':
            return 2

        if rules[0] == '-':
            return rules.fifteen * -1

        return rules.fifteen

    @_('":" digit digit digit', '":" digit digit', '":" digit')  # type: ignore
    def map(self, rules):

        digits = rules[1:]

        return int(''.join(digits))

    @_('"@"', '"@" "@"')  # type: ignore
    def chiral(self, rules):
        self.update_last_rule('chiral')
        pass

    @_('digit digit', 'digit')  # type: ignore
    def fifteen(self, rules):

        if len(rules) == 2:
            x = int(rules[0]+rules[1])

            if x > 15:
                raise "cannot execeed 15"

            return x

        return int(rules.digit)


parser = SmilesParser()
lexer = SmilesLex()


def validateSMILES(mol: str) -> bool:
    try:
        parser.parse(lexer.tokenize(mol))
        return True
    except Exception as e:
        return False


print(validateSMILES("[Cl+]"))
