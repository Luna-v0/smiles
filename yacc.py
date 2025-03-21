from sly import Parser
from lex import SmilesLex
import chem as chem
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
    use_only_grammar = False

    def error(self, t):
        raise Exception(f"Error on {str(t)}")

    @_('atom', 'atom chain_branch')  # type: ignore
    def line(self, rules):
        pass

    @_('chains', 'branch', 'chain_branch chains', 'chain_branch branch')  # type: ignore
    def chain_branch(self, rules):
        pass

    @_('chain', 'chain chains')  # type: ignore
    def chains(self, rules):
        pass

    @_('"[" internal_bracket "]"')  # type: ignore
    def bracket_atom(self, rules):
        return rules.internal_bracket


    @_(*generate_combinations('isotope? symbol chiral? hcount? charge? map?')) # type: ignore
    def internal_bracket(self, rules):
        
        mol = getAttributes(rules, ['isotope', 'symbol', 'chiral','hcount', 'charge', 'map'])
        
        isotope, symbol, chiral, hcount, charge, mol_map = mol

        if self.use_only_grammar: return

        if not chem.validate_valency_bracket(isotope, symbol, chiral, hcount, charge, mol_map): 
            raise Exception(f'Invalid valency in Bracket [{','.join([str(x) for x in mol if x is not None])}]')
        

    @_('dot_proxy', 'bond atom', 'bond rnum', 'atom', 'rnum')  # type: ignore
    def chain(self, rules):
        pass

    @_('"." atom')  # type: ignore
    def dot_proxy(self, rules):
        pass

    @_('semi_symbol', '"H"')  # type: ignore
    def symbol(self, rules):
        return rules[0]

    @_('"(" inner_branch ")"')  # type: ignore
    def branch(self, rules):
        pass

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
        if type(rules[0]) != str or len(rules[0]) == 1 or rules[0] in chem.organic_atoms: return rules[0]
        
        elem1, elem2 = rules[0]
        
        if elem1 in chem.organic_atoms and elem2 in chem.organic_atoms:
            return [elem1,elem2]
        
        raise Exception(f"Inorganic Atom Outside Bracket {rules[0]}")

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


def validateSMILES(mol: str,use_only_grammar=False) -> bool:
    try:
        parser.use_only_grammar = use_only_grammar
        parser.parse(lexer.tokenize(mol))
        return True
    except Exception as e:
        print(type(e))
        print("Error:",e)
        return False

