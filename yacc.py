from sly import Parser
from lex import SmilesLex
import chem

# https://depth-first.com/articles/2020/04/20/smiles-formal-grammar/

class SmilesParser(Parser):
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

    @_('chain', 'chains chain')  # type: ignore
    def chains(self, rules):
        self.update_last_rule('chains')
        pass

    @_('"[" opt_isotope symbol opt_chiral opt_hcount opt_charge opt_map "]"', #type: ignore
       '"[" opt_isotope opt_hcount symbol opt_chiral opt_charge opt_map "]"')
    def bracket_atom(self, rules):
        self.update_last_rule('bracket_atom')
        chem.validate_valency_mol(rules.opt_isotope, rules.symbol,
                                  rules.opt_chiral, rules.opt_hcount,
                                  rules.opt_charge, rules.opt_map)
        pass

    @_('"." atom', 'opt_bond atom', 'opt_bond rnum')  # type: ignore
    def chain(self, rules):
        self.update_last_rule('chain')
        pass

    @_('semi_symbol', 'organic_symbol')  # type: ignore
    def symbol(self, rules):
        self.update_last_rule('symbol')
        return rules[0] if len(rules) > 0 else None

    @_('"(" inner_branch ")"')  # type: ignore
    def branch(self, rules):
        self.update_last_rule('branch')
        pass

    @_('opt_bond_dot line', 'inner_branch opt_bond_dot line')  # type: ignore
    def inner_branch(self, rules):
        self.update_last_rule('inner_branch')
        pass

    @_('bond', 'empty')  # type: ignore
    def opt_bond(self, rules):
        self.update_last_rule('opt_bond')
        return rules[0] if len(rules) > 0 else None

    @_('isotope', 'empty')  # type: ignore
    def opt_isotope(self, rules):
        self.update_last_rule('opt_isotope')
        return rules[0] if len(rules) > 0 else None

    @_('chiral', 'empty')  # type: ignore
    def opt_chiral(self, rules):
        self.update_last_rule('opt_chiral')
        return rules[0] if len(rules) > 0 else None

    @_('hcount', 'empty')  # type: ignore
    def opt_hcount(self, rules):
        self.update_last_rule('opt_hcount')
        return rules[0] if len(rules) > 0 else None

    @_('bond', '"."', 'empty')  # type: ignore
    def opt_bond_dot(self, rules):
        self.update_last_rule('opt_bond_dot')
        return rules[0] if len(rules) > 0 else None

    @_('charge', 'empty')  # type: ignore
    def opt_charge(self, rules):
        self.update_last_rule('opt_charge')
        return rules[0] if len(rules) > 0 else None

    @_('map', 'empty')  # type: ignore
    def opt_map(self, rules):
        self.update_last_rule('opt_map')
        return rules[0] if len(rules) > 0 else None

    @_('digit', 'empty')  # type: ignore
    def opt_digit(self, rules):
        self.update_last_rule('opt_digit')
        return rules[0] if len(rules) > 0 else None

    @_('semi_bond', '"-"')  # type: ignore
    def bond(self, rules):
        self.update_last_rule('bond')
        return rules[0] if len(rules) > 0 else None

    @_('"H"', 'semi_organic_symbol')  # type: ignore
    def organic_symbol(self, rules):
        self.update_last_rule('organic_symbol')
        pass

    @_('organic_symbol', 'bracket_atom')  # type: ignore
    def atom(self, rules):
        self.update_last_rule('atom')
        pass

    @_('digit', '"%" digit digit ')  # type: ignore
    def rnum(self, rules):
        self.update_last_rule('rnum')
        pass

    @_('opt_digit opt_digit digit')  # type: ignore
    def isotope(self, rules):
        self.update_last_rule('isotope')
        pass

    @_('"H" opt_digit')  # type: ignore
    def hcount(self, rules):
        self.update_last_rule('hcount')
        pass

    @_('"+"', '"+" "+"', '"-"', '"-" "-"', '"-" fifteen', '"+" fifteen')  # type: ignore
    def charge(self, rules):
        self.update_last_rule('charge')
        pass

    @_('":" opt_digit opt_digit digit')  # type: ignore
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
    
    @_('')  # type: ignore
    def empty(self,rules):
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