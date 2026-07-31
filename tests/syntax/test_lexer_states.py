"""Regression tests for the bracket-state lexer (Phase 1).

The lexer enters ``BracketLex`` on ``[`` and returns to ``SmilesLex`` on
``]``.  Bracket contents must tokenize with their own rules (full periodic
table, no bare-organic-subset splitting) and an aborted parse inside a
bracket must not poison the next tokenization.
"""

import pytest

from syntax.lex import BracketLex, SmilesLex


@pytest.fixture
def lexer():
    return SmilesLex()


def _types_values(lexer, text):
    return [(t.type, t.value) for t in lexer.tokenize(text)]


class TestBracketState:
    def test_two_char_symbol_inside_bracket_not_split(self, lexer):
        toks = _types_values(lexer, "[Cn]")
        assert ("semi_symbol", "Cn") in toks

    def test_two_char_symbol_outside_bracket_is_split(self, lexer):
        toks = _types_values(lexer, "Cn1ccnc1")
        assert toks[0] == ("semi_symbol", "C")
        assert toks[1] == ("semi_symbol", "n")

    def test_brackets_delimit_state(self, lexer):
        types = [t for t, _ in _types_values(lexer, "C[CH3]C")]
        assert types == ["semi_symbol", "[", "semi_symbol", "H", "NUMBER", "]", "semi_symbol"]

    def test_bracket_digit_run_is_one_number(self, lexer):
        toks = _types_values(lexer, "[13CH4:1234]")
        assert ("NUMBER", 13) in toks
        assert ("NUMBER", 1234) in toks

    def test_chirality_is_single_token(self, lexer):
        toks = _types_values(lexer, "[C@TH1][C@@][S@AL2]")
        chirals = [v for t, v in toks if t == "CHIRAL"]
        assert chirals == ["@TH1", "@@", "@AL2"]

    def test_hydroxide_not_eaten_by_chirality_keyword(self, lexer):
        # OH is only a chirality class after '@'; [OH-] must stay O + H + '-'
        types = [t for t, _ in _types_values(lexer, "[OH-]")]
        assert types == ["[", "semi_symbol", "H", "-", "]"]

    def test_state_restored_after_bracket(self, lexer):
        # The trailing Sc must be split (outside brackets only Cl/Br are bare)
        toks = _types_values(lexer, "[Sc]Sc")
        assert toks[0:1] == [("[", "[")]
        assert ("semi_symbol", "Sc") in toks[:3]
        assert toks[-2:] == [("semi_symbol", "S"), ("semi_symbol", "c")]

    def test_aborted_bracket_does_not_poison_next_call(self, lexer):
        # Simulate a parse that dies while the lexer is inside a bracket
        with pytest.raises(Exception):
            for _ in lexer.tokenize("[C"):
                raise RuntimeError("abort mid-bracket")
        toks = _types_values(lexer, "CCO")
        assert [t for t, _ in toks] == ["semi_symbol"] * 3
