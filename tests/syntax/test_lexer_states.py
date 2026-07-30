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
        assert types == ["semi_symbol", "[", "semi_symbol", "H", "digit", "]", "semi_symbol"]

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
