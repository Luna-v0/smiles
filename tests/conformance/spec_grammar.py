"""Reference acceptor for the OpenSMILES formal grammar (§2.2).

A direct recursive-descent transcription of the BNF at
http://opensmiles.org/opensmiles.html — deliberately independent of the SLY
parser under test, so the two can be compared differentially.

This is a *syntax* acceptor: it decides derivability from the context-free
grammar only.  Ring-closure semantics (§3.4/§3.6) and chemistry are out of
scope here, which is exactly what makes it usable as the soundness oracle:
everything the validator accepts must at minimum be derivable from this
grammar (the validator's extra semantic tiers only ever *shrink* the
language).

Documented extensions
---------------------
With ``extensions=True`` (the default used by the differential tests) the
acceptor also admits the validator's known, intentional deviations from the
strict BNF:

- **E1 — bare atoms beyond the organic subset**: any single-character
  periodic-table symbol (e.g. ``K``, ``W``, ``U``, ``H``) plus aromatic
  ``se``/``as`` may appear outside brackets.  Strict BNF allows only
  ``B C N O S P F Cl Br I`` and ``b c n o p s``.
- **E2 — multi-digit hcount**: ``[CH12]``.  Strict BNF is ``'H' DIGIT?``.
- **E3 — charge digit width**: any digit run whose value is at most 14
  (e.g. ``[C+007]``).  Strict BNF is ``sign DIGIT? DIGIT``; the validator
  additionally caps the magnitude below 15.
- **E4 — ``@@`` with a chirality class**: ``[C@@TH2]``.  Strict BNF
  enumerates ``@`` ``@@`` ``@TH1`` … as alternatives, so ``@@TH2`` is not
  derivable; accepting it is mandated by the conformance plan (RDKit
  rejects it — an intentional disagreement).

Known restrictions (the validator rejects; strict BNF derives) are listed in
``test_differential.py`` — a restriction can never break soundness, only
completeness, and the generator avoids producing them.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src'))
from chem.chemistry import pt_symbols

ALIPHATIC_ORGANIC = ("Cl", "Br", "B", "C", "N", "O", "S", "P", "F", "I")
AROMATIC_ORGANIC = ("b", "c", "n", "o", "p", "s")
AROMATIC_SYMBOLS = ("se", "as", "b", "c", "n", "o", "p", "s")
BOND_CHARS = "-=#$:/\\"

# Full periodic table (case-sensitive), longest first for greedy matching.
ELEMENT_SYMBOLS = sorted(set(pt_symbols), key=len, reverse=True)
SINGLE_CHAR_ELEMENTS = tuple(s for s in pt_symbols if len(s) == 1)

CHIRALITY_RANGES = {"TH": 2, "AL": 2, "SP": 3, "TB": 20, "OH": 30}


class _Reject(Exception):
    """Internal signal: the input is not derivable."""


class SpecSmilesAcceptor:
    """Recursive-descent acceptor for the OpenSMILES grammar.

    Args:
        extensions: Also accept the documented extensions E1–E4 above.
    """

    def __init__(self, extensions: bool = True):
        self.extensions = extensions
        self.text = ""
        self.pos = 0

    # -- character helpers -------------------------------------------------
    def _peek(self) -> str:
        return self.text[self.pos] if self.pos < len(self.text) else ""

    def _take(self, s: str) -> bool:
        if self.text.startswith(s, self.pos):
            self.pos += len(s)
            return True
        return False

    def _expect(self, s: str) -> None:
        if not self._take(s):
            raise _Reject

    def _digits(self, min_len: int = 1, max_len: int | None = None) -> str:
        start = self.pos
        while self._peek().isdigit() and (max_len is None or self.pos - start < max_len):
            self.pos += 1
        if self.pos - start < min_len:
            raise _Reject
        return self.text[start:self.pos]

    # -- grammar rules -----------------------------------------------------
    def accepts(self, smiles: str) -> bool:
        """True iff ``smiles`` is derivable from the (extended) grammar."""
        self.text = smiles.strip()
        self.pos = 0
        if not self.text:
            return False
        try:
            self._chain()
        except _Reject:
            return False
        return self.pos == len(self.text)

    def _chain(self) -> None:
        # chain ::= branched_atom | chain branched_atom
        #         | chain bond branched_atom | chain dot branched_atom
        self._branched_atom()
        while True:
            save = self.pos
            try:
                if self._peek() == ".":
                    self.pos += 1
                elif self._peek() in BOND_CHARS:
                    self.pos += 1
                self._branched_atom()
            except _Reject:
                self.pos = save
                return

    def _branched_atom(self) -> None:
        # branched_atom ::= atom ringbond* branch*
        self._atom()
        while self._ringbond():
            pass
        while self._branch():
            pass

    def _ringbond(self) -> bool:
        # ringbond ::= bond? DIGIT | bond? '%' DIGIT DIGIT
        save = self.pos
        if self._peek() in BOND_CHARS:
            self.pos += 1
        if self._peek() == "%":
            self.pos += 1
            try:
                self._digits(2, 2)
                return True
            except _Reject:
                self.pos = save
                return False
        if self._peek().isdigit():
            self.pos += 1
            return True
        self.pos = save
        return False

    def _branch(self) -> bool:
        # branch ::= '(' chain ')' | '(' bond chain ')' | '(' dot chain ')'
        if self._peek() != "(":
            return False
        self.pos += 1
        if self._peek() in BOND_CHARS or self._peek() == ".":
            self.pos += 1
        self._chain()
        self._expect(")")
        return True

    def _atom(self) -> None:
        # atom ::= bracket_atom | aliphatic_organic | aromatic_organic | '*'
        if self._peek() == "[":
            self._bracket_atom()
            return
        if self._take("*"):
            return
        for sym in ALIPHATIC_ORGANIC + AROMATIC_ORGANIC:
            if self._take(sym):
                return
        if self.extensions:
            # E1: bare se/as and any single-character element symbol.
            for sym in ("se", "as") + SINGLE_CHAR_ELEMENTS + ("H",):
                if self._take(sym):
                    return
        raise _Reject

    def _bracket_atom(self) -> None:
        # bracket_atom ::= '[' isotope? symbol chiral? hcount? charge? class? ']'
        self._expect("[")
        if self._peek().isdigit():
            self._digits()  # isotope ::= NUMBER
        self._symbol()
        self._chiral()
        self._hcount()
        self._charge()
        if self._take(":"):
            self._digits()  # class ::= ':' NUMBER
        self._expect("]")

    def _symbol(self) -> None:
        # symbol ::= element_symbols | aromatic_symbols | '*'
        if self._take("*"):
            return
        for sym in AROMATIC_SYMBOLS:
            if self._take(sym):
                return
        for sym in ELEMENT_SYMBOLS:
            if self._take(sym):
                return
        raise _Reject

    def _chiral(self) -> None:
        # chiral ::= '@' | '@@' | '@TH1' … '@OH30'  (+E4: '@@' before a class)
        if self._peek() != "@":
            return
        self.pos += 1
        double = self._take("@")
        for keyword, limit in CHIRALITY_RANGES.items():
            if double and not self.extensions:
                break  # strict BNF: '@@' never precedes a class
            save = self.pos
            if self._take(keyword):
                try:
                    value = int(self._digits(1, 2))
                except _Reject:
                    self.pos = save
                    raise
                if not 1 <= value <= limit:
                    raise _Reject
                return

    def _hcount(self) -> None:
        # hcount ::= 'H' | 'H' DIGIT  (+E2: multi-digit)
        if not self._take("H"):
            return
        if self._peek().isdigit():
            self._digits(1, None if self.extensions else 1)

    def _charge(self) -> None:
        # charge ::= '-' | '+' | '-' DIGIT? DIGIT | '+' DIGIT? DIGIT
        #          | '--' | '++'  (deprecated)  (+E3: wider digit runs, value<15)
        sign = self._peek()
        if sign not in "+-":
            return
        self.pos += 1
        if self._take(sign):  # '--' / '++'
            return
        if self._peek().isdigit():
            digits = self._digits(1, None if self.extensions else 2)
            if self.extensions and int(digits) >= 15:
                raise _Reject
