# Smiles Validator

A validator for SMILE chemical language using SLY LaLR(1) parser. The grammar was based on the LL(1) parser from this [article](https://depth-first.com/articles/2020/04/20/smiles-formal-grammar/), and from the [OpenSMILES](https://opensmiles.org/opensmiles.html) Specification.

All the data related to the periodic table was retrived from the [Bowserinator/Periodic-Table-JSON repo](https://github.com/Bowserinator/Periodic-Table-JSON/tree/master).

## 📚 Documenatation

All the documentation is available in the [Docs Page](https://luna-v0.github.io/smiles/).

## 📦 Installation

For installation make sure you have [uv](https://github.com/astral-sh/uv) installed.

```bash
uv venv .venv
source .venv/bin/activate
uv pip install -e .
```

## 🚀 Usage

There is a single validation entry point, `validate_smiles`, which runs the
full pipeline: lexing, grammar, ring-closure semantics and chemistry
validation (valence + aromaticity). It is available both from the package
root and — for backwards compatibility — from `syntax.yacc`; the two are
the same function.

```py
from src import validate_smiles

is_valid, error = validate_smiles("c1ccccc1")  # (True, None)
```

For structured failure information (which tier rejected the input, at which
character position), use `validate_smiles_detailed`:

```py
from src import validate_smiles_detailed

result = validate_smiles_detailed("C1CCC")
result.valid     # False
result.tier      # "ring_semantics"
result.message   # "Unclosed ring numbers: 1"
result.as_tuple()  # legacy (is_valid, exception) contract
```

Internal errors (a missing chemistry backend, a non-string input) **raise**
instead of being reported as "invalid molecule".

## 🧪 Run Tests

```py
pytest
```
