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

```py
from validator.yacc import validate_smiles

print(validate_smiles("c1ccccc1"))
```

## 🧪 Run Tests

```py
pytest
```
