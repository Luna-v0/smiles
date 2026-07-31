# Welcome

🚧 **This project is currently under construction.** 🚧

Welcome to the documentation for the **Smiles Validator** library, an essential part of a broader initiative focused on generating molecules using Generative AI techniques. This library is designed with two main components: the Validator and the Chemistry Module. The Validator is responsible for ensuring that SMILES (Simplified Molecular Input Line Entry System) strings are syntactically correct, using a robust parsing mechanism. Meanwhile, the Chemistry Module implements domain-specific chemical rules to guarantee the validity and consistency of molecular structures.

The SMILES validation system is powered by a custom parser built with the [SLY](https://github.com/dabeaz/sly) LALR(1) toolkit. The grammar used is adapted from the LL(1) parser described in [this article](https://depth-first.com/articles/2020/04/20/smiles-formal-grammar/) and the [OpenSMILES](https://opensmiles.org/opensmiles.html) specification. For chemical data, such as atomic weights and valence information, the project retrieves information from the [Bowserinator/Periodic-Table-JSON](https://github.com/Bowserinator/Periodic-Table-JSON/tree/master) repository.

## Entry point

The single public entry point is `validation.validate_smiles` (re-exported
from the package root and, for backwards compatibility, from `syntax.yacc`).
It returns the legacy `(is_valid, exception)` tuple;
`validation.validate_smiles_detailed` returns a structured
`ValidationResult` with the rejecting tier (`lex`, `grammar`,
`ring_semantics` or `chemistry`), a message and the failure position.
Internal errors — e.g. a missing chemistry backend — raise instead of being
reported as invalid molecules.

## Documentation

- [**Validation Algorithms**](validation.md) - Hückel's rule (aromaticity) and octet rule (valency) validation
- [**RDKit Comparison**](rdkit_comparison.md) - Comparison with RDKit validation

## Next Steps

- [x] Implementing aromaticity
    * [x] Check if rings are closed
    * [x] Check if all benzene examples are working properly
    * [x] Complex polycyclic aromatic hydrocarbons with multiple fused rings (Acridine, Carbazole, etc.)
    * [x] Molecules using advanced ring numbering patterns (reusing ring numbers)
    * [x] Molecules with mixed aromatic/non-aromatic structures
- [x] Bracket atom valency validation using octet rule
- [ ] Analysing and Implementing [RD filters](https://github.com/PatWalters/rd_filters)
- [ ] Analysing more chemistry filters for the parser like: [Filters 1](https://practicalcheminformatics.blogspot.com/2023/07/a-simple-tool-for-exploring-functional.html) [Filters 2](https://practicalcheminformatics.blogspot.com/2024/05/generative-molecular-design-isnt-as.html)
- [ ] Removing all shift/reduce conflicts (currently 26)
