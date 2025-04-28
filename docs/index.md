# Welcome

🚧 **This project is currently under construction.** 🚧

Welcome to the documentation for the **Smiles Validator** library, an essential part of a broader initiative focused on generating molecules using Generative AI techniques. This library is designed with two main components: the Validator and the Chemistry Module. The Validator is responsible for ensuring that SMILES (Simplified Molecular Input Line Entry System) strings are syntactically correct, using a robust parsing mechanism. Meanwhile, the Chemistry Module implements domain-specific chemical rules to guarantee the validity and consistency of molecular structures.

The SMILES validation system is powered by a custom parser built with the [SLY](https://github.com/dabeaz/sly) LALR(1) toolkit. The grammar used is adapted from the LL(1) parser described in [this article](https://depth-first.com/articles/2020/04/20/smiles-formal-grammar/) and the [OpenSMILES](https://opensmiles.org/opensmiles.html) specification. For chemical data, such as atomic weights and valence information, the project retrieves information from the [Bowserinator/Periodic-Table-JSON](https://github.com/Bowserinator/Periodic-Table-JSON/tree/master) repository.

## Next Steps

- [ ] Implementing aromacity
    * [x] Check if rings are closed
    * [ ] Check if all benzene examples are working properly
    * [ ] Finish passing everything from the yacc to the parser manager
- [ ] Analysing and Implementing [RD filters](https://github.com/PatWalters/rd_filters)
- [ ] Analysing more chemistry filters for the parser like: [Filters 1](https://practicalcheminformatics.blogspot.com/2023/07/a-simple-tool-for-exploring-functional.html) [Filters 2](https://practicalcheminformatics.blogspot.com/2024/05/generative-molecular-design-isnt-as.html)
- [ ] Removing all shift/reduce conflicts
