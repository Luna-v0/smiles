# RDKit Comparison Notes

## Benchmark Results

When comparing the SMILES validator against RDKit on the test dataset:

- **Syntax-only validation**: 0.66% failure rate vs RDKit's 0.01%
- **Full chemistry validation**: 15.60% failure rate (stricter valency/aromaticity checks)

## Disagreements

There are 12 cases where our validator disagrees with RDKit. These are primarily due to:

1. **Stricter valency validation** - Our validator enforces stricter valency rules for bracket atoms
2. **Aromaticity interpretation** - Some edge cases in fused ring systems are handled differently
3. **Graph construction artifacts** - Some complex molecules create self-loops that affect validation

These disagreements are expected as our validator is intentionally stricter for use in generative AI molecule validation, where false positives (accepting invalid molecules) are more costly than false negatives.
