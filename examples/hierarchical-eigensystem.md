# Hierarchical substitution model eigensystems

Set `eigenSystem` on `hierarchicalSubstitutionModel` to choose how it obtains the eigensystem for conventional BEAGLE pruning:

```xml
<hierarchicalSubstitutionModel id="hrate" eigenSystem="numerical">
    <!-- Existing frequencies, character tree, and rate parameters. -->
</hierarchicalSubstitutionModel>
```

The accepted values are `analytical` and `numerical`, without case sensitivity. Omitting the attribute selects `analytical`, preserving existing behavior. Invalid values cause an XML parse error. The model's startup message reports the selected eigensystem.

- `analytical` uses the hierarchical closed-form eigenvalues and eigenvectors.
- `numerical` constructs the same rate matrix Q and numerically diagonalizes its reversible symmetric transform with the existing `RobustEigenDecomposition` solver. It preserves the supplied frequencies, including exactly uniform frequencies, and supports repeated eigenvalues. Stationary frequencies must be finite and strictly positive for this mode.

The selection controls `getEigenDecomposition()`. Both modes retain the model's existing closed-form `getTransitionProbabilities()` calculation. Rates, frequencies, character-tree changes, and state restoration invalidate the cached eigensystem normally.

To compare numerical and analytical eigensystems with the BEAGLE SSE backend, use the ordinary representation in the enclosing likelihood:

```xml
<treeDataLikelihood id="likelihood" useHierarchicalRepresentation="false">
    <!-- Existing patterns, phylogeny, branch model, and site-rate model. -->
</treeDataLikelihood>
```

Run BEAST with the desired SSE/double-precision options, such as `-beagle_SSE -beagle_double -beagle_threads 0`, and verify the selected implementation in its log. Change only `eigenSystem` between the two runs. The XML fragments above show the attributes and must be combined with the normal model and likelihood children.

With `useHierarchicalRepresentation="true"`, the hrate plugin uses its hierarchical model representation directly and does not request the dense eigensystem. Its aggregation setting remains `aggregate="on|off|auto"`. The eigensystem selector does not select an execution backend or change aggregation.

For timing, distinguish cached pruning from rebuilding a model. Changing only phylogeny branch lengths reuses the eigensystem; changing a relative clade rate invalidates it. Numerical and analytical eigenvectors can differ within repeated-eigenvalue subspaces, so compare reconstructed Q, transition probabilities, and likelihoods rather than eigenvector entries.
