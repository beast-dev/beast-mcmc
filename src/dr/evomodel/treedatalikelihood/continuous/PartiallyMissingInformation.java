/*
 * PartiallyMissingInformation.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.Tree;
import dr.math.matrixAlgebra.Vector;

import java.util.Arrays;
import java.util.List;

/**
 * @author Marc A. Suchard
 */

public class PartiallyMissingInformation {

    public PartiallyMissingInformation(Tree tree, ContinuousTraitPartialsProvider dataModel,
                                       ContinuousDataLikelihoodDelegate likelihoodDelegate) {
        this.tipCount = tree.getExternalNodeCount();
        this.numTraits = dataModel.getTraitCount(); //likelihoodDelegate.getTraitCount();
        this.dimTrait = dataModel.getTraitDimension(); //likelihoodDelegate.getTraitDim();

        this.rawMissingIndices = dataModel.getMissingIndices();

        final int length = tipCount * numTraits;
        anyMissing = new boolean[length];
        allMissing = new boolean[length];
        missingIndices = new HashedIntArray[length];

        setupIndices();
    }

    public boolean isPartiallyMissing(final int tip, final int trait) {
        return anyMissing[getIndex(tip, trait)];
    }

    public boolean isCompletelyMissing(final int tip, final int trait) {
        return allMissing[getIndex(tip, trait)];
    }

    public HashedIntArray getMissingIndices(final int tip, final int trait) {
        return missingIndices[getIndex(tip, trait)];
    }

    private void setupIndices() {

        for (int taxon = 0; taxon < tipCount; ++taxon) {

            for (int trait = 0; trait < numTraits; ++trait) {

                final int index = getIndex(taxon, trait);

                int count = 0;
                for (int dim = 0; dim < dimTrait; ++dim) {
                    if (isObservationMissing(index, dim)) {
                        ++count;
                    }
                }

                if (count > 0) {

                    int[] missing = new int[count];
                    int[] notMissing = new int[dimTrait - count];

                    int offsetMissing = 0;
                    int offsetNotMissing = 0;
                    for (int dim = 0; dim < dimTrait; ++dim) {
                        if (isObservationMissing(index, dim)) {
                            missing[offsetMissing] = dim;
                            ++offsetMissing;
                        } else {
                            notMissing[offsetNotMissing] = dim;
                            ++offsetNotMissing;
                        }
                    }

                    anyMissing[index] = true;
                    missingIndices[index] = new HashedIntArray(missing, notMissing);

                    if (count == dimTrait) {
                        allMissing[index] = true;
                    }
                } else {
                    anyMissing[index] = false;
                    allMissing[index] = false;
                    missingIndices[index] = null;
                }
            }
        }
    }

    private int getIndex(final int tip, final int trait) {
        return tip * numTraits + trait;
    }

    private boolean isObservationMissing(final int index, final int dim) {
        final int id = index * dimTrait + dim;
        return rawMissingIndices.contains(id);
    }

    public class HashedIntArray {
        final private int[] array;
        final private int[] complement;

        HashedIntArray(final int[] array, final int[] complement) {
            this.array = array;
            this.complement = complement;
        }

        public int[] getArray() {
            return array;
        }

        public int[] getComplement() {
            return complement;
        }

        public int get(int index) {
            return array[index];
        }

        public int getLength() {
            return array.length;
        }

        public int getComplementLength() {
            return complement.length;
        }

        @Override
        public int hashCode() {
            return Arrays.hashCode(array);
        }

        @Override
        public boolean equals(Object obj) {
            if (obj instanceof int[]) {
                return Arrays.equals(array, (int[]) obj);
            } else {
                return false;
            }
        }

        public String toString() {
            return new Vector(array).toString();
        }
    }

    final private int tipCount;
    final private int numTraits;
    final private int dimTrait;

    final private List<Integer> rawMissingIndices;

    final private boolean[] anyMissing;
    final private boolean[] allMissing;
    final private HashedIntArray[] missingIndices;
}

