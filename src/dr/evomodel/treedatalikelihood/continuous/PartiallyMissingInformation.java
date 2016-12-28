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

import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.Matrix;

import java.util.Arrays;

/**
 * @author Marc A. Suchard
 */

public class PartiallyMissingInformation {

    public class HashedIntArray {
        final private int[] array;
        final private int[] complement;

//        HashedIntArray(final int[] array, final int dim) {
//            this.array = array;
//            this.complement = makeComplement(array, dim);
//        }

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

        //        private int[] makeComplement(final int[] array, final int dim) {
//            int[] complemenet = new int[dim - array.length];
//
//        }
    }

//    public class HashedPairIntArray {
//        final private HashedIntArray missing;
//        final private HashedIntArray notMissing;
//
//        HashedPairIntArray(final HashedIntArray missing, final HashedIntArray notMissing) {
//            this.missing = missing;
//            this.notMissing = notMissing;
//        }
//    }

    final private int tipCount;
    final private int numTraits;
    final private int dimTrait;
    final private Parameter missingParameter;

    final private boolean[] anyMissing;
    final private HashedIntArray[] missingIndices;
//    final private HashedIntArray[] notMissingIndices;

    public PartiallyMissingInformation(int tipCount,
                                       int numTraits,
                                       int dimTrait,
                                       Parameter missingParameter) {
        this.tipCount = tipCount;
        this.numTraits = numTraits;
        this.dimTrait = dimTrait;
        this.missingParameter = missingParameter;

        assert (tipCount * numTraits * dimTrait == missingParameter.getDimension());

        final int length = tipCount * numTraits;
        anyMissing = new boolean[length];
        missingIndices = new HashedIntArray[length];
//        notMissingIndices = new HashedIntArray[length];

        setupIndices();
    }

    public boolean isPartiallyMissing(final int tip, final int trait) {
        return anyMissing[getIndex(tip, trait)];
    }

    public HashedIntArray getMissingIndices(final int tip, final int trait) {
        return missingIndices[getIndex(tip, trait)];
    }

//    public HashedIntArray getNotMissingIndices(final int tip, final int trait) {
//        return notMissingIndices[getIndex(tip, trait)];
//    }

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

                System.err.println("count = " + count + " for " + index);

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
//                    notMissingIndices[index] = new HashedIntArray(notMissing);
                } else {
                    anyMissing[index] = false;
                    missingIndices[index] = null;
//                    notMissingIndices[index] = null;
                }
            }
        }
    }

    private int getIndex(final int tip, final int trait) {
        return tip * numTraits + trait;
    }

    private boolean isObservationMissing(final int index, final int dim) {
        final int id = index * dimTrait + dim;
//        System.err.println("id = " + id + " " +  missingParameter.getParameterValue(id));
        return missingParameter.getParameterValue(id) == 1;
    }

//    public static Matrix extractPartialVarianceMatrix(final double[][] variance, final int[] indices) {
//
//        final int varianceLength = indices.length;
//        double[][] var = new double[varianceLength][varianceLength];
//
//        for (int i = 0; i < varianceLength; ++i) {
//
//            final double[] in = variance[indices[i]];
//            final double[] out = var[i];
//
//            for (int j = 0; j < varianceLength; ++i) {
//                out[j] = in[indices[j]];
//            }
//        }
//        return new Matrix(var);
//    }

//    public static Matrix extractPartialVariance(final Matrix variance, final int[] indices) {
//        Matrix partialVariance = variance.extractRowsAndColumns(indices, indices);
//        return partialVariance.inverse();
//    }
}

