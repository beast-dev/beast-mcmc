/*
 * BeagleBastaLikelihoodDelegate.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.coalescent.basta;

import beagle.Beagle;
import beagle.BeagleFlag;
import beagle.basta.BeagleBasta;
import beagle.basta.BastaFactory;
import dr.evolution.tree.Tree;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.BufferIndexHelper;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.math.matrixAlgebra.WrappedVector;

import java.util.Arrays;
import java.util.List;

import static beagle.basta.BeagleBasta.BASTA_OPERATION_SIZE;

/**
 * @author Marc A. Suchard
 */
public class BeagleBastaLikelihoodDelegate extends BastaLikelihoodDelegate.AbstractBastaLikelihoodDelegate {

    private static final int COALESCENT_PROBABILITY_INDEX = 0;

    private final BeagleBasta beagle;

    private final BufferIndexHelper eigenBufferHelper;
    private final OffsetBufferIndexHelper populationSizesBufferHelper;

    public BeagleBastaLikelihoodDelegate(String name,
                                         Tree tree,
                                         int stateCount,
                                         boolean transpose) {
        super(name, tree, stateCount, transpose);

        int partialsCount = maxNumCoalescentIntervals * (tree.getNodeCount() + 1); // TODO much too large
        int matricesCount = maxNumCoalescentIntervals; // TODO much too small (except for strict-clock)

        int coalescentBufferCount = 5; // E, F, G, H, probabilities

        long requirementFlags = 0L;
        requirementFlags |= BeagleFlag.EIGEN_COMPLEX.getMask();

        beagle = BastaFactory.loadBastaInstance(0, coalescentBufferCount, maxNumCoalescentIntervals,
                partialsCount, 0, stateCount,
                1, 2, matricesCount, 1,
                1, null, 0L, requirementFlags);

        eigenBufferHelper = new BufferIndexHelper(1, 0);
        populationSizesBufferHelper = new OffsetBufferIndexHelper(1, 0, 0);

        beagle.setCategoryRates(new double[] { 1.0 });
    }

    @Override
    protected void computeBranchIntervalOperations(List<Integer> intervalStarts,
                                                   List<BranchIntervalOperation> branchIntervalOperations) {

        int[] operations = new int[branchIntervalOperations.size() * BASTA_OPERATION_SIZE]; // TODO instantiate once
        int[] intervals = new int[intervalStarts.size()]; // TODO instantiate once
        double[] lengths = new double[intervalStarts.size() - 1]; // TODO instantiate once

        vectorizeBranchIntervalOperations(intervalStarts, branchIntervalOperations, operations, intervals, lengths);

        int populationSizeIndex = populationSizesBufferHelper.getOffsetIndex(0);

        beagle.updateBastaPartials(operations, branchIntervalOperations.size(),
                intervals, intervalStarts.size(),
                populationSizeIndex, COALESCENT_PROBABILITY_INDEX);

        if (PRINT_COMMANDS) {
            double[] partials = new double[stateCount];
            for (BranchIntervalOperation operation : branchIntervalOperations) {
                getPartials(operation.outputBuffer, partials);
                System.err.println(operation + " " + new WrappedVector.Raw(partials));
            }

            double[] coalescent = new double[maxNumCoalescentIntervals];
            beagle.getBastaBuffer(COALESCENT_PROBABILITY_INDEX, coalescent);
            System.err.println("coalescent: " + new WrappedVector.Raw(coalescent));
        }
    }

    @Override
    String getStamp() { return "beagle"; }


    @Override
    protected void computeTransitionProbabilityOperations(List<TransitionMatrixOperation> matrixOperations) {

        int[] transitionMatrixIndices = new int[matrixOperations.size()]; // TODO instantiate once
        double[] branchLengths = new double[matrixOperations.size()]; // TODO instantiate once

        vectorizeTransitionMatrixOperations(matrixOperations, transitionMatrixIndices, branchLengths);

        int eigenIndex = eigenBufferHelper.getOffsetIndex(0);

        beagle.updateTransitionMatrices(eigenIndex, transitionMatrixIndices, null, null,
                branchLengths, matrixOperations.size());

        if (PRINT_COMMANDS) {
            double[] matrix = new double[stateCount * stateCount];
            for (TransitionMatrixOperation operation : matrixOperations) {
                getMatrix(operation.outputBuffer, matrix);
                System.err.println(operation + " " + new WrappedVector.Raw(matrix));
            }
        }
    }

    @Override
    protected double computeCoalescentIntervalReduction(List<Integer> intervalStarts,
                                                        List<BranchIntervalOperation> branchIntervalOperations) {

        // TODO vectorize once per likelihood-calculation
        int[] operations = new int[branchIntervalOperations.size() * BASTA_OPERATION_SIZE]; // TODO instantiate once
        int[] intervals = new int[intervalStarts.size()]; // TODO instantiate once
        double[] lengths = new double[intervalStarts.size() - 1]; // TODO instantiate once

        vectorizeBranchIntervalOperations(intervalStarts, branchIntervalOperations, operations, intervals, lengths);

        int populationSizeIndex = populationSizesBufferHelper.getOffsetIndex(0);

        double[] result = new double[1];

        beagle.accumulateBastaPartials(operations, branchIntervalOperations.size(),
                intervals, intervalStarts.size(), lengths,
                populationSizeIndex, COALESCENT_PROBABILITY_INDEX,
                result);

        return result[0];
    }

    @Override
    protected void computeBranchIntervalOperationsGrad(List<Integer> intervalStarts, List<TransitionMatrixOperation> matrixOperations, List<BranchIntervalOperation> branchIntervalOperations) {
        int[] operations = new int[branchIntervalOperations.size() * BASTA_OPERATION_SIZE]; // TODO instantiate once
        int[] intervals = new int[intervalStarts.size()]; // TODO instantiate once
        double[] lengths = new double[intervalStarts.size() - 1]; // TODO instantiate once

        vectorizeBranchIntervalOperations(intervalStarts, branchIntervalOperations, operations, intervals, lengths);

        int populationSizeIndex = populationSizesBufferHelper.getOffsetIndex(0);

        beagle.updateBastaPartialsGrad(operations, branchIntervalOperations.size(),
                intervals, intervalStarts.size(),
                populationSizeIndex, COALESCENT_PROBABILITY_INDEX);
    }

    @Override
    protected void computeTransitionProbabilityOperationsGrad(List<TransitionMatrixOperation> matrixOperations) {
        int[] transitionMatrixIndices = new int[matrixOperations.size()]; // TODO instantiate once
        double[] branchLengths = new double[matrixOperations.size()]; // TODO instantiate once

        vectorizeTransitionMatrixOperations(matrixOperations, transitionMatrixIndices, branchLengths);

        beagle.updateTransitionMatricesGrad(transitionMatrixIndices,
                branchLengths, matrixOperations.size());
    }

    @Override
    protected double[][] computeCoalescentIntervalReductionGrad(List<Integer> intervalStarts, List<BranchIntervalOperation> branchIntervalOperations) {
        int[] operations = new int[branchIntervalOperations.size() * BASTA_OPERATION_SIZE]; // TODO instantiate once
        int[] intervals = new int[intervalStarts.size()]; // TODO instantiate once
        double[] lengths = new double[intervalStarts.size() - 1]; // TODO instantiate once

        vectorizeBranchIntervalOperations(intervalStarts, branchIntervalOperations, operations, intervals, lengths);

        int populationSizeIndex = populationSizesBufferHelper.getOffsetIndex(0);

        double[] out = new double[stateCount*stateCount];
        double[][] result = new double[stateCount][stateCount];
        beagle.accumulateBastaPartialsGrad(operations, branchIntervalOperations.size(),
                intervals, intervalStarts.size(), lengths,
                populationSizeIndex, COALESCENT_PROBABILITY_INDEX,
                out);

        for (int i = 0; i < stateCount; ++i) {
            for (int j = 0; j < stateCount; ++j) {
                result[i][j] = out[i*stateCount+j];
            }
        }

        return result;
    }

    @Override
    protected double[] computeCoalescentIntervalReductionGradPopSize(List<Integer> intervalStarts, List<BranchIntervalOperation> branchIntervalOperations) {
        return new double[0];
    }

    @Override
    public void setPartials(int index, double[] partials) {
        beagle.setPartials(index, partials);
    }

    @Override
    public void getPartials(int index, double[] partials) {
        assert index >= 0;
        assert partials != null;
        assert partials.length >= stateCount;

        beagle.getPartials(index, Beagle.NONE, partials);
    }

    public void getMatrix(int index, double[] matrix) {
        assert index >= 0;
        assert matrix != null;
        assert matrix.length >= stateCount * stateCount;

        beagle.getTransitionMatrix(index, matrix);
    }

    @Override
    public void updateEigenDecomposition(int index, EigenDecomposition decomposition, boolean flip) {
        if (flip) {
            eigenBufferHelper.flipOffset(0);
        }

        if (transpose) {
            decomposition = decomposition.transpose();
        }

        beagle.setEigenDecomposition(
                eigenBufferHelper.getOffsetIndex(0),
                decomposition.getEigenVectors(),
                decomposition.getInverseEigenVectors(),
                decomposition.getEigenValues());
    }

    @Override
    public void updatePopulationSizes(int index, double[] sizes, boolean flip) {
        if (flip) {
            populationSizesBufferHelper.flipOffset(0);
        }

        beagle.setStateFrequencies(populationSizesBufferHelper.getOffsetIndex(0), sizes);
    }

    private void vectorizeTransitionMatrixOperations(List<TransitionMatrixOperation> matrixOperations,
                                                     int[] transitionMatrixIndices,
                                                     double[] branchLengths) {
        int k = 0;
        for (TransitionMatrixOperation op : matrixOperations) {
            transitionMatrixIndices[k] = op.outputBuffer; // TODO double-buffer
            branchLengths[k] = op.time;
            ++k;
        }
    }

    int tipCount = -1;
    int[] map;
    int used;

    private static final boolean CACHE_FRIENDLY = true;

    int map(int buffer) {
        if (buffer < tipCount) {
            return buffer;
        } else {
            if (map[buffer] == -1) {
                map[buffer] = used;
                ++used;
            }
            return map[buffer];
        }
    }


    private void vectorizeBranchIntervalOperations(List<Integer> intervalStarts,
                                                   List<BranchIntervalOperation> branchIntervalOperations,
                                                   int[] operations,
                                                   int[] intervals,
                                                   double[] lengths) {

        if (CACHE_FRIENDLY) {
            this.tipCount = tree.getExternalNodeCount();
            if (this.map == null) {
                this.map = new int[maxNumCoalescentIntervals * (tree.getNodeCount() + 1)];
            }
            Arrays.fill(this.map, -1);
            this.used = this.tipCount;
        }

        // TODO double-buffer
        int k = 0;
        for (BranchIntervalOperation op : branchIntervalOperations) {

            if (CACHE_FRIENDLY) {
                operations[k] = map(op.outputBuffer);
                operations[k + 1] = map(op.inputBuffer1);
                operations[k + 2] = op.inputMatrix1;
                operations[k + 3] = map(op.inputBuffer2);
                operations[k + 4] = op.inputMatrix2;
                operations[k + 5] = map(op.accBuffer1);
                operations[k + 6] = map(op.accBuffer2);
                operations[k + 7] = op.intervalNumber;
            } else {
                operations[k] = op.outputBuffer;
                operations[k + 1] = op.inputBuffer1;
                operations[k + 2] = op.inputMatrix1;
                operations[k + 3] = op.inputBuffer2;
                operations[k + 4] = op.inputMatrix2;
                operations[k + 5] = op.accBuffer1;
                operations[k + 6] = op.accBuffer2;
                operations[k + 7] = op.intervalNumber;
            }

            k += BASTA_OPERATION_SIZE;
        }

//        int i = 0;
//        for (int start : intervalStarts) {
//            intervals[i] = start;
//            lengths[i] = branchIntervalOperations.get(start).intervalLength;
//            ++i;
//        }

        int i = 0;
        for (int end = intervalStarts.size() - 1; i < end; ++i) {
            int start = intervalStarts.get(i);
            intervals[i] = start;
            lengths[i] = branchIntervalOperations.get(start).intervalLength;
        }
        intervals[i] = intervalStarts.get(i);
    }

    private boolean releaseSingleton = true;

    private void releaseBeagle() throws Throwable {
        if (beagle != null && releaseSingleton) {
            beagle.finalize();
            releaseSingleton = false;
        }
    }

    public static void releaseBeagleBastaLikelihoodDelegate(BastaLikelihood treeDataLikelihood) throws Throwable {
        BastaLikelihoodDelegate likelihoodDelegate = treeDataLikelihood.getLikelihoodDelegate();
        if (likelihoodDelegate instanceof BeagleBastaLikelihoodDelegate) {
            BeagleBastaLikelihoodDelegate delegate = (BeagleBastaLikelihoodDelegate) likelihoodDelegate;
            delegate.releaseBeagle();
        }
    }

    public static void releaseAllBeagleBastaInstances() throws Throwable {
        for (Likelihood likelihood : dr.inference.model.Likelihood.FULL_LIKELIHOOD_SET) {
            if (likelihood instanceof BastaLikelihood) {
                releaseBeagleBastaLikelihoodDelegate((BastaLikelihood) likelihood);
            } else if (likelihood instanceof CompoundLikelihood) {
                for (Likelihood likelihood2: ((CompoundLikelihood) likelihood).getLikelihoods()) {
                    if (likelihood2 instanceof BastaLikelihood) {
                        releaseBeagleBastaLikelihoodDelegate((BastaLikelihood) likelihood2);
                    }
                }
            }
        }
    }

    static class OffsetBufferIndexHelper extends BufferIndexHelper {

        public OffsetBufferIndexHelper(int maxIndexValue, int minIndexValue, int bufferSetNumber) {
            super(maxIndexValue, minIndexValue, bufferSetNumber);
        }

        @Override
        protected int computeOffset(int offset) { return offset; }
    }
}
