/*
 * GenericBastaLikelihoodDelegate.java
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

import dr.evolution.tree.Tree;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.math.matrixAlgebra.WrappedVector;

import java.util.Arrays;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class GenericBastaLikelihoodDelegate extends BastaLikelihoodDelegate.AbstractBastaLikelihoodDelegate {

    private final double[] partials;
    private final double[] matrices;
    private final double[] sizes;
    private final double[] coalescent;

    private final double[] e;
    private final double[] f;
    private final double[] g;
    private final double[] h;

    private final double[] temp;

    private final EigenDecomposition[] decompositions; // TODO flatten?
    private final double[][][] partialsGrad;
    private final double[][][] matricesGrad;
    private final double[][][] coalescentGrad;
    private final double[][][] eGrad;
    private final double[][][] fGrad;
    private final double[][][] gGrad;
    private final double[][][] hGrad;
    private final double[][] partialsGradPopSize;
    private final double[][] coalescentGradPopSize;
    private final double[][] eGradPopSize;
    private final double[][] fGradPopSize;
    private final double[][] gGradPopSize;
    private final double[][] hGradPopSize;


    public GenericBastaLikelihoodDelegate(String name,
                                          Tree tree,
                                          int stateCount,
                                          boolean transpose) {
        super(name, tree, stateCount, transpose);

        this.partials = new double[maxNumCoalescentIntervals * (tree.getNodeCount() + 1) * stateCount]; // TODO much too large
        this.matrices = new double[maxNumCoalescentIntervals * stateCount * stateCount]; // TODO much too small (except for strict-clock)
        this.coalescent = new double[maxNumCoalescentIntervals];
        this.sizes = new double[2 * stateCount];
        this.decompositions = new EigenDecomposition[1];

        this.e = new double[maxNumCoalescentIntervals * stateCount];
        this.f = new double[maxNumCoalescentIntervals * stateCount];
        this.g = new double[maxNumCoalescentIntervals * stateCount];
        this.h = new double[maxNumCoalescentIntervals * stateCount];

        this.partialsGrad = new double[stateCount][stateCount][maxNumCoalescentIntervals * (tree.getNodeCount() + 1) * stateCount];
        this.matricesGrad = new double[stateCount][stateCount][maxNumCoalescentIntervals * stateCount * stateCount];
        this.coalescentGrad = new double[stateCount][stateCount][maxNumCoalescentIntervals];

        this.eGrad = new double[stateCount][stateCount][maxNumCoalescentIntervals * stateCount];
        this.fGrad = new double[stateCount][stateCount][maxNumCoalescentIntervals * stateCount];
        this.gGrad = new double[stateCount][stateCount][maxNumCoalescentIntervals * stateCount];
        this.hGrad = new double[stateCount][stateCount][maxNumCoalescentIntervals * stateCount];

        this.partialsGradPopSize = new double[stateCount][maxNumCoalescentIntervals * (tree.getNodeCount() + 1) * stateCount];
        this.coalescentGradPopSize = new double[stateCount][maxNumCoalescentIntervals];

        this.eGradPopSize = new double[stateCount][maxNumCoalescentIntervals * stateCount];
        this.fGradPopSize = new double[stateCount][maxNumCoalescentIntervals * stateCount];
        this.gGradPopSize = new double[stateCount][maxNumCoalescentIntervals * stateCount];
        this.hGradPopSize = new double[stateCount][maxNumCoalescentIntervals * stateCount];
        this.temp = new double[stateCount * stateCount];
    }

    @Override
    protected void computeBranchIntervalOperations(List<Integer> intervalStarts,
                                                   List<BranchIntervalOperation> branchIntervalOperations) {

        Arrays.fill(coalescent, 0.0);

        for (int interval = 0; interval < intervalStarts.size() - 1; ++interval) { // execute in series by intervalNumber
            // TODO try grouping by executionOrder (unclear if more efficient, same total #)
            int start = intervalStarts.get(interval);
            int end = intervalStarts.get(interval + 1);

            computeInnerBranchIntervalOperations(branchIntervalOperations, start, end);
        }

        if (PRINT_COMMANDS) {
            System.err.println("coalescent: " + new WrappedVector.Raw(coalescent));
        }
    }

    protected void computeInnerBranchIntervalOperations(List<BranchIntervalOperation> branchIntervalOperations,
                                                      int start, int end) {

        for (int i = start; i < end; ++i) {
            BranchIntervalOperation operation = branchIntervalOperations.get(i);

            peelPartials(
                    partials, operation.outputBuffer,
                    operation.inputBuffer1, operation.inputBuffer2,
                    matrices,
                    operation.inputMatrix1, operation.inputMatrix2,
                    operation.accBuffer1, operation.accBuffer2,
                    coalescent, operation.intervalNumber,
                    sizes, 0,
                    stateCount);

            if (PRINT_COMMANDS) {
                System.err.println(operation + " " +
                        new WrappedVector.Raw(partials, operation.outputBuffer * stateCount, stateCount));
            }
        }
    }

    @Override
    protected void computeTransitionProbabilityOperations(List<TransitionMatrixOperation> matrixOperations) {
        computeInnerTransitionProbabilityOperations(matrixOperations, 0, matrixOperations.size(), temp);
    }

    protected void computeInnerTransitionProbabilityOperations(List<TransitionMatrixOperation> matrixOperations,
                                                               int start, int end, double[] temp) {
        for (int i = start; i < end; ++i) {
            TransitionMatrixOperation operation = matrixOperations.get(i);
            computeTransitionProbabilities(
                    operation.time,
                    matrices, operation.outputBuffer * stateCount * stateCount,
                    decompositions[operation.decompositionBuffer],
                    stateCount, temp);

            if (PRINT_COMMANDS) {
                System.err.println(operation + " " +
                        new WrappedVector.Raw(matrices,
                                operation.outputBuffer * stateCount * stateCount, stateCount * stateCount));
            }
        }
    }

    @Override
    String getStamp() { return "generic"; }

    @Override
    protected double computeCoalescentIntervalReduction(List<Integer> intervalStarts,
                                                        List<BranchIntervalOperation> branchIntervalOperations) {

        Arrays.fill(e, 0.0);
        Arrays.fill(f, 0.0);
        Arrays.fill(g, 0.0);
        Arrays.fill(h, 0.0);

        for (int interval = 0; interval < intervalStarts.size() - 1; ++interval) { // TODO execute in parallel (no race conditions)
            int start = intervalStarts.get(interval);
            int end = intervalStarts.get(interval + 1);

            for (int i = start; i < end; ++i) { // TODO execute in parallel (has race conditions)
                BranchIntervalOperation operation = branchIntervalOperations.get(i);
                reduceWithinInterval(e, f, g, h, partials,
                        operation.inputBuffer1, operation.inputBuffer2,
                        operation.accBuffer1, operation.accBuffer2,
                        operation.intervalNumber,
                        stateCount);

            }
        }

        double logL = 0.0;
        for (int i = 0; i < intervalStarts.size() - 1; ++i) { // TODO execute in parallel
            BranchIntervalOperation operation = branchIntervalOperations.get(intervalStarts.get(i));

            logL += reduceAcrossIntervals(e, f, g, h,
                    operation.intervalNumber, operation.intervalLength,
                    sizes, coalescent, stateCount);
        }

        return logL;
    }

    @Override
    protected void computeBranchIntervalOperationsGrad(List<Integer> intervalStarts, List<TransitionMatrixOperation> matrixOperations, List<BranchIntervalOperation> branchIntervalOperations) {
        for (int interval = 0; interval < intervalStarts.size() - 1; ++interval) { // execute in series by intervalNumber
            // TODO try grouping by executionOrder (unclear if more efficient, same total #)
            int start = intervalStarts.get(interval);
            int end = intervalStarts.get(interval + 1);
            computeInnerBranchIntervalOperationsGrad(branchIntervalOperations, matrixOperations, start, end);
        }
    }

    protected void computeInnerBranchIntervalOperationsGrad(List<BranchIntervalOperation> branchIntervalOperations, List<TransitionMatrixOperation> matrixOperations, int start, int end) {
        for (int i = start; i < end; ++i) {
            BranchIntervalOperation operation = branchIntervalOperations.get(i);
            TransitionMatrixOperation Matrixoperation = matrixOperations.get(operation.intervalNumber);

            peelPartialsGrad(
                    partials, Matrixoperation.time, operation.outputBuffer,
                    operation.inputBuffer1, operation.inputBuffer2,
                    matrices,
                    operation.inputMatrix1, operation.inputMatrix2,
                    operation.accBuffer1, operation.accBuffer2,
                    coalescent, operation.intervalNumber,
                    sizes, 0,
                    stateCount);

        }
    }

    private void peelPartialsGrad(double[] partials,
                                  double distance, int resultOffset,
                                  int leftPartialOffset, int rightPartialOffset,
                                  double[] matrices,
                                  int leftMatrixOffset, int rightMatrixOffset,
                                  int leftAccOffset, int rightAccOffset,
                                  double[] probability, int probabilityOffset,
                                  double[] sizes, int sizesOffset,
                                  int stateCount) {

        resultOffset *= stateCount;

        // Handle left
        leftPartialOffset *= stateCount;
        leftMatrixOffset *= stateCount * stateCount;
        leftAccOffset *= stateCount;

        for (int a = 0; a < stateCount; ++a) {
            for (int b = 0; b < stateCount; ++b) {
                for (int i = 0; i < stateCount; ++i) {
                    double sum = 0.0;
                    if (transpose && i == b) {
                        sum += partials[leftAccOffset + a] * distance;
                    }
                    if (!transpose) {
                        for (int j = 0; j < stateCount; ++j) {
                            sum += matricesGrad[a][b][leftMatrixOffset + i * stateCount + j] * partials[leftPartialOffset + j];
                        }
                    }
                    for (int j = 0; j < stateCount; ++j) {
                        sum += matrices[leftMatrixOffset + i * stateCount + j] * partialsGrad[a][b][leftPartialOffset + j];
                    }
                    partialsGrad[a][b][resultOffset + i] = sum;

                    throw new RuntimeException("Function should not depend on `transpose`");
                }
            }
        }

        for (int a = 0; a < stateCount; ++a) {
            for (int i = 0; i < stateCount; ++i) {
                double sum = 0.0;
                for (int j = 0; j < stateCount; ++j) {
                    sum += matrices[leftMatrixOffset + i * stateCount + j] * partialsGradPopSize[a][leftPartialOffset + j];
                }
                partialsGradPopSize[a][resultOffset + i] = sum;
            }
        }

        if (rightPartialOffset >= 0) {
            // Handle right
            rightPartialOffset *= stateCount;
            rightMatrixOffset *= stateCount * stateCount;

            rightAccOffset *= stateCount;
            // TODO: check bug?
            sizesOffset *= sizesOffset * stateCount;

            for (int a = 0; a < stateCount; ++a) {
                for (int b = 0; b < stateCount; ++b) {
                    double J = probability[probabilityOffset];

                    // first half
                    double partial_J_ab = 0.0;
                    for (int i = 0; i < stateCount; ++i) {
                        double rightGrad = 0.0;
                        if (transpose && i == b) {
                            rightGrad += partials[rightAccOffset + a] * distance;
                        }

                        if (!transpose) {
                            for (int j = 0; j < stateCount; ++j) {
                                rightGrad += matricesGrad[a][b][rightMatrixOffset + i * stateCount + j] * partials[rightPartialOffset + j];
                            }
                        }
                        for (int j = 0; j < stateCount; ++j) {
                            rightGrad += matrices[rightMatrixOffset + i * stateCount + j] * partialsGrad[a][b][rightPartialOffset + j];
                        }
                        double leftGrad = partialsGrad[a][b][resultOffset + i];
                        double left = partials[leftAccOffset + i];
                        double right = partials[rightAccOffset + i];

                        double entry = (leftGrad * right + rightGrad * left) / sizes[sizesOffset + i];
                        partial_J_ab += entry;

                        partialsGrad[a][b][resultOffset + i] = entry / J;
                        partialsGrad[a][b][leftAccOffset + i] = leftGrad;
                        partialsGrad[a][b][rightAccOffset + i] = rightGrad;

                        throw new RuntimeException("Function should not depend on `transpose`");
                    }
                    // second half
                    for (int i = 0; i < stateCount; ++i) {
                        double entry = partials[resultOffset + i];
                        partialsGrad[a][b][resultOffset + i] -= partial_J_ab * entry / J;
                    }
                    coalescentGrad[a][b][probabilityOffset] = partial_J_ab;
                }
            }

            for (int a = 0; a < stateCount; ++a) {
                double J = probability[probabilityOffset];
                // first half
                double partial_J_ab_PopSize = 0.0;
                for (int i = 0; i < stateCount; ++i) {
                    double rightGradPopSize = 0.0;
                    for (int j = 0; j < stateCount; ++j) {
                        rightGradPopSize += matrices[rightMatrixOffset + i * stateCount + j] *  partialsGradPopSize[a][rightPartialOffset + j];
                    }
                    double leftGradPopSize = partialsGradPopSize[a][resultOffset + i];
                    double left = partials[leftAccOffset + i];
                    double right = partials[rightAccOffset + i];

                    double entry = (leftGradPopSize * right + rightGradPopSize * left) / sizes[sizesOffset + i];
                    if (i == a){
                        entry += left * right;
                    }
                    partial_J_ab_PopSize += entry;

                    partialsGradPopSize[a][resultOffset + i] = entry / J;
                    partialsGradPopSize[a][leftAccOffset + i] = leftGradPopSize;
                    partialsGradPopSize[a][rightAccOffset + i] = rightGradPopSize;
                }
                // second half
                for (int i = 0; i < stateCount; ++i) {
                    double entry = partials[resultOffset + i];
                    partialsGradPopSize[a][resultOffset + i] -= partial_J_ab_PopSize * entry / J;
                }
                coalescentGradPopSize[a][probabilityOffset] = partial_J_ab_PopSize;
            }
        }
    }


    @Override
    protected void computeTransitionProbabilityOperationsGrad(List<TransitionMatrixOperation> matrixOperations) {
        computeInnerTransitionProbabilityOperationsGrad(matrixOperations, 0, matrixOperations.size());
    }

    protected void computeInnerTransitionProbabilityOperationsGrad(List<TransitionMatrixOperation> matrixOperations,
                                                               int start, int end) {
        for (int i = start; i < end; ++i) {
            TransitionMatrixOperation operation = matrixOperations.get(i);
            computeTransitionProbabilitiesGrad(operation.time, operation.outputBuffer * stateCount * stateCount);
        }
    }

    private void computeTransitionProbabilitiesGrad(double distance, int matrixOffset) {
        for (int a = 0; a < stateCount; a++) {
            for (int b = 0; b < stateCount; b++) {
                for (int c = 0; c < stateCount; c++) {
                    for (int d = 0; d < stateCount; d++) { // TODO MAS: last loop unnecessary (also S^4 storage is unnecessary)
                        if (d == b) {
                            // TODO MAS: should these be cached at all? why not generate on the fly (t * matrices[])
                            matricesGrad[a][b][matrixOffset + c*stateCount + b] =  distance * matrices[matrixOffset + c*stateCount + a];
                        } else {
                            matricesGrad[a][b][matrixOffset + c*stateCount + d] = 0; // TODO MAS: avoid caching (many) zeros
                        }
                    }
                }
            }
        }
    }

    @Override
    protected double[][] computeCoalescentIntervalReductionGrad(List<Integer> intervalStarts, List<BranchIntervalOperation> branchIntervalOperations) {

        for (int interval = 0; interval < intervalStarts.size() - 1; ++interval) { // TODO execute in parallel (no race conditions)
            int start = intervalStarts.get(interval);
            int end = intervalStarts.get(interval + 1);

            for (int i = start; i < end; ++i) { // TODO execute in parallel (has race conditions)
                BranchIntervalOperation operation = branchIntervalOperations.get(i);
                reduceWithinIntervalGrad(partials, partialsGrad,
                        operation.inputBuffer1, operation.inputBuffer2,
                        operation.accBuffer1, operation.accBuffer2,
                        operation.intervalNumber,
                        stateCount);

            }
        }

        double[][] grad = new double[stateCount][stateCount];
        for (int i = 0; i < intervalStarts.size() - 1; ++i) { // TODO execute in parallel
            BranchIntervalOperation operation = branchIntervalOperations.get(intervalStarts.get(i));

            double[][] temp_grad =  reduceAcrossIntervalsGrad(e, f, g, h,
                    operation.intervalNumber, operation.intervalLength,
                    sizes, coalescent, stateCount);

            for (int a = 0; a < stateCount; a++) {
                for (int b = 0; b < stateCount; b++) {
                    grad[a][b] += temp_grad[a][b];
                }
            }
        }
        return grad;
    }

    @Override
    protected double[] computeCoalescentIntervalReductionGradPopSize(List<Integer> intervalStarts, List<BranchIntervalOperation> branchIntervalOperations) {

        Arrays.stream(eGradPopSize).forEach(a -> Arrays.fill(a, 0));
        Arrays.stream(fGradPopSize).forEach(a -> Arrays.fill(a, 0));
        Arrays.stream(gGradPopSize).forEach(a -> Arrays.fill(a, 0));
        Arrays.stream(hGradPopSize).forEach(a -> Arrays.fill(a, 0));

        for (int interval = 0; interval < intervalStarts.size() - 1; ++interval) { // TODO execute in parallel (no race conditions)
            int start = intervalStarts.get(interval);
            int end = intervalStarts.get(interval + 1);

            for (int i = start; i < end; ++i) { // TODO execute in parallel (has race conditions)
                BranchIntervalOperation operation = branchIntervalOperations.get(i);
                reduceWithinIntervalGrad(partials, partialsGrad,
                        operation.inputBuffer1, operation.inputBuffer2,
                        operation.accBuffer1, operation.accBuffer2,
                        operation.intervalNumber,
                        stateCount);
            }
        }

        double[] grad = new double[stateCount];
        for (int i = 0; i < intervalStarts.size() - 1; ++i) { // TODO execute in parallel
            BranchIntervalOperation operation = branchIntervalOperations.get(intervalStarts.get(i));

            double[] temp_grad =  reduceAcrossIntervalsGradPopSize(e, f, g, h,
                    operation.intervalNumber, operation.intervalLength,
                    sizes, coalescent, stateCount);

            for (int a = 0; a < stateCount; a++) {
                    grad[a] += temp_grad[a];
                }
            }
        return grad;
    }

    private double[][] reduceAcrossIntervalsGrad(double[] e, double[] f, double[] g, double[] h,
                                         int interval, double length,
                                         double[] sizes, double[] coalescent,
                                         int stateCount) {

        int offset = interval * stateCount;
        double[][] grad = new double[stateCount][stateCount];

        for (int a = 0; a < stateCount; a++) {
            for (int b = 0; b < stateCount; b++) {
                double sum = 0.0;
                for (int k = 0; k < stateCount; ++k) {
                    sum += (2 * e[offset + k] * eGrad[a][b][offset + k] - fGrad[a][b][offset + k] +
                            2 * g[offset + k] * gGrad[a][b][offset + k] - hGrad[a][b][offset + k]) / sizes[k];
                }

                grad[a][b] = -length * sum / 4;

                double J = coalescent[interval];
                if (J != 0.0) {
                    grad[a][b] += coalescentGrad[a][b][interval] / J;
                }
            }
        }
        return grad;
    }

    private double[] reduceAcrossIntervalsGradPopSize(double[] e, double[] f, double[] g, double[] h,
                                                 int interval, double length,
                                                 double[] sizes, double[] coalescent,
                                                 int stateCount) {

        int offset = interval * stateCount;
        double[] grad = new double[stateCount];

        for (int a = 0; a < stateCount; a++) {
                double sum = 0.0;
                for (int k = 0; k < stateCount; ++k) {
                    sum += (2 * e[offset + k] * eGradPopSize[a][offset + k] - fGradPopSize[a][offset + k] +
                            2 * g[offset + k] * gGradPopSize[a][offset + k] - hGradPopSize[a][offset + k]) / sizes[k];
                    if (k == a) {
                        sum += (e[offset + k] * e[offset + k]) - f[offset + k] +
                                (g[offset + k] * g[offset + k]) - h[offset + k];
                    }
                }


                grad[a] = -length * sum / 4;

                double J = coalescent[interval];
                if (J != 0.0) {
                    grad[a] += coalescentGradPopSize[a][interval] / J;
                }
            }
        return grad;
    }

    private void reduceWithinIntervalGrad(double[] partials, double[][][] partialsGrad,
                                          int startBuffer1, int startBuffer2,
                                          int endBuffer1, int endBuffer2,
                                          int interval, int stateCount)  {
        interval *= stateCount;

        startBuffer1 *= stateCount;
        endBuffer1 *= stateCount;

        for (int a = 0; a < stateCount; ++a) {
            for (int b = 0; b < stateCount; ++b) {
                for (int i = 0; i < stateCount; ++i) {
                    double startPGrad = partialsGrad[a][b][startBuffer1 + i];
                    double startP = partials[startBuffer1 + i];
                    eGrad[a][b][interval + i] += startPGrad;
                    fGrad[a][b][interval + i] += 2 * startP * startPGrad;

                    double endPGrad = partialsGrad[a][b][endBuffer1 + i];
                    double endP = partials[endBuffer1 + i];
                    gGrad[a][b][interval + i] += endPGrad;
                    hGrad[a][b][interval + i] += 2 * endP * endPGrad;
                }
            }
        }

        for (int a = 0; a < stateCount; ++a) {
                for (int i = 0; i < stateCount; ++i) {
                    double startPGradPopSize = partialsGradPopSize[a][startBuffer1 + i];
                    double startP = partials[startBuffer1 + i];
                    eGradPopSize[a][interval + i] += startPGradPopSize;
                    fGradPopSize[a][interval + i] += 2 * startP * startPGradPopSize;

                    double endPGradPopSize = partialsGradPopSize[a][endBuffer1 + i];
                    double endP = partials[endBuffer1 + i];
                    gGradPopSize[a][interval + i] += endPGradPopSize;
                    hGradPopSize[a][interval + i] += 2 * endP * endPGradPopSize;
                }
            }

        if (startBuffer2 >= 0) {
            startBuffer2 *= stateCount;
            endBuffer2 *= stateCount;
            for (int a = 0; a < stateCount; ++a) {
                for (int b = 0; b < stateCount; ++b) {
                    for (int i = 0; i < stateCount; ++i) {
                        double startPGrad = partialsGrad[a][b][startBuffer2 + i];
                        double startP = partials[startBuffer2 + i];
                        eGrad[a][b][interval + i] += startPGrad;
                        fGrad[a][b][interval + i] += 2 * startP * startPGrad;

                        double endPGrad = partialsGrad[a][b][endBuffer2 + i];
                        double endP = partials[endBuffer2 + i];
                        gGrad[a][b][interval + i] += endPGrad;
                        hGrad[a][b][interval + i] += 2 * endP * endPGrad;
                    }
                }
            }

            for (int a = 0; a < stateCount; ++a) {
                    for (int i = 0; i < stateCount; ++i) {
                        double startPGradPopSize = partialsGradPopSize[a][startBuffer2 + i];
                        double startP = partials[startBuffer2 + i];
                        eGradPopSize[a][interval + i] += startPGradPopSize;
                        fGradPopSize[a][interval + i] += 2 * startP * startPGradPopSize;

                        double endPGradPopSize = partialsGradPopSize[a][endBuffer2 + i];
                        double endP = partials[endBuffer2 + i];
                        gGradPopSize[a][interval + i] += endPGradPopSize;
                        hGradPopSize[a][interval + i] += 2 * endP * endPGradPopSize;
                    }

            }
        }
    }

    @Override
    public void setPartials(int index, double[] partials) {
        assert partials.length == stateCount;

        System.arraycopy(partials, 0, this.partials, index * stateCount, stateCount);
    }

    @Override
    public void updateEigenDecomposition(int index, EigenDecomposition decomposition, boolean flip) {

        if (transpose) {
            decomposition = decomposition.transpose();
        }

        decompositions[index] = decomposition;
    }

    @Override
    public void updatePopulationSizes(int index, double[] sizes, boolean flip) {
        assert sizes.length == stateCount;

        System.arraycopy(sizes, 0, this.sizes, index * stateCount, stateCount);
    }

    private static void peelPartials(double[] partials,
                                     int resultOffset,
                                     int leftPartialOffset, int rightPartialOffset,
                                     double[] matrices,
                                     int leftMatrixOffset, int rightMatrixOffset,
                                     int leftAccOffset, int rightAccOffset,
                                     double[] probability, int probabilityOffset,
                                     double[] sizes, int sizesOffset,
                                     int stateCount) {

        resultOffset *= stateCount;

        // Handle left
        leftPartialOffset *= stateCount;
        leftMatrixOffset *= stateCount * stateCount;

        for (int i = 0; i < stateCount; ++i) {
            double sum = 0.0;
            for (int j = 0; j < stateCount; ++j) {
                sum += matrices[leftMatrixOffset + i * stateCount + j] * partials[leftPartialOffset + j];
            }
            partials[resultOffset + i] = sum;
        }

        if (rightPartialOffset >= 0) {
            // Handle right
            rightPartialOffset *= stateCount;
            rightMatrixOffset *= stateCount * stateCount;

            leftAccOffset *= stateCount;
            rightAccOffset *= stateCount;

            sizesOffset *= sizesOffset * stateCount;

            double prob = 0.0;
            for (int i = 0; i < stateCount; ++i) {
                double right = 0.0;
                for (int j = 0; j < stateCount; ++j) {
                    right += matrices[rightMatrixOffset + i * stateCount + j] * partials[rightPartialOffset + j];
                }
                // entry = left * right * size
                double left = partials[resultOffset + i];
                double entry = left * right / sizes[sizesOffset + i];

                partials[resultOffset + i] = entry;
                partials[leftAccOffset + i] = left;
                partials[rightAccOffset + i] = right;

                prob += entry;
            }

            for (int i = 0; i < stateCount; ++i) {
                partials[resultOffset + i] /= prob;
            }

            probability[probabilityOffset] = prob;
        }

        // TODO rescale?
    }

    private static void computeTransitionProbabilities(double distance,
                                                       double[] matrix,
                                                       int matrixOffset,
                                                       EigenDecomposition eigen,
                                                       int stateCount,
                                                       double[] iexp) {

        assert matrix.length >= matrixOffset + stateCount * stateCount;
        assert iexp.length == stateCount * stateCount;

        boolean real = eigen.getEigenValues().length == stateCount;

        double[] Evec = eigen.getEigenVectors();
        double[] Eval = eigen.getEigenValues();
        double[] Ievc = eigen.getInverseEigenVectors();
        for (int i = 0; i < stateCount; i++) {

            if (real || Eval[stateCount + i] == 0) {
                // 1x1 block
                double temp = Math.exp(distance * Eval[i]);
                for (int j = 0; j < stateCount; j++) {
                    iexp[i * stateCount + j] = Ievc[i * stateCount + j] * temp;
                }
            } else {
                // 2x2 conjugate block
                // If A is 2x2 with complex conjugate pair eigenvalues a +/- bi, then
                // exp(At) = exp(at)*( cos(bt)I + \frac{sin(bt)}{b}(A - aI)).
                int i2 = i + 1;
                double b = Eval[stateCount + i];
                double expat = Math.exp(distance * Eval[i]);
                double expatcosbt = expat * Math.cos(distance * b);
                double expatsinbt = expat * Math.sin(distance * b);

                for (int j = 0; j < stateCount; j++) {
                    iexp[i * stateCount + j] = expatcosbt * Ievc[i * stateCount + j] +
                            expatsinbt * Ievc[i2 * stateCount + j];
                    iexp[i2 * stateCount + j] = expatcosbt * Ievc[i2 * stateCount + j] -
                            expatsinbt * Ievc[i * stateCount + j];
                }
                i++; // processed two conjugate rows
            }
        }

        int u = matrixOffset;
            for (int i = 0; i < stateCount; i++) {
                for (int j = 0; j < stateCount; j++) {
                    double temp = 0.0;
                    for (int k = 0; k < stateCount; k++) {
                        temp += Evec[i * stateCount + k] * iexp[k * stateCount + j];
                    }
                    matrix[u] = Math.abs(temp);
                    u++;
                }
            }
    }

    private static double reduceAcrossIntervals(double[] e, double[] f, double[] g, double[] h,
                                                int interval, double length,
                                                double[] sizes, double[] coalescent,
                                                int stateCount) {

        int offset = interval * stateCount;

        double sum = 0.0;
        for (int k = 0; k < stateCount; ++k) {
            sum += (e[offset + k] * e[offset + k] - f[offset + k] +
                    g[offset + k] * g[offset + k] - h[offset + k]) / sizes[k];
        }

        double logL = -length * sum / 4;

        double prob = coalescent[interval];
        if (prob != 0.0) {
            logL += Math.log(prob);
        }

        return logL;
    }

    private static void reduceWithinInterval(double[] e, double[] f, double[] g, double[] h,
                                             double[] partials,
                                             int startBuffer1, int startBuffer2,
                                             int endBuffer1, int endBuffer2,
                                             int interval, int stateCount) {

        interval *= stateCount;

        startBuffer1 *= stateCount;
        endBuffer1 *= stateCount;

        for (int i = 0; i < stateCount; ++i) {
            double startP = partials[startBuffer1 + i];
            e[interval + i] += startP;
            f[interval + i] += startP * startP;

            double endP = partials[endBuffer1 + i];
            g[interval + i] += endP;
            h[interval + i] += endP * endP;
        }

        if (startBuffer2 >= 0) {
            startBuffer2 *= stateCount;
            endBuffer2 *= stateCount;

            for (int i = 0; i < stateCount; ++i) {
                double startP = partials[startBuffer2 + i];
                e[interval + i] += startP;
                f[interval + i] += startP * startP;

                double endP = partials[endBuffer2 + i];
                g[interval + i] += endP;
                h[interval + i] += endP * endP;
            }
        }
    }
}

