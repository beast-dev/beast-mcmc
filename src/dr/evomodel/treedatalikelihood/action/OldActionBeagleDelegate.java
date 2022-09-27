/*
 * ActionBeagleDelegate.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.action;

import beagle.Beagle;
import beagle.InstanceDetails;
import dr.evolution.tree.Tree;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import org.newejml.data.DMatrixRMaj;
import org.newejml.data.DMatrixSparseCSC;
import org.newejml.dense.row.CommonOps_DDRM;
import org.newejml.dense.row.NormOps_DDRM;
import org.newejml.sparse.csc.CommonOps_DSCC;
import org.newejml.sparse.csc.misc.ImplCommonOps_DSCC;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.IntStream;

/**
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class OldActionBeagleDelegate implements Beagle {

    public static final boolean DEBUG = false;
    public static final boolean SCALING = true;
    public static final int SCALING_FACTOR_COUNT = 254;
    public static final int SCALING_FACTOR_OFFSET = 126;
    private static final int SCALING_EXPONENT_THRESHOLD = 2;
    protected final int tipCount;
    protected final int partialsBufferCount;
    protected final int stateCount;
    protected final int patternCount;
    protected final int categoryCount;
    protected int partialsSize;
    protected int matrixSize;
    protected double[][] cMatrices;
    protected double[] stateFrequencies;
    protected double[] categoryRates;
    protected double[] categoryWeights;
    protected DMatrixRMaj patternWeights;
    protected DMatrixRMaj[][] partials;
    protected int[][] scalingFactorCounts;
    protected double[][] matrices;
    double[] tmpPartials;
    protected DMatrixRMaj[] scalingFactors;
    protected double[] logScalingFactors;
    protected DMatrixSparseCSC[] instantaneousMatrices;
    private final ActionEvolutionaryProcessDelegate evolutionaryProcessDelegate;
    private final PartialsRescalingScheme rescalingScheme;
    private final boolean[] activeScalingFactors;
    private final DMatrixRMaj[] autoScalingBuffers;
    private final Tree treeModel;

    public OldActionBeagleDelegate(Tree treeModel,
                                   int partialsBufferCount,
                                   int patternCount,
                                   int stateCount,
                                   int categoryCount,
                                   int partialsSize,
                                   int scaleBufferCount,
                                   PartialsRescalingScheme rescalingScheme,
                                   ActionEvolutionaryProcessDelegate evolutionaryProcessDelegate) {
        this.treeModel = treeModel;
        this.tipCount = treeModel.getExternalNodeCount();
        this.partialsBufferCount = partialsBufferCount + 1;
        this.patternCount = patternCount;
        this.stateCount = stateCount;
        this.categoryCount = categoryCount;
        this.partialsSize = partialsSize;
        this.categoryWeights = new double[categoryCount];
        this.categoryRates = new double[categoryCount];
        this.patternWeights = new DMatrixRMaj(1, patternCount);
        this.stateFrequencies = new double[stateCount];
        this.partials = new DMatrixRMaj[partialsBufferCount][categoryCount];
        this.scalingFactors = new DMatrixRMaj[scaleBufferCount];
        this.autoScalingBuffers = new DMatrixRMaj[partialsBufferCount];
        this.activeScalingFactors = new boolean[partialsBufferCount];
        for (int i = 0; i < partialsBufferCount; i++) {
            for (int j = 0; j < categoryCount; j++) {
                partials[i][j] = new DMatrixRMaj(patternCount, stateCount);
            }
            autoScalingBuffers[i] = new DMatrixRMaj(1, patternCount);
        }
        for (int i = 0; i < scaleBufferCount; i++) {
            scalingFactors[i] = new DMatrixRMaj(1, patternCount);
        }
        this.instantaneousMatrices = instantaneousMatrices;
        this.evolutionaryProcessDelegate = evolutionaryProcessDelegate;
        this.rescalingScheme = rescalingScheme;
        this.identity = CommonOps_DSCC.identity(stateCount);
        this.A = new DMatrixSparseCSC(stateCount, stateCount);
        this.rescaleRows = IntStream.range(0, patternCount).toArray();
        this.tmpLogScalingFactors = new DMatrixRMaj(1, patternCount);
    }

    @Override
    public void finalize() throws Throwable {

    }

    @Override
    public void setCPUThreadCount(int i) {

    }

    @Override
    public void setPatternWeights(double[] doubles) {
        patternWeights.setData(doubles);
    }

    @Override
    public void setPatternPartitions(int i, int[] ints) {
        throw new UnsupportedOperationException("Not implemented yet");
    }

    @Override
    public void setTipStates(int i, int[] ints) {
        double[] translatedTipPartials = new double[ints.length * stateCount * categoryCount];
        for (int category = 0; category < categoryCount; category++) {
            for (int j = 0; j < ints.length; j++) {
                if (ints[j] < stateCount) {
                    translatedTipPartials[category * stateCount * patternCount + stateCount * j + ints[j]] = 1.0;
                } else {
                    Arrays.fill(translatedTipPartials, category * stateCount * patternCount + stateCount * j,
                            category * stateCount * patternCount + stateCount * (j + 1), 1.0);
                }
            }
        }
        setTipPartials(i, translatedTipPartials);
    }

    @Override
    public void getTipStates(int i, int[] ints) {
        throw new UnsupportedOperationException("Not implemented yet");
    }

    @Override
    public void setTipPartials(int i, double[] doubles) {
        assert(i >= 0 && i < tipCount);
        assert(doubles.length == partialsSize);
        if (partials[i] == null) {
            for (int j = 0; j < categoryCount; j++) {
                partials[i][j] = new DMatrixRMaj(patternCount, stateCount);
            }
        }

        double[] categoryData = new double[patternCount * stateCount];
        for (int j = 0; j < categoryCount; j++) {
            System.arraycopy(doubles, j * patternCount * stateCount, categoryData, 0, patternCount * stateCount);
            partials[i][j].set(patternCount, stateCount, true, categoryData);
        }

    }
//
//    private int getPartialIndex(int i) {  // TODO: use BufferIndexHelper by XJ
//        return i * categoryCount * patternCount * stateCount * partialsPerNode;
//    }

    @Override
    public void setRootPrePartials(int[] ints, int[] ints1, int i) {
        throw new RuntimeException("Should not be called with action likelihood");
    }

    @Override
    public void setPartials(int i, double[] doubles) { //TODO: check for double buffering by XJ
        assert this.partials != null;
        assert(doubles.length == partialsSize);

        if (partials[i] == null) {
            for (int j = 0; j < categoryCount; j++) {
                partials[i][j] = new DMatrixRMaj(patternCount, stateCount);
            }
        }

        double[] categoryData = new double[patternCount * stateCount];
        for (int j = 0; j < categoryCount; j++) {
            System.arraycopy(doubles, j * patternCount * stateCount, categoryData, 0, patternCount * stateCount);
            partials[i][j].set(patternCount, stateCount, true, categoryData);
        }

    }

    @Override
    public void getPartials(int i, int i1, double[] doubles) {
        for (int j = 0; j < categoryCount; j++) {
            System.arraycopy(partials[i][j].getData(), 0, doubles, j * patternCount * stateCount, patternCount * stateCount);
        }
    }

    @Override
    public void getLogScaleFactors(int i, double[] doubles) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public void setEigenDecomposition(int i, double[] doubles, double[] doubles1, double[] doubles2) {
        throw new RuntimeException("Should not be called with action likelihood");
    }

    @Override
    public void setStateFrequencies(int i, double[] doubles) {
        System.arraycopy(doubles, 0, this.stateFrequencies, 0, this.stateCount);
    }

    @Override
    public void setCategoryWeights(int i, double[] doubles) {
        System.arraycopy(doubles, 0, this.categoryWeights, 0, this.categoryCount);
    }

    @Override
    public void setCategoryRates(double[] doubles) {
        System.arraycopy(doubles, 0, this.categoryRates, 0, this.categoryRates.length);
    }

    @Override
    public void setCategoryRatesWithIndex(int i, double[] doubles) {
        throw new UnsupportedOperationException("setCategoryRatesWithIndex not implemented in ActionBeagleDelegate.");
    }

    @Override
    public void convolveTransitionMatrices(int[] ints, int[] ints1, int[] ints2, int i) {
        throw new UnsupportedOperationException("convolveTransitionMatrices not implemented in ActionBeagleDelegate.");
    }

    @Override
    public void addTransitionMatrices(int[] ints, int[] ints1, int[] ints2, int i) {
        throw new UnsupportedOperationException("addTransitionMatrices not implemented in ActionBeagleDelegate.");
    }

    @Override
    public void transposeTransitionMatrices(int[] ints, int[] ints1, int i) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void updateTransitionMatrices(int i, int[] ints, int[] ints1, int[] ints2, double[] doubles, int i1) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void updateTransitionMatricesWithMultipleModels(int[] ints, int[] ints1, int[] ints2, int[] ints3, int[] ints4, double[] doubles, int i) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void setTransitionMatrix(int i, double[] doubles, double v) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void setDifferentialMatrix(int i, double[] doubles) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void getTransitionMatrix(int i, double[] doubles) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void updatePrePartials(int[] ints, int i, int i1) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void updatePrePartialsByPartition(int[] ints, int i) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void calculateEdgeDerivative(int[] ints, int[] ints1, int i, int[] ints2, int[] ints3, int i1, int i2, int i3, int[] ints4, int i4, double[] doubles, double[] doubles1) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void calculateEdgeDifferentials(int[] ints, int[] ints1, int[] ints2, int[] ints3, int i, double[] doubles, double[] doubles1, double[] doubles2) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void calculateCrossProductDifferentials(int[] ints, int[] ints1, int[] ints2, int[] ints3, double[] doubles, int i, double[] doubles1, double[] doubles2) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void calculateRootLogLikelihoods(int[] ints, int[] ints1, int[] ints2, int[] ints3, int i, double[] doubles) {
        DMatrixRMaj colSums = new DMatrixRMaj(1, patternCount, true, new double[patternCount]);
        final int cumulateScaleBufferIndex = ints3[0];
        for (int j = 0; j < categoryCount; j++) {
            DMatrixRMaj rootPartial = CommonOps_DDRM.transpose(partials[ints[0]][j], null);

            CommonOps_DDRM.scale(categoryWeights[j], rootPartial);

            CommonOps_DDRM.multRows(stateFrequencies, rootPartial);
            DMatrixRMaj singleCategoryColSum = CommonOps_DDRM.sumCols(rootPartial, null);
            CommonOps_DDRM.add(singleCategoryColSum, colSums, colSums);
        }
        DMatrixRMaj siteLogL = CommonOps_DDRM.elementLog(colSums, null);

        if (cumulateScaleBufferIndex >= 0) {
            DMatrixRMaj cumulativeScaleFactors = scalingFactors[cumulateScaleBufferIndex];
            CommonOps_DDRM.add(siteLogL, cumulativeScaleFactors, siteLogL);
        }

        CommonOps_DDRM.elementMult(siteLogL, patternWeights, siteLogL);
        doubles[0] = CommonOps_DDRM.elementSum(siteLogL);
    }

    final private int operationSize = 7;

    @Override
    public void updatePartials(int[] ints, int i, int i1) {
        for (int operation = 0; operation < i; operation ++) {
            final int destinationPartialIndex = ints[operation * operationSize];
            final int writeScaleIndex = ints[operation * operationSize + 1];
            final int readScaleIndex = ints[operation * operationSize + 2];
            final int firstChildPartialIndex = ints[operation * operationSize + 3];
            final int firstChildSubstitutionMatrixIndex = ints[operation * operationSize + 4];
            final int secondChildPartialIndex = ints[operation * operationSize + 5];
            final int secondChildSubstitutionMatrixIndex = ints[operation * operationSize + 6];

            // consistent with Beagle cases for easier future reference
            int rescale = Beagle.NONE;
            if (writeScaleIndex != Beagle.NONE) {
                rescale = 1;
            } else if (readScaleIndex != Beagle.NONE) {
                rescale = 0;
            }

            for (int j = 0; j < categoryCount; j++) {
                DMatrixRMaj leftPartial = partials[firstChildPartialIndex][j];
                DMatrixRMaj rightPartial = partials[secondChildPartialIndex][j];

                DMatrixSparseCSC leftGeneratorMatrix = evolutionaryProcessDelegate.getScaledInstantaneousMatrix(firstChildSubstitutionMatrixIndex, categoryRates[j]);
                DMatrixSparseCSC rightGeneratorMatrix = evolutionaryProcessDelegate.getScaledInstantaneousMatrix(secondChildSubstitutionMatrixIndex, categoryRates[j]);

                DMatrixRMaj parentLeftPostPartial = simpleAction(leftGeneratorMatrix, leftPartial);
                DMatrixRMaj parentRightPostPartial = simpleAction(rightGeneratorMatrix, rightPartial);

                CommonOps_DDRM.elementMult(parentLeftPostPartial, parentRightPostPartial, partials[destinationPartialIndex][j]);

                if (rescale == 0) {
                    CommonOps_DDRM.elementDiv(partials[destinationPartialIndex][j], scalingFactors[readScaleIndex], partials[destinationPartialIndex][j]);
                }
            }
            if (rescale == 1 && i1 != Beagle.NONE) {
                rescalePartials(partials[destinationPartialIndex], scalingFactors[writeScaleIndex], scalingFactors[i1]);
            }

        }
    }

    private final int[] rescaleRows;
    private final DMatrixRMaj tmpLogScalingFactors;

    private void rescalePartials(DMatrixRMaj[] destPartials, DMatrixRMaj scalingFactors, DMatrixRMaj cumulativeScaleBuffer) {
        DMatrixRMaj categoryCombined = new DMatrixRMaj(patternCount, stateCount * categoryCount);
        final int rowSize = patternCount;
        final int colSize = stateCount;
        for (int j = 0; j < categoryCount; j++) {
            CommonOps_DDRM.insert(destPartials[j], categoryCombined,
                    rescaleRows, rowSize,
                    IntStream.range(j * colSize, j * colSize + rowSize).toArray(), colSize);
        }
        CommonOps_DDRM.maxRows(categoryCombined, scalingFactors);
        CommonOps_DDRM.transpose(scalingFactors);
        CommonOps_DDRM.elementLog(scalingFactors, tmpLogScalingFactors);
        CommonOps_DDRM.add(cumulativeScaleBuffer, tmpLogScalingFactors, cumulativeScaleBuffer);

        for (int j = 0; j < categoryCount; j++) {
            CommonOps_DDRM.divideRows(scalingFactors.getData(), destPartials[j]);
        }
    }

    private final DMatrixSparseCSC identity;
    private final DMatrixSparseCSC A;

    private DMatrixRMaj simpleAction(DMatrixSparseCSC matrix, DMatrixRMaj partials) {

        // Algorithm 3.2 in Al-Mohy and Higham (2011), balance = false, t = 1.0
        assert(matrix.numCols == matrix.numRows);

        DMatrixRMaj B = CommonOps_DDRM.transpose(partials, null);
        assert(matrix.numCols == B.numRows);

        final double tol = Math.pow(2, -53);
        final double t = 1.0;
        final int nCol = B.getNumCols();
        final double mu = CommonOps_DSCC.trace(matrix)/((double) matrix.numCols);

//        CommonOps_DSCC.add(1.0, matrix, -mu, identity, A, null, null);
        ImplCommonOps_DSCC.add(1.0, matrix, -mu, identity, A, null, null);

        final double A1Norm = normP1(A);
        int m, s;
        if (t * A1Norm == 0.0) {
            m = 0;
            s = 1;
        } else {
            TaylorSeriesStatistics taylorSeriesStatistics = new TaylorSeriesStatistics(A1Norm, A, t, nCol);
            m = taylorSeriesStatistics.getM();
            s = taylorSeriesStatistics.getS();
        }

        DMatrixRMaj F = B.copy();
        final double eta = Math.exp(t * mu / ((double) s));
        double c1,c2;
        for (int i = 0; i < s; i++) {
            c1 = NormOps_DDRM.normPInf(B);
            for (int j = 1; j < m + 1; j++) {
                B = CommonOps_DSCC.mult(A, B, null);
                CommonOps_DDRM.scale(t / ((double) s * j), B);
                c2 = NormOps_DDRM.normPInf(B);
                CommonOps_DDRM.add(F, B, F);
                if (c1 + c2 <= tol * NormOps_DDRM.normPInf(F)) {
                    break;
                }
                c1 = c2;
            }
            CommonOps_DDRM.scale(eta, F);
            B.setTo(F);
        }

        return CommonOps_DDRM.transpose(F, null);
    }

    private class TaylorSeriesStatistics {
        private int m;
        private int s;
        private final int mMax = 55;
        private final double pMax;

        TaylorSeriesStatistics(double A1Norm, DMatrixSparseCSC A, double t, int nCol) {
            this.thetaConstants = new HashMap<>();
            this.d = new HashMap<>();
            this.powerMatrices = new HashMap<>();
            powerMatrices.put(1, A);
            highestPower = 1;
            this.pMax = getPMax(mMax);
            cachePowerSeries((int) pMax);
            setThetaConstants();


            if (t * A1Norm == 0.0) {
                m = 0;
                s = 1;
            } else {
                setStatistics(A1Norm, A, t, nCol);
            }
        }

        public int getS() {
            return s;
        }

        public int getM() {
            return m;
        }

        private void setStatistics(double A1Norm, DMatrixSparseCSC A, double t, int nCol) {
            // Code fragment 3.1 in Al-Mohy and Higham
            if (t != 1.0) {
                System.err.println("Not yet implemented!");
//                CommonOps_DSCC.scale(t, A, null);  // TODO: make sure the where the scaling happens, should be at updateMatrix()
            }
            int bestM = Integer.MAX_VALUE;
            int bestS = Integer.MAX_VALUE;
            if (conditionFragment313(A1Norm, nCol, mMax)) {
                for (int thisM : thetaConstants.keySet()) {
                    final double thisS = Math.ceil(A1Norm/thetaConstants.get(thisM));
                    if (bestM == Integer.MAX_VALUE || ((double) thisM) * thisS < bestM * bestS) {
                        bestS = (int) thisS;
                        bestM = thisM;
                    }
                }
                this.s = bestS;
            } else {
                for (int p = 2; p < pMax; p++) {
                    for (int thisM = p * (p - 1) - 1; thisM < mMax + 1; thisM++) {
                        if (thetaConstants.containsKey(thisM)) {
                            // part of equation 3.10
                            final double thisS = Math.ceil(getAlpha(p) / thetaConstants.get(thisM));
                            if (bestM == Integer.MAX_VALUE || ((double) thisM) * thisS < bestM * bestS) {
                                bestS = (int) thisS;
                                bestM = thisM;
                            }
                        }
                    }
                }
                this.s = Math.max(bestS, 1);
            }
            this.m = bestM;
        }

        private double getAlpha(int p) {
            // equation 3.7 in Al-Mohy and Higham
            return Math.max(getDValue(p), getDValue(p + 1));
        }



        private double getDValue(int p) {
            // equation 3.7 in Al-Mohy and Higham
            if (!d.containsKey(p)) {
                cachePowerSeries(p);
                DMatrixSparseCSC powerPMatrix = powerMatrices.get(p);
                d.put(p, Math.pow(normP1(powerPMatrix), 1.0 / ((double) p)));
            }
            return d.get(p);
        }

        private void cachePowerSeries(int p) {
            if (highestPower < p) {
                for (int i = highestPower; i < p; i++) {
                    DMatrixSparseCSC currentPowerMatrix = powerMatrices.get(highestPower);
                    DMatrixSparseCSC nextPowerMatrix = CommonOps_DSCC.mult(currentPowerMatrix, powerMatrices.get(1), null);
                    powerMatrices.put(i + 1, nextPowerMatrix);
                    highestPower++;
                }
            }
        }

        private final Map<Integer, Double> d;
        private final Map<Integer, DMatrixSparseCSC> powerMatrices;
        private int highestPower;

        private boolean conditionFragment313(double A1Norm, int nCol, int mMax) {
            // using l = 1 as in equation 3.13
            final double theta = thetaConstants.get(mMax);

            return A1Norm <= 2.0 * theta / ((double) nCol * mMax) * pMax * (pMax + 3);
        }
        private double getPMax(int mMax) {
            // pMax is the largest positive integer such that p*(p-1) <= mMax + 1
            final double pMax = Math.floor(0.5 + 0.5 * Math.sqrt(5.0 + 4.0 * mMax));

            return pMax;
        }


        private final Map<Integer, Double> thetaConstants;
        private void setThetaConstants() {
            //The first 30 values are from table A.3 of  Computing Matrix Functions.
            // For double precision, tol = 2^(-53)
            // TODO: maybe calculate this
            thetaConstants.put(1, 2.29E-16);
            thetaConstants.put(2, 2.58E-8);
            thetaConstants.put(3, 1.39E-5);
            thetaConstants.put(4, 3.40E-4);
            thetaConstants.put(5, 2.40E-3);
            thetaConstants.put(6, 9.07E-3);
            thetaConstants.put(7, 2.38E-2);
            thetaConstants.put(8, 5.00E-2);
            thetaConstants.put(9, 8.96E-2);
            thetaConstants.put(10, 1.44E-1);
            thetaConstants.put(11, 2.14E-1);
            thetaConstants.put(12, 3.00E-1);
            thetaConstants.put(13, 4.00E-1);
            thetaConstants.put(14, 5.14E-1);
            thetaConstants.put(15, 6.41E-1);
            thetaConstants.put(16, 7.81E-1);
            thetaConstants.put(17, 9.31E-1);
            thetaConstants.put(18, 1.09);
            thetaConstants.put(19, 1.26);
            thetaConstants.put(20, 1.44);
            thetaConstants.put(21, 1.62);
            thetaConstants.put(22, 1.82);
            thetaConstants.put(23, 2.01);
            thetaConstants.put(24, 2.22);
            thetaConstants.put(25, 2.43);
            thetaConstants.put(26, 2.64);
            thetaConstants.put(27, 2.86);
            thetaConstants.put(28, 3.08);
            thetaConstants.put(29, 3.31);
            thetaConstants.put(30, 3.54);
            //The rest are from table 3.1 of Computing the Action of the Matrix Exponential.
            thetaConstants.put(35, 4.7);
            thetaConstants.put(40, 6.0);
            thetaConstants.put(45, 7.2);
            thetaConstants.put(50, 8.5);
            thetaConstants.put(55, 9.9);
        }

    }

    private double normP1( DMatrixSparseCSC A ) {
        DMatrixSparseCSC absA = A.copy();
        for (int i = 0; i < absA.numCols; i++) {
            absA.set(i, i, Math.abs(absA.get(i, i)));
        }
        DMatrixRMaj colSums = CommonOps_DSCC.sumCols(absA, null);

        return CommonOps_DDRM.elementMax(colSums);
    }

    @Override
    public void updatePartialsByPartition(int[] ints, int i) {
        throw new RuntimeException("Not yet implemented");
    }

    private static final double kLn2 = Math.log(2);

    @Override
    public void accumulateScaleFactors(int[] scalingIndices, int count, int cumulativeScalingIndex) {
        if (rescalingScheme == PartialsRescalingScheme.AUTO) {
            DMatrixRMaj cumulativeScaleBuffer = scalingFactors[0];
            cumulativeScaleBuffer.fill(0.0);
            for (int j = 0; j < count; j++) {
                final int sIndex = scalingIndices[j] - tipCount;
                if (activeScalingFactors[sIndex]) {
                    DMatrixRMaj scaleBuffer = autoScalingBuffers[sIndex];
                    CommonOps_DDRM.add(scaleBuffer, kLn2, cumulativeScaleBuffer);
                }
            }
            throw new RuntimeException("Not tested yet.");
        } else {
            DMatrixRMaj cumulativeScaleBuffer = scalingFactors[cumulativeScalingIndex];
            for (int j = 0; j < count; j++) {
                DMatrixRMaj scaleBuffer = scalingFactors[scalingIndices[j]];
                //TODO: enable BEAGLE_FLAG_SCALERS_LOG switch
                CommonOps_DDRM.add(cumulativeScaleBuffer, scaleBuffer, cumulativeScaleBuffer);
            }
        }
    }

    @Override
    public void accumulateScaleFactorsByPartition(int[] ints, int i, int i1, int i2) {
        throw new RuntimeException("Not yet impelmented!");
    }

    @Override
    public void removeScaleFactors(int[] ints, int i, int i1) {

    }

    @Override
    public void removeScaleFactorsByPartition(int[] ints, int i, int i1, int i2) {
        throw new RuntimeException("Not yet impelmented!");
    }

    @Override
    public void copyScaleFactors(int i, int i1) {

    }

    @Override
    public void resetScaleFactors(int i) {
        scalingFactors[i].fill(0.0);
    }

    @Override
    public void resetScaleFactorsByPartition(int i, int i1) {

    }

    @Override
    public void calculateRootLogLikelihoodsByPartition(int[] ints, int[] ints1, int[] ints2, int[] ints3, int[] ints4, int i, int i1, double[] doubles, double[] doubles1) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public void getSiteLogLikelihoods(double[] doubles) {

    }

    @Override
    public InstanceDetails getDetails() {
        return null;
    }
}
