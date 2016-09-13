/*
 * ProductChainSubstitutionModel.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.substmodel;

import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.inference.model.Model;
import dr.math.KroneckerOperation;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.*;
import java.util.logging.Logger;

//import dr.math.matrixAlgebra.Vector;

/**
 * @author Marc A. Suchard
 * @author Vladimir Minin
 *         <p/>
 *         A class for implementing a kronecker sum of CTMC models in BEAST using BEAGLE
 *         This work is supported by NSF grant 0856099
 *         <p/>
 *         O'Brien JD, Minin VN and Suchard MA (2009) Learning to count: robust estimates for labeled distances between
 *         molecular sequences. Molecular Biology and Evolution, 26, 801-814
 */
public class ProductChainSubstitutionModel extends BaseSubstitutionModel implements Citable {

    public ProductChainSubstitutionModel(String name, List<SubstitutionModel> baseModels) {
        this(name, baseModels, null);
    }

    public ProductChainSubstitutionModel(String name,
                                         List<SubstitutionModel> baseModels,
                                         List<SiteRateModel> rateModels) {
        this(name, baseModels, rateModels, false);
    }

    public ProductChainSubstitutionModel(String name,
                                         List<SubstitutionModel> baseModels,
                                         List<SiteRateModel> rateModels,
                                         boolean forceAverageModel) {
        super(name);

        this.baseModels = baseModels;
        this.rateModels = rateModels;
        this.forceAverageModel = forceAverageModel;
        numBaseModel = baseModels.size();

        if (numBaseModel == 0) {
            throw new RuntimeException("May not construct ProductChainSubstitutionModel with 0 base models");
        }

        if (rateModels != null) {
            for(SiteRateModel rateModel : rateModels) {
                if (rateModel.getCategoryCount() > 1) {
                    throw new RuntimeException("ProductChainSubstitutionModels with multiple categories not yet implemented");
                }
            }
        }

        List<FrequencyModel> freqModels = new ArrayList<FrequencyModel>();
        stateSizes = new int[numBaseModel];
        stateCount = 1;
        for (int i = 0; i < numBaseModel; i++) {
            freqModels.add(baseModels.get(i).getFrequencyModel());
            DataType dataType = baseModels.get(i).getDataType();
            stateSizes[i] = dataType.getStateCount();
            stateCount *= dataType.getStateCount();
            addModel(baseModels.get(i));
            addModel(rateModels.get(i));
        }

        pcFreqModel = new ProductChainFrequencyModel("pc",freqModels);
        addModel(pcFreqModel);

        String[] codeStrings = getCharacterStrings();

        dataType = new GeneralDataType(codeStrings);

        updateMatrix = true;

    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.SUBSTITUTION_MODELS;
    }

    @Override
    public String getDescription() {
        return "Product chain substitution model";
    }

    public List<Citation> getCitations() {
        return Collections.singletonList(CommonCitations.OBRIEN_2009_LEARNING);
    }

    public EigenDecomposition getEigenDecomposition() {
        synchronized (this) {
            if (updateMatrix) {
                computeKroneckerSumsAndProducts();
            }
        }
        return eigenDecomposition;
    }

    private String[] getCharacterStrings() {
        String[] strings = null;
        for (int i = numBaseModel - 1; i >= 0; i--) {
            strings = recursivelyAppendCharacterStates(baseModels.get(i).getDataType(), strings);
        }

        return strings;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        super.handleModelChangedEvent(model, object, index);
        // Propogate change to higher models
        fireModelChanged(model);
        averageModel = null;
    }

    private String[] recursivelyAppendCharacterStates(DataType dataType, String[] inSubStates) {

        String[] subStates = inSubStates;
        if (subStates == null) {
            subStates = new String[]{""};
        }

        final int previousStateCount = subStates.length;
        final int inStateCount = dataType.getStateCount();
        String[] states = new String[previousStateCount * inStateCount];

        for (int i = 0; i < inStateCount; i++) {
            String code = dataType.getCode(i);
            for (int j = 0; j < previousStateCount; j++) {
                states[i * previousStateCount + j] = code + subStates[j];
            }
        }
        return states;
    }

    public void getInfinitesimalMatrix(double[] out) {
        getEigenDecomposition(); // Updates rate matrix if necessary
        System.arraycopy(rateMatrix, 0, out, 0, stateCount * stateCount);
    }

    double getRateForModel(int index) {
        if (!forceAverageModel) {
            return rateModels.get(index).getRateForCategory(0);
        } else {
            double total = 0.0;
            for (int i = 0; i < rateModels.size(); ++i) {
                total += rateModels.get(i).getRateForCategory(0);
            }
            return total / rateModels.size();
        }
    }

    protected double[] scaleForProductChain(double[] in, int model) {
        if (rateModels == null) {
            return in;
        }
        final double scalar = getRateForModel(model);
        if (scalar == 1.0) {
            return in;
        }
        final int len = in.length;
        double[] out = new double[len];
        for (int i = 0; i < len; i++) {
            out[i] = scalar * in[i];
        }
        return out;
    }

    private SubstitutionProcess computeAverageModel() {

        return new SubstitutionProcess() {
            @Override
            public void getTransitionProbabilities(double distance, double[] matrix) {
                throw new RuntimeException("Should not be called");
            }

            @Override
            public EigenDecomposition getEigenDecomposition() {
                if (eigenDecomposition == null) {

//                    System.err.println("Statecount = " + stateSizes[0]);
//                    System.exit(-1);

                    double[][] mat = new double[stateSizes[0]][stateSizes[0]];
                    double[] vec = new double[stateSizes[0] * stateSizes[0]];
                    getInfinitesimalMatrix(vec);
                    for (int i = 0; i < stateSizes[0]; ++i) {
                        System.arraycopy(vec, i * stateSizes[0], mat[i], 0, stateSizes[0]);
                    }
                    eigenDecomposition = getDefaultEigenSystem(stateSizes[0]).decomposeMatrix(mat);
                }
                return eigenDecomposition;
            }

            @Override
            public FrequencyModel getFrequencyModel() {
                throw new RuntimeException("Should not be called");
            }

            @Override
            public void getInfinitesimalMatrix(double[] matrix) {
                if (averageMatrix == null) {
                    final int dim = matrix.length;
                    averageMatrix = new double[dim];

                    double[][] allMatrices = new double[baseModels.size()][dim];
                    for (int i = 0; i < baseModels.size(); ++i) {
                        baseModels.get(i).getInfinitesimalMatrix(allMatrices[i]);
                    }

                    for (int i = 0; i < dim; ++i) {
                        double total = 0.0;
                        for (int j = 0; j < baseModels.size(); ++j) {
                            total += allMatrices[j][i];
                        }
                        total /= baseModels.size();
                        averageMatrix[i] = total;
                    }
                }
//                System.err.println("averageMatrix.length " + averageMatrix.length);
//                System.err.println("matrix.length " + matrix.length);
                System.arraycopy(averageMatrix, 0, matrix, 0, averageMatrix.length);
            }

            @Override
            public DataType getDataType() {
                throw new RuntimeException("Should not be called");
            }

            @Override
            public boolean canReturnComplexDiagonalization() {
                throw new RuntimeException("Should not be called");
            }

            private double[] averageMatrix = null;
            private EigenDecomposition eigenDecomposition = null;
        };
    }

//    private SubstitutionProcess getAverageModel() {
//        if (!forceAverageModel) {
//            throw new RuntimeException("Error getting averaged model with non-averaged product chain");
//        }
//        if (averageModel == null) {
//            averageModel = computeAverageModel();
//        }
//        return averageModel;
//    }

    private SubstitutionProcess getBaseModel(int index) {
        if (!forceAverageModel) {
            return baseModels.get(index);
        } else {
            if (averageModel == null) {
                averageModel = computeAverageModel();
            }
            return averageModel;
        }
    }

    private void computeKroneckerSumsAndProducts() {

//        if (forceAverageModel) {
//            if (averageModel == null) {
//                averageModel = computeAverageModel();
//            }
//            eigenDecomposition = averageModel.getEigenDecomposition();
//            if (rateMatrix == null) {
//                rateMatrix = new double[stateCount * stateCount];
//            }
//            averageModel.getInfinitesimalMatrix(rateMatrix);
//            updateMatrix = false;
//
//            System.err.println("lambda = " + new Vector(rateMatrix));
//
//            return;
//        }

        int currentStateSize = stateSizes[0];
        double[] currentRate = new double[currentStateSize * currentStateSize];
        getBaseModel(0).getInfinitesimalMatrix(currentRate);
        currentRate = scaleForProductChain(currentRate, 0);

        EigenDecomposition currentED = getBaseModel(0).getEigenDecomposition();
        double[] currentEval = scaleForProductChain(currentED.getEigenValues(), 0);
        double[] currentEvec = currentED.getEigenVectors();
        double[] currentIevcT = transpose(currentED.getInverseEigenVectors(), currentStateSize);

        for (int i = 1; i < numBaseModel; i++) {
            SubstitutionProcess nextModel = getBaseModel(i);
            int nextStateSize = stateSizes[i];
            double[] nextRate = new double[nextStateSize * nextStateSize];
            nextModel.getInfinitesimalMatrix(nextRate);
            nextRate = scaleForProductChain(nextRate, i);
            currentRate = KroneckerOperation.sum(currentRate, currentStateSize, nextRate, nextStateSize);

            EigenDecomposition nextED = nextModel.getEigenDecomposition();
            double[] nextEval = scaleForProductChain(nextED.getEigenValues(), i);
            double[] nextEvec = nextED.getEigenVectors();
            double[] nextIevcT = transpose(nextED.getInverseEigenVectors(), nextStateSize);

            currentEval = KroneckerOperation.sum(currentEval, nextEval);

            currentEvec = KroneckerOperation.product(
                    currentEvec, currentStateSize, currentStateSize,
                    nextEvec, nextStateSize, nextStateSize);

            currentIevcT = KroneckerOperation.product(
                    currentIevcT, currentStateSize, currentStateSize,
                    nextIevcT, nextStateSize, nextStateSize);
            currentStateSize *= nextStateSize;

        }

        rateMatrix = currentRate;

        eigenDecomposition = new EigenDecomposition(
                currentEvec,
                transpose(currentIevcT, currentStateSize),
                currentEval);
        updateMatrix = false;
    }

//   private static void printSquareMatrix(double[] A, int dim) {
//        double[] row = new double[dim];
//        for (int i = 0; i < dim; i++) {
//            System.arraycopy(A, i * dim, row, 0, dim);
//            System.err.println(new Vector(row));
//        }
//    }

    // transposes a square matrix
    private static double[] transpose(double[] mat, int dim) {
        double[] out = new double[dim * dim];
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                out[j * dim + i] = mat[i * dim + j];
            }
        }
        return out;
    }

    public FrequencyModel getFrequencyModel() {
        return pcFreqModel;
    }

    protected void frequenciesChanged() {
        // Do nothing
    }

    protected void ratesChanged() {
        // Do nothing
    }

    protected void setupRelativeRates(double[] rates) {
        // Do nothing
    }

    protected final int numBaseModel;
    protected final List<SubstitutionModel> baseModels;
    protected final List<SiteRateModel> rateModels;
    protected final int[] stateSizes;
    protected final ProductChainFrequencyModel pcFreqModel;
    protected double[] rateMatrix = null;
    private final boolean forceAverageModel;
    private SubstitutionProcess averageModel = null;

}
