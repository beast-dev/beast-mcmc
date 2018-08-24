/*
 * MarkovModulatedSubstitutionModel.java
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

import dr.evolution.datatype.DataType;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Marc A. Suchard
 */
public class MarkovModulatedSubstitutionModel extends ComplexSubstitutionModel implements Citable, Loggable {

    private List<SubstitutionModel> baseModels;
    private final int numBaseModel;
    private final int baseStateCount;
    //    private final int stateCount;
    private final Parameter switchingRates;

    private static final boolean IGNORE_RATES = false;
    private static final boolean DEBUG = false;
    private static final boolean NEW_STORE_RESTORE = true;

    private final double[] baseMatrix;
    private Parameter rateScalar;

    private boolean birthDeathModel;
    private boolean geometricRates;

    private final SiteRateModel gammaRateModel;

    private EigenDecomposition storedEigenDecomposition;
    private boolean storedUpdateMatrix;

    public MarkovModulatedSubstitutionModel(String name,
                                            List<SubstitutionModel> baseModels,
                                            Parameter switchingRates,
                                            DataType dataType,
                                            EigenSystem eigenSystem) {
        this(name, baseModels, switchingRates, dataType, eigenSystem, null, false, null);
    }

    public MarkovModulatedSubstitutionModel(String name,
                                            List<SubstitutionModel> baseModels,
                                            Parameter switchingRates,
                                            DataType dataType,
                                            EigenSystem eigenSystem,
                                            Parameter rateScalar,
                                            boolean geometricRates,
                                            SiteRateModel gammaRateModel) {
//        super(name, dataType, null, eigenSystem);
        super(name, dataType, null, null);

        this.baseModels = baseModels;
        numBaseModel = baseModels.size();

        if (numBaseModel == 0) {
            throw new RuntimeException("May not construct MarkovModulatedSubstitutionModel with 0 base models");
        }

        this.switchingRates = switchingRates;
        addVariable(switchingRates);

        if (switchingRates.getDimension() != 2 * (numBaseModel - 1)
                && switchingRates.getDimension() != numBaseModel * (numBaseModel - 1)
                ) {
            throw new RuntimeException("Wrong switching rate dimensions");
        }

        List<FrequencyModel> freqModels = new ArrayList<FrequencyModel>();
        int stateSizes = 0;

        baseStateCount = baseModels.get(0).getFrequencyModel().getFrequencyCount();
        baseMatrix = new double[baseStateCount * baseStateCount];

        for (int i = 0; i < numBaseModel; i++) {
            addModel(baseModels.get(i));
            freqModels.add(baseModels.get(i).getFrequencyModel());
            addModel(baseModels.get(i).getFrequencyModel());
            DataType thisDataType = baseModels.get(i).getDataType();
            stateSizes += thisDataType.getStateCount();
        }

        // This constructor also checks that all models have the same base stateCount
        freqModel = new MarkovModulatedFrequencyModel("mm", freqModels, switchingRates);
        addModel(freqModel);

        if (stateCount != stateSizes) {
            throw new RuntimeException("Incompatible state counts in " + getModelName() + " (currently: " + stateCount + "). Models add up to " + stateSizes + ".");
        }

        birthDeathModel = true;
        this.geometricRates = geometricRates;
        // Check switching rate dimension
        if (numBaseModel > 1) {
            if (switchingRates.getDimension() != 2 * (numBaseModel - 1)) {
                birthDeathModel = false;
//                throw new RuntimeException("Wrong dimension of switching rates in MarkovModulatedSubstitutionModel " + switchingRates.getDimension() + " " + 2 * (numBaseModel - 1) + " " + numBaseModel);
            }
        }

        if (gammaRateModel != null) {
            addModel(gammaRateModel);

            if (gammaRateModel.getCategoryCount() != numBaseModel && numBaseModel % gammaRateModel.getCategoryCount()  != 0) {
                throw new RuntimeException("Wrong discretized gamma dimension");
            }
        }
        this.gammaRateModel = gammaRateModel;

        if (rateScalar != null) {
            addVariable(rateScalar);

            if (rateScalar.getDimension() != 1 && rateScalar.getDimension() != numBaseModel) {
                throw new RuntimeException("Wrong rate scalar dimensions");
            }
        }
        this.rateScalar = rateScalar;

        setDoNormalization(false);

        updateMatrix = true;

        Logger.getLogger("dr.app.beagle").info("\tConstructing a Markov-modulated Markov chain substitution model with " + stateCount + " states;  please cite:\n"
                + Citable.Utils.getCitationString(this));
    }

    public int getNumBaseModel() {
        return numBaseModel;
    }

    public double getModelRateScalar(int model) {
        if (gammaRateModel != null) {
            model = model % gammaRateModel.getCategoryCount();
            if (DEBUG) {
                System.err.println("M" + model + " = " + gammaRateModel.getRateForCategory(model));
            }
            return gammaRateModel.getRateForCategory(model);
        }
        if (rateScalar == null) {
            return 1.0;
        } else {
            if (rateScalar.getDimension() == 1) {
                return rateScalar.getParameterValue(0);
            } else {
                return rateScalar.getParameterValue(model);
            }
        }
    }

    protected void storeState() {
        if (DEBUG) {
            System.err.println("MMSM.sS");
        }

        if (NEW_STORE_RESTORE) {
            if (eigenDecomposition != null) {
                storedEigenDecomposition = eigenDecomposition.copy();
            }
            storedUpdateMatrix = updateMatrix;
        } else {
            super.storeState();
        }
    }

    protected void restoreState() {
        if (DEBUG) {
            System.err.println("MMSM.rS");
        }

        if (NEW_STORE_RESTORE) {
            EigenDecomposition tmp = storedEigenDecomposition;
            storedEigenDecomposition = eigenDecomposition;
            eigenDecomposition = tmp;

            updateMatrix = storedUpdateMatrix;
        } else {
            super.restoreState();
        }
    }

    protected void setupQMatrix(double[] rates, double[] pi, double[][] matrix) {

//        System.err.println("MMSM.sQM");
        // Zero matrix
        for (int i = 0; i < matrix.length; ++i) {
            Arrays.fill(matrix[i], 0.0);
        }
        // Set the instantaneous rate matrix
        for (int m = 0; m < numBaseModel; ++m) {
            final int offset = m * baseStateCount;
            baseModels.get(m).getInfinitesimalMatrix(baseMatrix);
            if (DEBUG) {
                System.err.println("m " + m + " : " + new dr.math.matrixAlgebra.Vector(baseMatrix));
            }
            final double rateScalar = getModelRateScalar(m);
            int k = 0;
            for (int i = 0; i < baseStateCount; i++) {
                for (int j = 0; j < baseStateCount; j++) {
                    matrix[offset + i][offset + j] = rateScalar * baseMatrix[k];
                    k++;
                }
            }
        }

        // Add switching rates to matrix
        if (!IGNORE_RATES && numBaseModel > 1) {
            double[] swRates = switchingRates.getParameterValues();
            if (DEBUG) {
                System.err.print("Switching rates: ");
                for (int i = 0; i < swRates.length; i++) {
                    System.err.print(swRates[i] + " ");
                }
                System.err.println();
            }
            double totalRate = 0.0;
            for (double rate : swRates) {
                totalRate += rate;
            }
            int sw = 0;
            for (int g = 0; g < numBaseModel; ++g) {
                for (int h = 0; h < numBaseModel; ++h) { // from g -> h
                    boolean valid = birthDeathModel ? Math.abs(g - h) == 1 : g != h;
                    if (valid) {
                        double rate = swRates[sw];
                        if (geometricRates) {
                            rate *= getModelRateScalar(numBaseModel - h - 1) /// numBaseModel; // TODO Why not: "/ numBaseModel" ??
                                    / totalRate;
                        }
                        for (int i = 0; i < baseStateCount; ++i) {
                            matrix[g * baseStateCount + i][h * baseStateCount + i] = rate;
                        }
                        sw++;
                    }
                }
            }
        }

//        if (DEBUG) {
//            System.err.println(new Matrix(matrix));
//        }
    }

//    protected double setupMatrix() {
////        System.err.println("In MM.setupMatrix");
////        setupRelativeRates(relativeRates);
////        double[] pi = freqModel.getFrequencies();
//        setupQMatrix(null, null, q);
////        makeValid(q, stateCount);
//        return 1.0;
//    }

//    public FrequencyModel getFrequencyModel() {
//        return pcFreqModel;
//    }

    // TODO Remove
    public EigenDecomposition getEigenDecomposition() {
        if (DEBUG) {
            System.err.println("MMSM.getED");
        }
        EigenDecomposition ed = super.getEigenDecomposition();

        if (DEBUG) {
            double[][] q = getQCopy();
            System.err.println(new Matrix(q));
            System.err.println("");
            System.err.println(new dr.math.matrixAlgebra.Vector(ed.getEigenValues()));
            System.err.println("");
            double[] tp = new double[q.length * q.length];
            //the code below gets the matrix entries (after possible renormalization)
            getTransitionProbabilities(1.0, tp, ed);
            System.err.println(new Vector(tp));

            double[] infinitesimal = new double[q.length * q.length];
            getInfinitesimalMatrix(infinitesimal);
            double[][] matrixForm = new double[q.length][q.length];
            int i = 0,j = 0;
            for (int k = 0; k < infinitesimal.length; k++) {
                matrixForm[i][j] = infinitesimal[k];
                j++;
                if ((k+1) % q.length == 0) {
                    i++;
                    j = 0;
                }
            }
            System.err.println("Infinitesimal Q matrix:\n" + new Matrix(matrixForm));

            //System.exit(0);
        }

        return ed;
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.SUBSTITUTION_MODELS;
    }

    @Override
    public String getDescription() {
        return "Markov modulated substitution model";
    }

    public List<Citation> getCitations() {
        return Collections.singletonList(CommonCitations.SUCHARD_2012);
    }

    @Override
    protected void frequenciesChanged() {
        // Do nothing
    }

    @Override
    protected void ratesChanged() {
        updateMatrix = true;  // Lazy recompute relative rates
    }

    @Override
    protected void setupRelativeRates(double[] rates) {
        // Do nothing
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (DEBUG) {
            System.err.println("MMSM.hMCE");
        }
        // base substitution model changed!
        updateMatrix = true;
//        frequenciesChanged();
//        System.err.println("Model " + model.getId() + " changed");
        fireModelChanged();
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
//        System.err.println("Variable " + variable.getId() + " changed");
        if (variable == switchingRates || variable == rateScalar) {
            // Update rates
            updateMatrix = true;
//            if (variable == rateScalar)
                fireModelChanged(); // TODO Determine if necessary for ExposeRateCategoriesWrapper
        }
        // else do nothing, action taken care of at individual base models
    }

    public LogColumn[] getColumns() {

        List<LogColumn> columns = new ArrayList<LogColumn>();
        for (LogColumn parentColumn : super.getColumns()) {
            columns.add(parentColumn);
        }

        for (int i = 0; i < numBaseModel; ++i) {
            String label = "rateScalar." + i;
            columns.add(new RateColumn(label, i));
        }

        return columns.toArray(new LogColumn[0]);
    }

    private class RateColumn extends NumberColumn {

        private final int index;

        public RateColumn(String label, int index) {
            super(label);
            this.index = index;
        }

        /**
         * Returns the current value as a double.
         */
        @Override
        public double getDoubleValue() {
            return getModelRateScalar(index);
        }
    }
}
