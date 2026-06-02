/*
 * AbstractLogAdditiveSubstitutionModelGradient.java
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

package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.substmodel.*;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.ProcessSimulationDelegate;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.loggers.Loggable;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.ModelListener;
import dr.util.Citable;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.List;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;

/**
 * @author Marc A. Suchard
 */

public abstract class AbstractLogAdditiveSubstitutionModelGradient implements
        GradientWrtParameterProvider, ModelListener, Reportable, Loggable, Citable {

    protected final TreeDataLikelihood treeDataLikelihood;
    protected final TreeTrait treeTraitProvider;
    protected final Tree tree;
    protected final BranchModel branchModel;

    protected final ComplexSubstitutionModel substitutionModel;
    protected final int stateCount;
    protected final List<Integer> crossProductAccumulationMap;

    private final ApproximationMode mode;

    public enum ApproximationMode {

        FIRST_ORDER("firstOrder") {
            @Override
            public String getInfo() {
                return "a first-order";
            }

            @Override
            CorrectionTermCache createCache(SubstitutionModel model,
                                            List<Integer> accumulateMap) { return null; }

            @Override
            void emptyCache(CorrectionTermCache cache) { }

            @Override
            double computeCorrection(int i, int j, double[] crossProducts, int stateCount,
                                     CorrectionTermCache correctionTermCache) {
                return 0.0;
            }
        },
        AFFINE_CORRECTED("affineCorrected") {
            @Override
            public String getInfo() {
                return "an affine-corrected";
            }

            @Override
            CorrectionTermCache createCache(SubstitutionModel model,
                                            List<Integer> accumulateMap) {
                if (accumulateMap.size() > 1) {
                    throw new RuntimeException("Not yet implemented");
                }
                return new CorrectionTermCache(model);
            }

            @Override
            void emptyCache(CorrectionTermCache cache) {
                cache.clear();
            }

            @Override
            double computeCorrection(int i, int j, double[] crossProducts, int stateCount,
                                     CorrectionTermCache correctionTermCache) {

                double[] affineMatrix = correctionTermCache.getAffineMatrix(i, j);

                double correction = 0.0;
                for (int m = 0; m < stateCount; ++m) {
                    for (int n = 0; n < stateCount; ++n) {
                        correction += crossProducts[m * stateCount + n] *
                                affineMatrix[m * stateCount + n];
                    }
                }

                return correction;
            }

        };

        ApproximationMode(String label) {
            this.label = label;
        }

        abstract String getInfo();

        abstract CorrectionTermCache createCache(SubstitutionModel model, List<Integer> accumulateMap);

        abstract void emptyCache(CorrectionTermCache cache);

        abstract double computeCorrection(int i, int j, double[] crossProducts, int stateCount,
                                          CorrectionTermCache correctionTermCache);

        public String getLabel() {
            return label;
        }

        private final String label;

        public static ApproximationMode factory(String label) {
            for (ApproximationMode mode : ApproximationMode.values()) {
                if (mode.getLabel().equalsIgnoreCase(label)) {
                    return mode;
                }
            }
            throw new IllegalArgumentException("Unknown approximation mode");
        }
    }

    private static final ApproximationMode DEFAULT_MODE = ApproximationMode.FIRST_ORDER;
    
    public AbstractLogAdditiveSubstitutionModelGradient(String traitName,
                                                        TreeDataLikelihood treeDataLikelihood,
                                                        BeagleDataLikelihoodDelegate likelihoodDelegate,
                                                        ComplexSubstitutionModel substitutionModel,
                                                        ApproximationMode mode) {
        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = treeDataLikelihood.getTree();
        this.branchModel = likelihoodDelegate.getBranchModel();
        this.substitutionModel = substitutionModel;
        this.stateCount = substitutionModel.getDataType().getStateCount();
        this.crossProductAccumulationMap = createCrossProductAccumulationMap(likelihoodDelegate.getBranchModel(),
                substitutionModel);
        this.mode = mode;
        this.correctionTermCache = mode.createCache(substitutionModel, crossProductAccumulationMap);

        String name = SubstitutionModelCrossProductDelegate.getName(traitName);

        if (treeDataLikelihood.getTreeTrait(name) == null) {
            ProcessSimulationDelegate gradientDelegate = new SubstitutionModelCrossProductDelegate(traitName,
                    treeDataLikelihood.getTree(),
                    likelihoodDelegate,
                    treeDataLikelihood.getBranchRateModel(),
                    substitutionModel.getDataType().getStateCount());
            TreeTraitProvider traitProvider = new ProcessSimulation(treeDataLikelihood, gradientDelegate);
            treeDataLikelihood.addTraits(traitProvider.getTreeTraits());
        }

        treeTraitProvider = treeDataLikelihood.getTreeTrait(name);
        assert (treeTraitProvider != null);

        this.branchModel.addModelListener(this);
        this.substitutionModel.addModelListener(this);
        
        Logger.getLogger("dr.evomodel.treedatalikelihood.discrete").info(
                "Gradient wrt " + traitName + " using " + mode.getInfo() + " approximation");
    }

    protected abstract double preProcessNormalization(double[] differentials, double[] generator, boolean normalize);

    abstract double processSingleGradientDimension(int dim,
                                                   double[] differentials, double[] generator, double[] pi,
                                                   boolean normalize, double normalizationConstant);

    @Override
    public double[] getGradientLogDensity() {

        long startTime;
        if (COUNT_TOTAL_OPERATIONS) {
            startTime = System.nanoTime();
        }

        double[] crossProducts = (double[]) treeTraitProvider.getTrait(tree, null);
        double[] generator = new double[crossProducts.length];

        accumulateAcrossSubstitutionModelInstances(crossProducts);

        substitutionModel.getInfinitesimalMatrix(generator);

        double[] pi = substitutionModel.getFrequencyModel().getFrequencies();

        double normalizationConstant = preProcessNormalization(crossProducts, generator,
                substitutionModel.getNormalization());

        final double[] gradient = new double[getParameter().getDimension()];
        for (int i = 0; i < getParameter().getDimension(); ++i) {
            gradient[i] = processSingleGradientDimension(i, crossProducts, generator, pi,
                    substitutionModel.getNormalization(),
                    normalizationConstant);
        }

        if (COUNT_TOTAL_OPERATIONS) {
            ++gradientCount;
            long endTime = System.nanoTime();
            totalGradientTime += (endTime - startTime) / 1000000;
        }

        return gradient;
    }

    double correction(int i, int j, double[] crossProducts) {
        return mode.computeCorrection(i, j, crossProducts, stateCount, correctionTermCache);
    }

    private void accumulateAcrossSubstitutionModelInstances(double[] crossProducts) {
        final int length = stateCount * stateCount;

        int firstModel = crossProductAccumulationMap.get(0);
        if (firstModel > 0) {
            // Copy first set of entries
            System.arraycopy(
                    crossProducts, firstModel * length,
                    crossProducts, 0, length);
        }

        for (int i = 1; i < crossProductAccumulationMap.size(); ++i) {
            int nextModel = crossProductAccumulationMap.get(i);
            for (int j = 0; j < length; ++j) {
                crossProducts[j] += crossProducts[nextModel * length + j];
            }
        }
    }

    private List<Integer> createCrossProductAccumulationMap(BranchModel branchModel,
                                                            ComplexSubstitutionModel substitutionModel) {

        List<SubstitutionModel> substitutionModels = branchModel.getSubstitutionModels();
        List<Integer> map = new ArrayList<>();
        
        for (int i = 0; i < substitutionModels.size(); ++i) {
            if (substitutionModel == substitutionModels.get(i)) {
                map.add(i);
            }
        }

        return map;
    }

    @Override
    public Likelihood getLikelihood() {
        return treeDataLikelihood;
    }

    @Override
    public int getDimension() {
        return getParameter().getDimension();
    }

    @Override
    public String getReport() {

        StringBuilder sb = new StringBuilder();

        String message = GradientWrtParameterProvider.getReportAndCheckForError(this, 0.0, Double.POSITIVE_INFINITY, getReportTolerance());
        sb.append(message);

        if (COUNT_TOTAL_OPERATIONS) {
            sb.append("\n\tgetCrossProductGradientCount = ").append(gradientCount);
            sb.append("\n\taverageGradientTime = ");
            if (gradientCount > 0) {
                sb.append(totalGradientTime / gradientCount);
            } else {
                sb.append("NA");
            }
            sb.append("\n");
        }

        return  sb.toString();
    }


    Double getReportTolerance() {
        return null;
    }

    @SuppressWarnings("unused")
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == branchModel) {
            throw new RuntimeException("Not yet implemented");
        } else if (model == substitutionModel) {
            mode.emptyCache(correctionTermCache);
        } else {
            throw new RuntimeException("Unknown model");
        }
    }

    public void modelChangedEvent(Model model, Object object, int index) {

    }

    public void modelRestored(Model model) {

    }

    static class CorrectionTermCache {

        private final SubstitutionModel model;
        private final Map<Integer, double[]> map;
        private final int stateCount;
        private double[] qQPlus;

        CorrectionTermCache(SubstitutionModel model) {
            this.model = model;
            this.map = new HashMap<>();
            this.stateCount = model.getDataType().getStateCount();
            this.qQPlus = null;
        }

        private int index12(int i, int j) {
            return i * stateCount + j;
        }

        private double[] getQQPlus() {
            if (qQPlus == null) {
                EigenDecomposition ed = model.getEigenDecomposition();
                qQPlus = DifferentiableSubstitutionModelUtil.getQQPlus(ed.getEigenVectors(),
                        ed.getInverseEigenVectors(), ed.getEigenValues(), stateCount);
            }
            return qQPlus;
        }

        double[] getAffineMatrix(int i, int j) {

            double[] affineMatrix = map.get(i * stateCount + j);

            if (affineMatrix == null) {

                affineMatrix = new double[stateCount * stateCount]; // TODO there are only stateCount unique values

                double[] qQPlus = getQQPlus();

                for (int m = 0; m < stateCount; ++m) {
                    for (int n = 0; n < stateCount; n++) {
                        affineMatrix[index12(m, n)] = (m == i) ?
                                (qQPlus[index12(m, i)] - 1.0) * qQPlus[index12(j, n)] :
                                qQPlus[index12(m, i)] * qQPlus[index12(j, n)];
                    }
                }

                map.put(i * stateCount + j, affineMatrix);
            }

            return affineMatrix;
        }

        public void clear() {
            map.clear();
            qQPlus = null;
        }
    }

    private final CorrectionTermCache correctionTermCache;

    protected static final boolean COUNT_TOTAL_OPERATIONS = false;
    protected long gradientCount = 0;
    protected long totalGradientTime = 0;
}
