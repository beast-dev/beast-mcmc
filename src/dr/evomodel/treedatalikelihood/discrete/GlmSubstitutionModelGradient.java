/*
 * DiscreteTraitBranchRateGradient.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.substmodel.OldGLMSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.ProcessSimulationDelegate;
import dr.inference.distribution.GeneralizedLinearModel;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * @author Marc A. Suchard
 */

@SuppressWarnings("deprecation")
public class GlmSubstitutionModelGradient implements GradientWrtParameterProvider, Reportable,
        Loggable, Citable {

    protected final TreeDataLikelihood treeDataLikelihood;
    protected final TreeTrait treeTraitProvider;
    protected final Tree tree;

    protected final OldGLMSubstitutionModel substitutionModel;
    protected final GeneralizedLinearModel glm;
    protected final int stateCount;

    private final ParameterMap parameterMap;
    private final int whichSubstitutionModel;

    public GlmSubstitutionModelGradient(String traitName,
                                        TreeDataLikelihood treeDataLikelihood,
                                        BeagleDataLikelihoodDelegate likelihoodDelegate,
                                        OldGLMSubstitutionModel substitutionModel) {

        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = treeDataLikelihood.getTree();
        this.substitutionModel = substitutionModel;
        this.glm = substitutionModel.getGeneralizedLinearModel();
        this.parameterMap = makeParameterMap(glm);
        this.stateCount = substitutionModel.getDataType().getStateCount();
        this.whichSubstitutionModel = determineSubstitutionNumber(
                likelihoodDelegate.getBranchModel(), substitutionModel);

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
    }

    private int determineSubstitutionNumber(BranchModel branchModel,
                                            OldGLMSubstitutionModel substitutionModel) {

        List<SubstitutionModel> substitutionModels = branchModel.getSubstitutionModels();
        for (int i = 0; i < substitutionModels.size(); ++i) {
            if (substitutionModel == substitutionModels.get(i)) {
                return i;
            }
        }
        throw new IllegalArgumentException("Unknown substitution model");
    }

    ParameterMap makeParameterMap(GeneralizedLinearModel glm) {

        final List<Integer> whichBlock = new ArrayList<>();
        final List<Integer> whichIndex = new ArrayList<>();

        CompoundParameter cp = new CompoundParameter("fixedEffects");
        boolean multi = glm.getNumberOfFixedEffects() > 1;

        for (int i = 0; i < glm.getNumberOfFixedEffects(); ++i) {
            Parameter p = glm.getFixedEffect(i);
            if (multi) {
                cp.addParameter(p);
            }
            for (int j = 0; j < p.getDimension(); ++j) {
                whichBlock.add(i);
                whichIndex.add(j);
            }

            if (glm.getFixedEffectIndicator(i) != null) {
                throw new IllegalArgumentException("GLM fixed effects gradients do not currently work with indicator variables");
            }
        }

        final Parameter whichParameter = multi ? cp : glm.getFixedEffect(0);
        return new ParameterMap() {

            @Override
            public double[] getCovariateColumn(int i) {
                return glm.getDesignMatrix(whichBlock.get(i)).getColumnValues(whichIndex.get(i));
            }

            @Override
            public Parameter getParameter() { return whichParameter; }
        };
    }

    @Override
    public Likelihood getLikelihood() {
        return treeDataLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return parameterMap.getParameter();
    }

    @Override
    public int getDimension() {
        return getParameter().getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {

        long startTime;
        if (COUNT_TOTAL_OPERATIONS) {
            startTime = System.nanoTime();
        }

        double[] differentials = (double[]) treeTraitProvider.getTrait(tree, null);
        double[] generator = new double[differentials.length];

        if (whichSubstitutionModel > 0) {
            final int length = stateCount * stateCount;
            System.arraycopy(
                    differentials, whichSubstitutionModel * length,
                    differentials, 0, length);
        }

        if (DEBUG_CROSS_PRODUCTS) {
            savedDifferentials = differentials.clone();
        }

        substitutionModel.getInfinitesimalMatrix(generator);
        double[] pi = substitutionModel.getFrequencyModel().getFrequencies();

        double normalizationConstant = preProcessNormalization(differentials, generator,
                substitutionModel.getNormalization());

        final double[] gradient = new double[getParameter().getDimension()];
        for (int i = 0; i < getParameter().getDimension(); ++i) {
            gradient[i] = processSingleGradientDimension(i, differentials, generator, pi,
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

    double preProcessNormalization(double[] differentials, double[] generator,
                                   boolean normalize) {
        return 0.0;
    }

    double processSingleGradientDimension(int i,
                                          double[] differentials, double[] generator, double[] pi,
                                          boolean normalize, double normalizationConstant) {

        double[] covariate = parameterMap.getCovariateColumn(i);
        return calculateCovariateDifferential(generator, differentials, covariate, pi, normalize);
    }

    private double calculateCovariateDifferential(double[] generator, double[] differential,
                                                  double[] covariate, double[] pi,
                                                  boolean doNormalization) {

        double normalization = 0.0;
        double total = 0.0;

        int k = 0;
        for (int i = 0; i < stateCount; ++i) {
            for (int j = i + 1; j < stateCount; ++j) {

                double xij = covariate[k++];
                double element = xij * generator[index(i,j)];

                total += differential[index(i,j)] * element;
                total -= differential[index(i,i)] * element;

                normalization += element * pi[i];
            }
        }

        for (int j = 0; j < stateCount; ++j) {
            for (int i = j + 1; i < stateCount; ++i) {

                double xij = covariate[k++];
                double element = xij * generator[index(i,j)];

                total += differential[index(i,j)] * element;
                total -= differential[index(i,i)] * element;

                normalization += element * pi[i];
            }
        }

        if (doNormalization) {
            for (int i = 0; i < stateCount; ++i) {
                for (int j = 0; j < stateCount; ++j) {
                    total -= differential[index(i,j)] * generator[index(i,j)] * normalization;
                }
            }
        }

        return total;
    }

    int index(int i, int j) {
        return i * stateCount + j;
    }

    @Override
    public String getReport() {

        StringBuilder sb = new StringBuilder();

        String message = GradientWrtParameterProvider.getReportAndCheckForError(this, 0.0, Double.POSITIVE_INFINITY, null);
        sb.append(message);

        if (DEBUG_CROSS_PRODUCTS) {
            sb.append("\n\tdifferentials: ").append(new WrappedVector.Raw(savedDifferentials, 0, savedDifferentials.length));
        }

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

    private static final boolean COUNT_TOTAL_OPERATIONS = true;
    private static final boolean DEBUG_CROSS_PRODUCTS = true;

    private double[] savedDifferentials;

    private long gradientCount = 0;
    private long totalGradientTime = 0;

    @Override
    public LogColumn[] getColumns() {
        return Loggable.getColumnsFromReport(this, "gradient report");
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.FRAMEWORK;
    }

    @Override
    public String getDescription() {
        return "Using linear-time differential calculations for all substitution generator elements";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    private static final Citation CITATION = new Citation(
            new Author[]{
                    new Author( "A", "Magee"),
                    new Author("P", "Lemey"),
                    new Author("MA", "Suchard"),
            },
            "Phylo-geographic GLM random effects",
            "",
            Citation.Status.IN_PREPARATION);

    interface ParameterMap {
        double[] getCovariateColumn(int i);
        Parameter getParameter();
    }
}
