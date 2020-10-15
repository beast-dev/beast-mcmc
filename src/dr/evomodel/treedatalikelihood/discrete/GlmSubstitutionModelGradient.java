/*
 * DiscreteTraitBranchRateGradient.java
 *
 * Copyright (c) 2002-2020 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evomodel.substmodel.OldGLMSubstitutionModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.ProcessSimulationDelegate;
import dr.inference.distribution.GeneralizedLinearModel;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.CompoundParameter;
import dr.inference.model.DesignMatrix;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.Reportable;

import java.util.Collections;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class GlmSubstitutionModelGradient implements GradientWrtParameterProvider, Reportable, Loggable, Citable {

    private final TreeDataLikelihood treeDataLikelihood;
    private final TreeTrait treeTraitProvider;
    private final Tree tree;

    private final OldGLMSubstitutionModel substitutionModel;
    private final GeneralizedLinearModel glm;
    private final Parameter allCoefficients;
    private final int stateCount;


    public GlmSubstitutionModelGradient(String traitName,
                                        TreeDataLikelihood treeDataLikelihood,
                                        BeagleDataLikelihoodDelegate likelihoodDelegate,
                                        OldGLMSubstitutionModel substitutionModel) {

        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = treeDataLikelihood.getTree();
        this.substitutionModel = substitutionModel;
        this.glm = substitutionModel.getGeneralizedLinearModel();
        this.allCoefficients = makeCompoundParameter(glm);
        this.stateCount = substitutionModel.getDataType().getStateCount();

        String name = SubstitutionModelCrossProductDelegate.getName(traitName);
        TreeTrait test = treeDataLikelihood.getTreeTrait(name);

        if (test == null) {
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

    public static Parameter makeCompoundParameter(GeneralizedLinearModel glm) {
        CompoundParameter parameter = new CompoundParameter("test");
        for (int i = 0; i < glm.getNumberOfFixedEffects(); ++i) {
            parameter.addParameter(glm.getFixedEffect(i));
        }
        return parameter;
    }

    @Override
    public Likelihood getLikelihood() {
        return treeDataLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return allCoefficients;
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

        substitutionModel.getInfinitesimalMatrix(generator);

        System.err.println("D: " + new WrappedVector.Raw(differentials));
        System.err.println("G: " + new WrappedVector.Raw(generator));

//        System.err.println("length = " + differentials.length);

        double[] pi = substitutionModel.getFrequencyModel().getFrequencies();

        double[] gradient = new double[1];

        DesignMatrix designMatrix = glm.getDesignMatrix(0);
        double[] covariate = designMatrix.getColumnValues(0);
        double[] covariateT = transposeOffDiagonal(covariate);

//        double[] chainDifferential = new double[stateCount * stateCount];
//
//        chainDifferential[0]  = -generator[2];
//        chainDifferential[2]  =  generator[2];
//        chainDifferential[5]  = -generator[7] - generator[6];
//        chainDifferential[6]  =  generator[6];
//        chainDifferential[7]  =  generator[7];
//        chainDifferential[8]  =  generator[8];
//        chainDifferential[10] = -generator[8];
//        chainDifferential[13] =  generator[13];
//        chainDifferential[15] = -generator[13];
//
//        double normalizationDifferential = calculateNormalizationDifferential(generator, covariate, pi);
//
//        for (int i = 0; i < stateCount * stateCount; ++i) {
//            chainDifferential[i] -= generator[i] * normalizationDifferential;
//        }

//        double[] chainDifferential = calculateCovariateDifferential(generator, differentials, covariate, pi,   false);

//        gradient[0] = dotProduct(differentials, chainDifferential);
        double x = calculateCovariateDifferential(generator, differentials, covariate, pi, substitutionModel.getNormalization());

        System.err.println("x1: " + x);
        System.err.println("x2: " + calculateCovariateDifferential(
                transposeAll(generator), transposeAll(differentials), transposeOffDiagonal(covariate), pi, false));
        System.err.println("x3: " + calculateCovariateDifferential(
                generator, transposeAll(differentials), transposeOffDiagonal(covariate), pi, false));
        System.err.println("x4: " + calculateCovariateDifferential(
                 transposeAll(generator), differentials, transposeOffDiagonal(covariate), pi, false));
        System.err.println("x5: " + calculateCovariateDifferential(
                 transposeAll(generator), transposeAll(differentials), covariate, pi, false));
        System.err.println("x6: " + calculateCovariateDifferential(
                 transposeAll(generator), differentials, covariate, pi, false));
        System.err.println("x7: " + calculateCovariateDifferential(
                 generator, transposeAll(differentials), covariate, pi, false));
        System.err.println("x8: " + calculateCovariateDifferential(
                 generator, differentials, transposeOffDiagonal(covariate), pi, false));

        gradient[0] = x;


//        double normalization = g

//        gradient[0] = grad;

        if (COUNT_TOTAL_OPERATIONS) {
            ++gradientCount;
            long endTime = System.nanoTime();
            totalGradientTime += (endTime - startTime) / 1000000;
        }

        return gradient;
    }

    private double[] transposeOffDiagonal(double[] x) {
        double[] result = new double[x.length];
        int half = x.length / 2;
        System.arraycopy(x, half, result, 0, half);
        System.arraycopy(x, 0, result, half, half);

        return result;
    }

    private double[] transposeAll(double[] x) {
        double[] result = new double[stateCount * stateCount];

        for (int i = 0; i < stateCount; ++i) {
            for (int j = 0; j < stateCount; ++j) {
                result[index(j,i)] = x[index(i,j)];
            }
        }

        return result;
    }


    private double calculateCovariateDifferential(double[] generator, double[] differential, double[] covariate, double[] pi,
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

    private int index(int i, int j) {
        return i * stateCount + j;
    }

//    private int idx(int i, int j) {
//        if (j > i) {
//            return
//        }
//    }

    private double calculateNormalizationDifferential(double[] generator, double[] covariate, double[] pi) {

        double total = 0.0;

        int k = 0;
        for (int i = 0; i < stateCount; ++i) {
            for (int j = i + 1; j < stateCount; ++j) {
                double xij = covariate[k++];
                total += xij * generator[i * stateCount + j] * pi[i];
            }
        }

        for (int j = 0; j < stateCount; ++j) {
            for (int i = j + 1; i < stateCount; ++i) {
                double xij = covariate[k++];
                total += xij * generator[i * stateCount + j] * pi[i];
            }
        }

        return total;
    }

    private static double dotProduct(double[] x, double[] y) {
        assert x.length == y.length;

        double total = 0.0;
        for (int i = 0; i < x.length; ++i) {
            total += x[i] * y[i];
        }

        return total;
    }


    @Override
    public String getReport() {

        StringBuilder sb = new StringBuilder();
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

        String message = GradientWrtParameterProvider.getReportAndCheckForError(this, 0.0, Double.POSITIVE_INFINITY, null);
        sb.append(message);

        return  sb.toString();
    }

    private static final boolean COUNT_TOTAL_OPERATIONS = true;
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
                    new Author("P", "Lemey"),
                    new Author("MA", "Suchard"),
            },
            "Phylogeographic GLM random effects",
            "",
            Citation.Status.IN_PREPARATION);
}
