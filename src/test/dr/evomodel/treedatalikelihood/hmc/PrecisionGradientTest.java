/*
 * PrecisionGradientTest.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package test.dr.evomodel.treedatalikelihood.hmc;

import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodel.treedatalikelihood.hmc.CorrelationPrecisionGradient;
import dr.evomodel.treedatalikelihood.hmc.DiagonalPrecisionGradient;
import dr.evomodel.treedatalikelihood.hmc.GradientWrtPrecisionProvider;
import dr.evomodel.treedatalikelihood.preorder.ConditionalOnTipsRealizedDelegate;
import dr.evomodel.treedatalikelihood.preorder.MultivariateConditionalOnTipsRealizedDelegate;
import dr.evomodel.treedatalikelihood.preorder.ProcessSimulationDelegate;
import dr.evomodel.treedatalikelihood.preorder.TipRealizedValuesViaFullConditionalDelegate;
import dr.inference.model.CachedMatrixInverse;
import dr.inference.model.CompoundParameter;
import dr.inference.model.CompoundSymmetricMatrix;
import dr.inference.model.Parameter;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * @author Paul Bastide
 * @author Marc Suchard
 */

public class PrecisionGradientTest extends TraceCorrelationAssert {

    private int dim;
    private CompoundSymmetricMatrix precisionMatrix;
    private CachedMatrixInverse precisionMatrixInv;

    private MultivariateDiffusionModel diffusionModel;
    private ContinuousTraitPartialsProvider dataModel;
    private ConjugateRootTraitPrior rootPrior;

    private NumberFormat format = NumberFormat.getNumberInstance(Locale.ENGLISH);

    public PrecisionGradientTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        format.setMaximumFractionDigits(2);

        dim = 6;

        Parameter offDiagonal = new Parameter.Default(new double[]{0.12, -0.13, 0.14, -0.15, 0.16,
                -0.12, 0.13, -0.14, 0.15,
                0.12, -0.13, 0.14,
                -0.12, 0.13,
                0.12});

        Parameter diagonal = new Parameter.Default(new double[]{0.1, 0.2, 0.3, 0.4, 0.5, 0.6});

        precisionMatrix = new CompoundSymmetricMatrix(diagonal, offDiagonal, true, true);

        Parameter diagonalVar = new Parameter.Default(new double[]{1.0, 20.0, 30.0, 40.0, 50.0, 60.0});

        precisionMatrixInv = new CachedMatrixInverse("var",
                new CompoundSymmetricMatrix(diagonalVar, offDiagonal, true, false));

        // Data
        Parameter[] dataTraits = new Parameter[6];
        dataTraits[0] = new Parameter.Default("human", new double[]{-1.0, 2.0, 3.0, 4.0, 5.0, -6.0});
        dataTraits[1] = new Parameter.Default("chimp", new double[]{10.0, 12.0, 14.0, 16.0, 18.0, 20.0});
        dataTraits[2] = new Parameter.Default("bonobo", new double[]{0.5, -2.0, 5.5, -5.2, 3.1, 1.1});
        dataTraits[3] = new Parameter.Default("gorilla", new double[]{2.0, 5.0, -8.0, -4.0, 3.2, 3.4});
        dataTraits[4] = new Parameter.Default("orangutan", new double[]{11.0, 1.0, -1.5, 2.4, -4.2, 6.0});
        dataTraits[5] = new Parameter.Default("siamang", new double[]{1.0, 2.5, 4.0, 4.0, -5.2, 1.0});
        CompoundParameter traitParameter = new CompoundParameter("trait", dataTraits);

        List<Integer> missingIndices = new ArrayList<Integer>();
        traitParameter.setParameterValue(2, 0);

        // Tree
        createAlignment(PRIMATES_TAXON_SEQUENCE, Nucleotides.INSTANCE);
        treeModel = createPrimateTreeModel();

        //// Standard Model //// ***************************************************************************************

        PrecisionType precisionType = PrecisionType.FULL;

        // Root prior
//        rootPrior = new ConjugateRootTraitPrior(new double[]{-1.0, -3.0, 2.5, -2.5, 1.3, 4.0}, 10.0, true);

        // Data Model
        dataModel = new ContinuousTraitDataModel("dataModel",
                traitParameter,
                missingIndices, true,
                6, precisionType);
    }

    public void testGradientPrecision() {
        System.out.println("\nTest gradient precision.");

        // Diffusion
        diffusionModel = new MultivariateDiffusionModel(precisionMatrix);
        DiffusionProcessDelegate diffusionProcessDelegate
                = new HomogeneousDiffusionModelDelegate(treeModel, diffusionModel);

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel, true);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihood = new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel);

        ProcessSimulationDelegate simulationDelegate =
                likelihoodDelegate.getPrecisionType() == PrecisionType.SCALAR ?
                        new ConditionalOnTipsRealizedDelegate("trait", treeModel,
                                diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate) :
                        new MultivariateConditionalOnTipsRealizedDelegate("trait", treeModel,
                                diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);

        TreeTraitProvider traitProvider = new ProcessSimulation(dataLikelihood, simulationDelegate);

        dataLikelihood.addTraits(traitProvider.getTreeTraits());

        ProcessSimulationDelegate fullConditionalDelegate = new TipRealizedValuesViaFullConditionalDelegate(
                "trait", treeModel, diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);

        dataLikelihood.addTraits(new ProcessSimulation(dataLikelihood, fullConditionalDelegate).getTreeTraits());


        // Wishart Statistic
        ContinuousDataLikelihoodDelegate cdld = (ContinuousDataLikelihoodDelegate) dataLikelihood.getDataLikelihoodDelegate();
        WishartStatisticsWrapper wishartStatistics
                = new WishartStatisticsWrapper("wishart", "trait", dataLikelihood,
                cdld);

        GradientWrtPrecisionProvider gPPWiwhart = new GradientWrtPrecisionProvider.WishartGradientWrtPrecisionProvider(wishartStatistics);

        // Branch Specific
        ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient traitGradient =
                new ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient(
                        dim, treeModel, cdld,
                        ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter.WRT_PRECISION);
        BranchSpecificGradient branchSpecificGradient =
                new BranchSpecificGradient("trait", dataLikelihood, cdld, traitGradient, precisionMatrix);

        GradientWrtPrecisionProvider gPPBranchSpecific = new GradientWrtPrecisionProvider.BranchSpecificGradientWrtPrecisionProvider(branchSpecificGradient);


        // Correlation Gradient
        CorrelationPrecisionGradient gradientProviderWishart = new CorrelationPrecisionGradient(gPPWiwhart, dataLikelihood, precisionMatrix);

        String sW = gradientProviderWishart.getReport();
        System.err.println(sW);
        double[] gradientAnalyticalW = parseGradient(sW, "analytic");
        double[] gradientNumeric = parseGradient(sW, "numeric (with Cholesky):");

        assertEquals("Sizes", gradientAnalyticalW.length, gradientNumeric.length);

        for (int k = 0; k < gradientAnalyticalW.length; k++) {
            assertEquals("gradient correlation k=" + k,
                    format.format(gradientAnalyticalW[k]),
                    format.format(gradientNumeric[k]));
        }

        // Correlation Gradient Branch Specific
        CorrelationPrecisionGradient gradientProviderBranchSpecific = new CorrelationPrecisionGradient(gPPBranchSpecific, dataLikelihood, precisionMatrix);

        String sBS = gradientProviderBranchSpecific.getReport();
        System.err.println(sBS);
        double[] gradientAnalyticalBS = parseGradient(sBS, "analytic");

        assertEquals("Sizes", gradientAnalyticalBS.length, gradientNumeric.length);

        for (int k = 0; k < gradientAnalyticalBS.length; k++) {
            assertEquals("gradient correlation k=" + k,
                    format.format(gradientAnalyticalBS[k]),
                    format.format(gradientAnalyticalW[k]));
        }

        // Diagonal Gradient
        DiagonalPrecisionGradient gradientDiagonalProviderW = new DiagonalPrecisionGradient(gPPWiwhart, dataLikelihood, precisionMatrix);

        String sDiagW = gradientDiagonalProviderW.getReport();
        System.err.println(sDiagW);
        double[] gradientDiagonalAnalyticalW = parseGradient(sDiagW, "analytic");
        double[] gradientDiagonalNumeric = parseGradient(sDiagW, "numeric:");

        assertEquals("Sizes", gradientDiagonalAnalyticalW.length, gradientDiagonalNumeric.length);

        for (int k = 0; k < gradientDiagonalAnalyticalW.length; k++) {
            assertEquals("gradient diagonal k=" + k,
                    format.format(gradientDiagonalAnalyticalW[k]),
                    format.format(gradientDiagonalNumeric[k]));
        }

        // Diagonal Gradient Branch Specific
        DiagonalPrecisionGradient gradientDiagonalProviderBS = new DiagonalPrecisionGradient(gPPBranchSpecific, dataLikelihood, precisionMatrix);

        String sDiagBS = gradientDiagonalProviderBS.getReport();
        System.err.println(sDiagBS);
        double[] gradientDiagonalAnalyticalBS = parseGradient(sDiagBS, "analytic");

        assertEquals("Sizes", gradientDiagonalAnalyticalBS.length, gradientDiagonalAnalyticalW.length);

        for (int k = 0; k < gradientDiagonalAnalyticalBS.length; k++) {
            assertEquals("gradient correlation k=" + k,
                    format.format(gradientDiagonalAnalyticalBS[k]),
                    format.format(gradientDiagonalAnalyticalW[k]));
        }
    }

    public void testGradientVariance() {
        System.out.println("\nTest gradient variance.");

        // Diffusion
        diffusionModel = new MultivariateDiffusionModel(precisionMatrixInv);
        DiffusionProcessDelegate diffusionProcessDelegate
                = new HomogeneousDiffusionModelDelegate(treeModel, diffusionModel);

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel, true);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihood = new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel);

        ProcessSimulationDelegate simulationDelegate =
                likelihoodDelegate.getPrecisionType() == PrecisionType.SCALAR ?
                        new ConditionalOnTipsRealizedDelegate("trait", treeModel,
                                diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate) :
                        new MultivariateConditionalOnTipsRealizedDelegate("trait", treeModel,
                                diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);

        TreeTraitProvider traitProvider = new ProcessSimulation(dataLikelihood, simulationDelegate);

        dataLikelihood.addTraits(traitProvider.getTreeTraits());

        ProcessSimulationDelegate fullConditionalDelegate = new TipRealizedValuesViaFullConditionalDelegate(
                "trait", treeModel, diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);

        dataLikelihood.addTraits(new ProcessSimulation(dataLikelihood, fullConditionalDelegate).getTreeTraits());


        // Wishart Statistic
        ContinuousDataLikelihoodDelegate cdld = (ContinuousDataLikelihoodDelegate) dataLikelihood.getDataLikelihoodDelegate();
        WishartStatisticsWrapper wishartStatistics
                = new WishartStatisticsWrapper("wishart", "trait", dataLikelihood, cdld);

        GradientWrtPrecisionProvider gradientWrtPrecisionProvider = new GradientWrtPrecisionProvider.WishartGradientWrtPrecisionProvider(wishartStatistics);

        // Branch Specific
        ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient traitGradient =
                new ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient(
                        dim, treeModel, cdld,
                        ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter.WRT_PRECISION);
        BranchSpecificGradient branchSpecificGradient =
                new BranchSpecificGradient("trait", dataLikelihood, cdld, traitGradient, precisionMatrix);

        GradientWrtPrecisionProvider gPPBranchSpecific = new GradientWrtPrecisionProvider.BranchSpecificGradientWrtPrecisionProvider(branchSpecificGradient);

        // Correlation Gradient
        CorrelationPrecisionGradient gradientProviderW = new CorrelationPrecisionGradient(gradientWrtPrecisionProvider, dataLikelihood, precisionMatrixInv);

        String sW = gradientProviderW.getReport();
        System.err.println(sW);
        double[] gradientAnalyticalW = parseGradient(sW, "analytic");
        double[] gradientNumeric = parseGradient(sW, "numeric (with Cholesky):");

        assertEquals("Sizes", gradientAnalyticalW.length, gradientNumeric.length);

        for (int k = 0; k < gradientAnalyticalW.length; k++) {
            assertEquals("gradient correlation k=" + k,
                    format.format(gradientAnalyticalW[k]),
                    format.format(gradientNumeric[k]));
        }

        // Correlation Gradient Branch Specific
        CorrelationPrecisionGradient gradientProviderBranchSpecific = new CorrelationPrecisionGradient(gPPBranchSpecific, dataLikelihood, precisionMatrixInv);

        String sBS = gradientProviderBranchSpecific.getReport();
        System.err.println(sBS);
        double[] gradientAnalyticalBS = parseGradient(sBS, "analytic");

        assertEquals("Sizes", gradientAnalyticalBS.length, gradientNumeric.length);

        for (int k = 0; k < gradientAnalyticalBS.length; k++) {
            assertEquals("gradient correlation k=" + k,
                    format.format(gradientAnalyticalBS[k]),
                    format.format(gradientAnalyticalW[k]));
        }

        // Diagonal Gradient
        DiagonalPrecisionGradient gradientDiagonalProviderW = new DiagonalPrecisionGradient(gradientWrtPrecisionProvider, dataLikelihood, precisionMatrixInv);

        String sDiagW = gradientDiagonalProviderW.getReport();
        System.err.println(sDiagW);
        double[] gradientDiagonalAnalyticalW = parseGradient(sDiagW, "analytic");
        double[] gradientDiagonalNumeric = parseGradient(sDiagW, "numeric:");

        assertEquals("Sizes", gradientDiagonalAnalyticalW.length, gradientDiagonalNumeric.length);

        for (int k = 0; k < gradientDiagonalAnalyticalW.length; k++) {
            assertEquals("gradient diagonal k=" + k,
                    format.format(gradientDiagonalAnalyticalW[k]),
                    format.format(gradientDiagonalNumeric[k]));
        }

        // Diagonal Gradient Branch Specific
        DiagonalPrecisionGradient gradientDiagonalProviderBS = new DiagonalPrecisionGradient(gPPBranchSpecific, dataLikelihood, precisionMatrixInv);

        String sDiagBS = gradientDiagonalProviderBS.getReport();
        System.err.println(sDiagBS);
        double[] gradientDiagonalAnalyticalBS = parseGradient(sDiagBS, "analytic");

        assertEquals("Sizes", gradientDiagonalAnalyticalBS.length, gradientDiagonalAnalyticalW.length);

        for (int k = 0; k < gradientDiagonalAnalyticalBS.length; k++) {
            assertEquals("gradient correlation k=" + k,
                    format.format(gradientDiagonalAnalyticalBS[k]),
                    format.format(gradientDiagonalAnalyticalW[k]));
        }
    }

    private double[] parseGradient(String s, String name) {
        int indBeg = s.indexOf(name) + name.length() + 3;
        int indEnd = s.indexOf("]", indBeg);
        return parseVector(s.substring(indBeg, indEnd), ",");
    }

    private double[] parseVector(String s, String sep) {
        String[] vectorString = s.split(sep);
        double[] gradient = new double[vectorString.length];
        for (int i = 0; i < gradient.length; i++) {
            gradient[i] = Double.parseDouble(vectorString[i]);
        }
        return gradient;
    }

    public void testGradientPrecisionWithDrift() {
        System.out.println("\nTest gradient precision.");

        // Diffusion
        diffusionModel = new MultivariateDiffusionModel(precisionMatrix);
        List<BranchRateModel> driftModels = new ArrayList<BranchRateModel>();
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{0.0})));
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{200.0})));
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.3", new double[]{-200.0})));
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.4", new double[]{1.0})));
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.5", new double[]{200.0})));
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.6", new double[]{-200.0})));
        DiffusionProcessDelegate diffusionProcessDelegate
                = new DriftDiffusionModelDelegate(treeModel, diffusionModel, driftModels);

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel, true);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihood = new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel);

        ProcessSimulationDelegate simulationDelegate =
                likelihoodDelegate.getPrecisionType() == PrecisionType.SCALAR ?
                        new ConditionalOnTipsRealizedDelegate("trait", treeModel,
                                diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate) :
                        new MultivariateConditionalOnTipsRealizedDelegate("trait", treeModel,
                                diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);

        TreeTraitProvider traitProvider = new ProcessSimulation(dataLikelihood, simulationDelegate);

        dataLikelihood.addTraits(traitProvider.getTreeTraits());

        ProcessSimulationDelegate fullConditionalDelegate = new TipRealizedValuesViaFullConditionalDelegate(
                "trait", treeModel, diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);

        dataLikelihood.addTraits(new ProcessSimulation(dataLikelihood, fullConditionalDelegate).getTreeTraits());

        // Branch Specific
        ContinuousDataLikelihoodDelegate cdld = (ContinuousDataLikelihoodDelegate) dataLikelihood.getDataLikelihoodDelegate();

        ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient traitGradient =
                new ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient(
                        dim, treeModel, cdld,
                        ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter.WRT_PRECISION);
        BranchSpecificGradient branchSpecificGradient =
                new BranchSpecificGradient("trait", dataLikelihood, cdld, traitGradient, precisionMatrix);

        GradientWrtPrecisionProvider gPPBranchSpecific = new GradientWrtPrecisionProvider.BranchSpecificGradientWrtPrecisionProvider(branchSpecificGradient);

        // Correlation Gradient Branch Specific
        CorrelationPrecisionGradient gradientProviderBranchSpecific = new CorrelationPrecisionGradient(gPPBranchSpecific, dataLikelihood, precisionMatrix);

        String sBS = gradientProviderBranchSpecific.getReport();
        System.err.println(sBS);
        double[] gradientAnalyticalBS = parseGradient(sBS, "analytic");
        double[] gradientNumeric = parseGradient(sBS, "numeric (with Cholesky):");

        assertEquals("Sizes", gradientAnalyticalBS.length, gradientNumeric.length);

        for (int k = 0; k < gradientAnalyticalBS.length; k++) {
            assertEquals("gradient correlation k=" + k,
                    format.format(gradientAnalyticalBS[k]),
                    format.format(gradientNumeric[k]));
        }

        // Diagonal Gradient Branch Specific
        DiagonalPrecisionGradient gradientDiagonalProviderBS = new DiagonalPrecisionGradient(gPPBranchSpecific, dataLikelihood, precisionMatrix);

        String sDiagBS = gradientDiagonalProviderBS.getReport();
        System.err.println(sDiagBS);
        double[] gradientDiagonalAnalyticalBS = parseGradient(sDiagBS, "analytic");
        double[] gradientDiagonalNumeric = parseGradient(sDiagBS, "numeric:");

        assertEquals("Sizes", gradientDiagonalAnalyticalBS.length, gradientDiagonalNumeric.length);

        for (int k = 0; k < gradientDiagonalAnalyticalBS.length; k++) {
            assertEquals("gradient correlation k=" + k,
                    format.format(gradientDiagonalAnalyticalBS[k]),
                    format.format(gradientDiagonalNumeric[k]));
        }
    }
}
