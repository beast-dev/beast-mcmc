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
import dr.evomodel.continuous.MultivariateElasticModel;
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
import dr.inference.model.*;
import dr.xml.AttributeParser;
import dr.xml.XMLParser;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.io.StringReader;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Paul Bastide
 * @author Marc Suchard
 */

public class PrecisionGradientTest extends TraceCorrelationAssert {

    private int dim;
    private CompoundSymmetricMatrix precisionMatrix;
    private CachedMatrixInverse precisionMatrixInv;

    private MultivariateDiffusionModel diffusionModel;
    private MultivariateDiffusionModel diffusionModelVar;
    private ContinuousTraitPartialsProvider dataModel;
    private ContinuousTraitPartialsProvider dataModelMissing;
    private ConjugateRootTraitPrior rootPrior;

    private ContinuousRateTransformation rateTransformation;
    private BranchRateModel rateModel;

    private Boolean fixedRoot = false;

    private double delta;

    public PrecisionGradientTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        delta = 1E-3;

        dim = 6;

        Parameter offDiagonal = new Parameter.Default(new double[]{0.12, -0.13, 0.14, -0.15, 0.16,
                -0.12, 0.13, -0.14, 0.15,
                0.12, -0.13, 0.14,
                -0.12, 0.13,
                0.12});

        Parameter diagonal = new Parameter.Default(new double[]{0.1, 0.2, 0.3, 0.4, 0.5, 0.6});

        precisionMatrix = new CompoundSymmetricMatrix(diagonal, offDiagonal, true, true);

        Parameter diagonalVar = new Parameter.Default(new double[]{10.0, 20.0, 30.0, 40.0, 50.0, 60.0});

        precisionMatrixInv = new CachedMatrixInverse("var",
                new CompoundSymmetricMatrix(diagonalVar, offDiagonal, true, false));

        diffusionModel = new MultivariateDiffusionModel(precisionMatrix);
        diffusionModelVar = new MultivariateDiffusionModel(precisionMatrixInv);

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
        missingIndices.add(6);
        missingIndices.add(7);
        missingIndices.add(8);
        missingIndices.add(9);
        missingIndices.add(10);
        missingIndices.add(11);
        missingIndices.add(13);
        missingIndices.add(15);
        missingIndices.add(25);
        missingIndices.add(29);


        // Tree
        createAlignment(PRIMATES_TAXON_SEQUENCE, Nucleotides.INSTANCE);
        treeModel = createPrimateTreeModel();

        // Rates
        rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        rateModel = new DefaultBranchRateModel();

        //// Standard Model //// ***************************************************************************************

        PrecisionType precisionType = PrecisionType.FULL;

        // Root prior
        final String rootVal = fixedRoot ? "Infinity" : "0.1";
        String s = "<beast>\n" +
                "    <conjugateRootPrior>\n" +
                "        <meanParameter>\n" +
                "            <parameter id=\"meanRoot\"  value=\"-1.0 -3.0 2.5 -2.5 1.3 4.0\"/>\n" +
                "        </meanParameter>\n" +
                "        <priorSampleSize>\n" +
                "            <parameter id=\"sampleSizeRoot\" value=\"" + rootVal + "\"/>\n" +
                "        </priorSampleSize>\n" +
                "    </conjugateRootPrior>\n" +
                "</beast>";
        XMLParser parser = new XMLParser(true, true, true, null);
        parser.addXMLObjectParser(new AttributeParser());
        parser.addXMLObjectParser(new ParameterParser());
        parser.parse(new StringReader(s), true);
        rootPrior = ConjugateRootTraitPrior.parseConjugateRootTraitPrior(parser.getRoot(), dim);
//        rootPrior = new ConjugateRootTraitPrior(new double[]{-1.0, -3.0, 2.5, -2.5, 1.3, 4.0}, 10.0, true);

        // Data Model
        dataModelMissing = new ContinuousTraitDataModel("dataModel",
                traitParameter,
                missingIndices, true,
                6, precisionType);

        dataModel = new ContinuousTraitDataModel("dataModel",
                traitParameter,
                missingIndices, false,
                6, precisionType);

    }

    public void testGradientBMWithMissing() {

        // Wrt Precision
        DiffusionProcessDelegate diffusionProcessDelegate
                = new HomogeneousDiffusionModelDelegate(treeModel, diffusionModel);
        System.out.println("\nTest gradient precision.");
        testGradient(diffusionModel, diffusionProcessDelegate, dataModel, precisionMatrix, true);
        System.out.println("\nTest gradient precision with missing.");
        testGradient(diffusionModel, diffusionProcessDelegate, dataModelMissing, precisionMatrix, false);

        // Wrt Variance
        DiffusionProcessDelegate diffusionProcessDelegateVariance
                = new HomogeneousDiffusionModelDelegate(treeModel, diffusionModelVar);
        System.out.println("\nTest gradient variance.");
        testGradient(diffusionModelVar, diffusionProcessDelegateVariance, dataModel, precisionMatrixInv, true);
        System.out.println("\nTest gradient variance with missing.");
        testGradient(diffusionModelVar, diffusionProcessDelegateVariance, dataModelMissing, precisionMatrixInv, false);
    }

    public void testGradientDriftWithMissing() {
        System.out.println("\nTest gradient precision.");

        // Diffusion
        List<BranchRateModel> driftModels = new ArrayList<BranchRateModel>();
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{0.0})));
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{200.0})));
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.3", new double[]{-200.0})));
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.4", new double[]{1.0})));
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.5", new double[]{200.0})));
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.6", new double[]{-200.0})));

        // Wrt Precision
        DiffusionProcessDelegate diffusionProcessDelegate
                = new DriftDiffusionModelDelegate(treeModel, diffusionModel, driftModels);
        System.out.println("\nTest drift gradient precision.");
        testGradient(diffusionModel, diffusionProcessDelegate, dataModel, precisionMatrix, false);
        System.out.println("\nTest drift gradient precision with missing.");
        testGradient(diffusionModel, diffusionProcessDelegate, dataModelMissing, precisionMatrix, false);

        // Wrt Variance
        DiffusionProcessDelegate diffusionProcessDelegateVariance
                = new DriftDiffusionModelDelegate(treeModel, diffusionModelVar, driftModels);
        System.out.println("\nTest drift gradient variance.");
        testGradient(diffusionModelVar, diffusionProcessDelegateVariance, dataModel, precisionMatrixInv, false);
        System.out.println("\nTest drift gradient variance with missing.");
        testGradient(diffusionModelVar, diffusionProcessDelegateVariance, dataModelMissing, precisionMatrixInv, false);
    }

    public void testGradientOUWithMissing() {
        System.out.println("\nTest gradient variance with missing.");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{1.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{2.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.3", new double[]{-2.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.4", new double[]{1.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.5", new double[]{2.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.6", new double[]{-2.0})));

        Parameter[] strengthOfSelectionParameters = new Parameter[6];
        strengthOfSelectionParameters[0] = new Parameter.Default(new double[]{10, 1.0, 0.0, 0.0, 0.0, 2.0});
        strengthOfSelectionParameters[1] = new Parameter.Default(new double[]{1.0, 20, 0.0, 0.0, 0.0, 0.0});
        strengthOfSelectionParameters[2] = new Parameter.Default(new double[]{0.0, 0.0, 30, 0.0, 0.0, 0.0});
        strengthOfSelectionParameters[3] = new Parameter.Default(new double[]{0.0, 0.0, 0.0, 40, 3.0, 0.0});
        strengthOfSelectionParameters[4] = new Parameter.Default(new double[]{0.0, 0.0, 0.0, 3.0, 50, 0.0});
        strengthOfSelectionParameters[5] = new Parameter.Default(new double[]{2.0, 0.0, 0.0, 0.0, 0.0, 60});
        MatrixParameter strengthOfSelectionMatrixParam
                = new MatrixParameter("strengthOfSelectionMatrix", strengthOfSelectionParameters);

        // Wrt Precision
        DiffusionProcessDelegate diffusionProcessDelegate
                = new OUDiffusionModelDelegate(treeModel, diffusionModel,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParam));
        System.out.println("\nTest OU gradient precision.");
        testGradient(diffusionModel, diffusionProcessDelegate, dataModel, precisionMatrix, false);
        System.out.println("\nTest OU gradient precision with missing.");
        testGradient(diffusionModel, diffusionProcessDelegate, dataModelMissing, precisionMatrix, false);

        // Wrt Variance
        DiffusionProcessDelegate diffusionProcessDelegateVariance
                = new OUDiffusionModelDelegate(treeModel, diffusionModelVar,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        System.out.println("\nTest OU gradient variance.");
        testGradient(diffusionModelVar, diffusionProcessDelegateVariance, dataModel, precisionMatrixInv, false);
        System.out.println("\nTest OU gradient variance with missing.");
        testGradient(diffusionModelVar, diffusionProcessDelegateVariance, dataModelMissing, precisionMatrixInv, false);

    }

    private void testGradient(MultivariateDiffusionModel diffusionModel,
                              DiffusionProcessDelegate diffusionProcessDelegate,
                              ContinuousTraitPartialsProvider dataModel,
                              MatrixParameterInterface precision,
                              Boolean wishart) {
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
                        ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter.WRT_VARIANCE);
        BranchSpecificGradient branchSpecificGradient =
                new BranchSpecificGradient("trait", dataLikelihood, cdld, traitGradient, precision);

        GradientWrtPrecisionProvider gPPBranchSpecific = new GradientWrtPrecisionProvider.BranchSpecificGradientWrtPrecisionProvider(branchSpecificGradient);

        // Correlation Gradient Branch Specific
        CorrelationPrecisionGradient gradientProviderBranchSpecific = new CorrelationPrecisionGradient(gPPBranchSpecific, dataLikelihood, precision);

        String sBS = gradientProviderBranchSpecific.getReport();
        System.err.println(sBS);
        double[] gradientAnalyticalBS = parseGradient(sBS, "analytic");
        double[] gradientNumeric = parseGradient(sBS, "numeric (with Cholesky):");

        assertEquals("Sizes", gradientAnalyticalBS.length, gradientNumeric.length);

        for (int k = 0; k < gradientAnalyticalBS.length; k++) {
            assertEquals("gradient correlation k=" + k,
                    gradientAnalyticalBS[k],
                    gradientNumeric[k],
                    delta);
        }

        // Diagonal Gradient Branch Specific
        DiagonalPrecisionGradient gradientDiagonalProviderBS = new DiagonalPrecisionGradient(gPPBranchSpecific, dataLikelihood, precision);

        String sDiagBS = gradientDiagonalProviderBS.getReport();
        System.err.println(sDiagBS);
        double[] gradientDiagonalAnalyticalBS = parseGradient(sDiagBS, "analytic");
        double[] gradientDiagonalNumeric = parseGradient(sDiagBS, "numeric:");

        assertEquals("Sizes", gradientDiagonalAnalyticalBS.length, gradientDiagonalNumeric.length);

        for (int k = 0; k < gradientDiagonalAnalyticalBS.length; k++) {
            assertEquals("gradient correlation k=" + k,
                    gradientDiagonalAnalyticalBS[k],
                    gradientDiagonalNumeric[k],
                    delta);
        }

        if (wishart) {
            // Wishart Statistic
            WishartStatisticsWrapper wishartStatistics
                    = new WishartStatisticsWrapper("wishart", "trait", dataLikelihood,
                    cdld);

            GradientWrtPrecisionProvider gPPWiwhart = new GradientWrtPrecisionProvider.WishartGradientWrtPrecisionProvider(wishartStatistics);

            // Correlation Gradient
            CorrelationPrecisionGradient gradientProviderWishart = new CorrelationPrecisionGradient(gPPWiwhart, dataLikelihood, precision);

            String sW = gradientProviderWishart.getReport();
            System.err.println(sW);
            double[] gradientAnalyticalW = parseGradient(sW, "analytic");
            assertEquals("Sizes", gradientAnalyticalW.length, gradientNumeric.length);

            for (int k = 0; k < gradientAnalyticalW.length; k++) {
                assertEquals("gradient correlation k=" + k,
                        gradientAnalyticalW[k],
                        gradientAnalyticalBS[k],
                        delta);
            }

            // Diagonal Gradient
            DiagonalPrecisionGradient gradientDiagonalProviderW = new DiagonalPrecisionGradient(gPPWiwhart, dataLikelihood, precision);

            String sDiagW = gradientDiagonalProviderW.getReport();
            System.err.println(sDiagW);
            double[] gradientDiagonalAnalyticalW = parseGradient(sDiagW, "analytic");

            assertEquals("Sizes", gradientDiagonalAnalyticalW.length, gradientDiagonalNumeric.length);

            for (int k = 0; k < gradientDiagonalAnalyticalW.length; k++) {
                assertEquals("gradient diagonal k=" + k,
                        gradientDiagonalAnalyticalW[k],
                        gradientDiagonalAnalyticalBS[k],
                        delta);
            }
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
}
