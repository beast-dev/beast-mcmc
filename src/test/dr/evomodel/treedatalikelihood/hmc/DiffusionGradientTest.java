/*
 * DiffusionGradientTest.java
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
import dr.evomodel.treedatalikelihood.hmc.AbstractDiffusionGradient;
import dr.evomodel.treedatalikelihood.hmc.CorrelationPrecisionGradient;
import dr.evomodel.treedatalikelihood.hmc.DiagonalPrecisionGradient;
import dr.evomodel.treedatalikelihood.hmc.GradientWrtPrecisionProvider;
import dr.evomodel.treedatalikelihood.preorder.*;
import dr.inference.model.*;
import test.dr.evomodel.treedatalikelihood.continuous.ContinuousTraitTest;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static dr.evomodel.treedatalikelihood.continuous.ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient;
import static dr.evomodel.treedatalikelihood.continuous.ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter;
import static dr.evomodel.treedatalikelihood.hmc.AbstractDiffusionGradient.ParameterDiffusionGradient.createDiagonalAttenuationGradient;
import static dr.evomodel.treedatalikelihood.hmc.AbstractDiffusionGradient.ParameterDiffusionGradient.createDriftGradient;

public class DiffusionGradientTest extends ContinuousTraitTest {

    private CompoundSymmetricMatrix precisionMatrix;
    private CachedMatrixInverse precisionMatrixInv;

    //    protected List<Integer> missingIndices = new ArrayList<Integer>();
    protected boolean[] missingIndicators;

    private MultivariateDiffusionModel diffusionModelVar;
    private ContinuousTraitPartialsProvider dataModelMissing;

    private MatrixParameter samplingPrecision;
    private CachedMatrixInverse samplingPrecisionInv;
    private RepeatedMeasuresTraitDataModel dataModelRepeatedMeasures;
    private RepeatedMeasuresTraitDataModel dataModelRepeatedMeasuresInv;

    private CompoundParameter meanRoot;

    private Boolean fixedRoot = true;

    private double delta;

    public DiffusionGradientTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        delta = 1E-2;

        dimTrait = 6;

        Parameter offDiagonal = new Parameter.Default(new double[]{0.12, -0.13, 0.14, -0.15, 0.16,
                -0.12, 0.13, -0.14, 0.15,
                0.12, -0.13, 0.14,
                -0.12, 0.13,
                0.12});
        offDiagonal.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, dimTrait * (dimTrait - 1) / 2));

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
        traitParameter = new CompoundParameter("trait", dataTraits);

        this.missingIndicators = new boolean[traitParameter.getDimension()];

        traitParameter.setParameterValue(2, 0);

        missingIndicators[6] = true;
        missingIndicators[7] = true;
        missingIndicators[8] = true;
        missingIndicators[9] = true;
        missingIndicators[10] = true;
        missingIndicators[11] = true;
        missingIndicators[13] = true;
        missingIndicators[15] = true;
        missingIndicators[25] = true;
        missingIndicators[29] = true;


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
        final double rootVal = fixedRoot ? Double.POSITIVE_INFINITY : 0.1;
        Parameter rootVar = new Parameter.Default("rootVar", rootVal);
        Parameter[] meanRootList = new Parameter[dimTrait];
        double[] meanRootVal = new double[]{-1.0, -3.0, 2.5, -2.5, 1.3, 4.0};
        for (int i = 0; i < dimTrait; i++) {
            meanRootList[i] = new Parameter.Default("meanRoot." + (i + 1), meanRootVal[i]);
        }
        meanRoot = new CompoundParameter("rootMean", meanRootList);
        rootPrior = new ConjugateRootTraitPrior(meanRoot, rootVar);

        // Data Model
        dataModelMissing = new ContinuousTraitDataModel("dataModel",
                traitParameter,
                missingIndicators, true,
                6, precisionType);

        dataModel = new ContinuousTraitDataModel("dataModel",
                traitParameter,
                missingIndicators, false,
                6, precisionType);

        //// Factor Model //// *****************************************************************************************
        // Error model
        Parameter factorPrecisionParameters = new Parameter.Default("factorPrecision", new double[]{1.0, 5.0, 0.5, 0.1, 0.2, 0.3});

        // Loadings
        Parameter[] loadingsParameters = new Parameter[2];
        loadingsParameters[0] = new Parameter.Default(new double[]{1.0, 2.0, 3.0, 0.0, 0.0, 0.0});
        loadingsParameters[1] = new Parameter.Default(new double[]{0.0, 0.0, 0.0, 0.0, 0.5, 1.0});
        MatrixParameterInterface loadingsMatrixParameters = new MatrixParameter("loadings", loadingsParameters);

        dataModelFactor = new IntegratedFactorAnalysisLikelihood("dataModelFactors",
                traitParameter,
                missingIndicators,
                loadingsMatrixParameters,
                factorPrecisionParameters, 0.0, null,
                IntegratedFactorAnalysisLikelihood.CacheProvider.NO_CACHE);

        //// Repeated Measures Model //// *****************************************************************************************
        Parameter offDiagonalSampling = new Parameter.Default(new double[]{0.16, -0.15, 0.14, -0.13, 0.12,
                -0.15, 0.14, -0.13, 0.12,
                0.0, 0.0, 0.0,
                -0.13, 0.12,
                0.0});

        offDiagonalSampling.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, dimTrait * (dimTrait - 1) / 2));

        Parameter diagonalSampling = new Parameter.Default(new double[]{0.1, 0.2, 0.3, 0.2, 0.1, 0.01});

        samplingPrecision = new CompoundSymmetricMatrix(diagonalSampling, offDiagonalSampling, true, true);

        Parameter diagonalVarSampling = new Parameter.Default(new double[]{10.0, 20.0, 30.0, 20.0, 10.0, 100.0});

        samplingPrecisionInv = new CachedMatrixInverse("samplingVar",
                new CompoundSymmetricMatrix(diagonalVarSampling, offDiagonalSampling, true, false));

        dataModelRepeatedMeasures = new RepeatedMeasuresTraitDataModel("dataModelRepeatedMeasures",
                dataModel,
                traitParameter,
                missingIndicators,
                true,
                dimTrait,
                1,
                samplingPrecision,
                PrecisionType.FULL);

        dataModelRepeatedMeasuresInv = new RepeatedMeasuresTraitDataModel("dataModelRepeatedMeasuresInv",
                dataModel,
                traitParameter,
                missingIndicators,
                true,
                dimTrait,
                1,
                samplingPrecisionInv,
                PrecisionType.FULL);

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

        // Factor Model
        DiffusionProcessDelegate diffusionProcessDelegateFactor
                = new HomogeneousDiffusionModelDelegate(treeModel, diffusionModelFactor);
        System.out.println("\nTest gradient precision factor.");
        testGradient(diffusionModelFactor, diffusionProcessDelegateFactor, dataModelFactor, rootPriorFactor, rootMeanFactor, precisionMatrixFactor, false);

        // Repeated Measures Model
        System.out.println("\nTest gradient precision repeated measures.");
        testGradient(diffusionModel, diffusionProcessDelegate, dataModelRepeatedMeasures, rootPrior, meanRoot, precisionMatrix, false, null, null, samplingPrecision);
        testGradient(diffusionModel, diffusionProcessDelegate, dataModelRepeatedMeasuresInv, rootPrior, meanRoot, precisionMatrix, false, null, null, samplingPrecisionInv);
    }

    public void testGradientDriftWithMissing() {

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
        testGradient(diffusionModel, diffusionProcessDelegate, dataModel, precisionMatrix);
        System.out.println("\nTest drift gradient precision with missing.");
        testGradient(diffusionModel, diffusionProcessDelegate, dataModelMissing, precisionMatrix);

        // Wrt Variance
        DiffusionProcessDelegate diffusionProcessDelegateVariance
                = new DriftDiffusionModelDelegate(treeModel, diffusionModelVar, driftModels);
        System.out.println("\nTest drift gradient variance.");
        testGradient(diffusionModelVar, diffusionProcessDelegateVariance, dataModel, precisionMatrixInv);
        System.out.println("\nTest drift gradient variance with missing.");
        testGradient(diffusionModelVar, diffusionProcessDelegateVariance, dataModelMissing, precisionMatrixInv);

        // Factor Model
        // Diffusion
        List<BranchRateModel> driftModelsFactor = new ArrayList<BranchRateModel>();
        driftModelsFactor.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{0.0})));
        driftModelsFactor.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{200.0})));

        DiffusionProcessDelegate diffusionProcessDelegateFactor
                = new DriftDiffusionModelDelegate(treeModel, diffusionModelFactor, driftModelsFactor);
        System.out.println("\nTest gradient precision.");
        testGradient(diffusionModelFactor, diffusionProcessDelegateFactor, dataModelFactor, rootPriorFactor, rootMeanFactor, precisionMatrixFactor, false);

        // Repeated Measures Model
        System.out.println("\nTest gradient precision repeated measures.");
        testGradient(diffusionModel, diffusionProcessDelegate, dataModelRepeatedMeasures, rootPrior, meanRoot, precisionMatrix, false, null, null, samplingPrecision);
        testGradient(diffusionModel, diffusionProcessDelegate, dataModelRepeatedMeasuresInv, rootPrior, meanRoot, precisionMatrix, false, null, null, samplingPrecisionInv);
    }

    public void testGradientSingleDriftWithMissing() {

        // Diffusion
        List<BranchRateModel> driftModels = new ArrayList<BranchRateModel>();
        CompoundParameter driftParam = new CompoundParameter("drift");
        for (int i = 0; i < dimTrait; i++) {
            Parameter rate = new Parameter.Default("rate." + (i + 1), new double[]{2.2 * i + 3.1});
            driftParam.addParameter(rate);
            driftModels.add(new StrictClockBranchRates(rate));
        }

        // Wrt Precision
        DiffusionProcessDelegate diffusionProcessDelegate
                = new DriftDiffusionModelDelegate(treeModel, diffusionModel, driftModels);
        System.out.println("\nTest single drift gradient drift.");
        testGradient(diffusionModel, diffusionProcessDelegate, dataModel, precisionMatrix, driftParam);
        System.out.println("\nTest single drift gradient precision with missing.");
        testGradient(diffusionModel, diffusionProcessDelegate, dataModelMissing, precisionMatrix, driftParam);

        // Wrt Variance
        DiffusionProcessDelegate diffusionProcessDelegateVariance
                = new DriftDiffusionModelDelegate(treeModel, diffusionModelVar, driftModels);
        System.out.println("\nTest single drift gradient variance.");
        testGradient(diffusionModelVar, diffusionProcessDelegateVariance, dataModel, precisionMatrixInv, driftParam);
        System.out.println("\nTest single drift gradient variance with missing.");
        testGradient(diffusionModelVar, diffusionProcessDelegateVariance, dataModelMissing, precisionMatrixInv, driftParam);

        // Repeated Measures Model
        System.out.println("\nTest gradient precision repeated measures.");
        testGradient(diffusionModel, diffusionProcessDelegate, dataModelRepeatedMeasures, rootPrior, meanRoot, precisionMatrix, false, null, driftParam, samplingPrecision);
        testGradient(diffusionModel, diffusionProcessDelegate, dataModelRepeatedMeasuresInv, rootPrior, meanRoot, precisionMatrix, false, null, driftParam, samplingPrecisionInv);
    }

    public void testGradientSingleDriftSameMeanWithMissing() {

        // Diffusion
        List<BranchRateModel> driftModels = new ArrayList<BranchRateModel>();
        for (int i = 0; i < dimTrait; i++) {
            driftModels.add(new StrictClockBranchRates(meanRoot.getParameter(i)));
        }

        // Wrt Precision
        DiffusionProcessDelegate diffusionProcessDelegate
                = new DriftDiffusionModelDelegate(treeModel, diffusionModel, driftModels);
        System.out.println("\nTest single drift same root gradient drift.");
        testGradient(diffusionModel, diffusionProcessDelegate, dataModel, precisionMatrix, meanRoot);
        System.out.println("\nTest single drift same root gradient precision with missing.");
        testGradient(diffusionModel, diffusionProcessDelegate, dataModelMissing, precisionMatrix, meanRoot);

        // Wrt Variance
        DiffusionProcessDelegate diffusionProcessDelegateVariance
                = new DriftDiffusionModelDelegate(treeModel, diffusionModelVar, driftModels);
        System.out.println("\nTest single drift same root gradient variance.");
        testGradient(diffusionModelVar, diffusionProcessDelegateVariance, dataModel, precisionMatrixInv, meanRoot);
        System.out.println("\nTest single drift same root gradient variance with missing.");
        testGradient(diffusionModelVar, diffusionProcessDelegateVariance, dataModelMissing, precisionMatrixInv, meanRoot);

        // Repeated Measures Model
        System.out.println("\nTest gradient precision repeated measures.");
        testGradient(diffusionModel, diffusionProcessDelegate, dataModelRepeatedMeasures, rootPrior, meanRoot, precisionMatrix, false, null, meanRoot, samplingPrecision);
        testGradient(diffusionModel, diffusionProcessDelegate, dataModelRepeatedMeasuresInv, rootPrior, meanRoot, precisionMatrix, false, null, meanRoot, samplingPrecisionInv);
    }

    public void testGradientOUWithMissing() {

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
        testGradient(diffusionModel, diffusionProcessDelegate, dataModel, precisionMatrix);
        System.out.println("\nTest OU gradient precision with missing.");
        testGradient(diffusionModel, diffusionProcessDelegate, dataModelMissing, precisionMatrix);

        // Wrt Variance
        DiffusionProcessDelegate diffusionProcessDelegateVariance
                = new OUDiffusionModelDelegate(treeModel, diffusionModelVar,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        System.out.println("\nTest OU gradient variance.");
        testGradient(diffusionModelVar, diffusionProcessDelegateVariance, dataModel, precisionMatrixInv);
        System.out.println("\nTest OU gradient variance with missing.");
        testGradient(diffusionModelVar, diffusionProcessDelegateVariance, dataModelMissing, precisionMatrixInv);

        // Factor Model
        List<BranchRateModel> optimalTraitsModelsFactor = new ArrayList<BranchRateModel>();
        optimalTraitsModelsFactor.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{0.0})));
        optimalTraitsModelsFactor.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{2.0})));

        Parameter[] strengthOfSelectionParametersFactor = new Parameter[2];
        strengthOfSelectionParametersFactor[0] = new Parameter.Default(new double[]{10, 1.0});
        strengthOfSelectionParametersFactor[1] = new Parameter.Default(new double[]{1.0, 20});
        MatrixParameter strengthOfSelectionMatrixParamFactor
                = new MatrixParameter("strengthOfSelectionMatrix", strengthOfSelectionParametersFactor);

        DiffusionProcessDelegate diffusionProcessDelegateFactor
                = new OUDiffusionModelDelegate(treeModel, diffusionModelFactor,
                optimalTraitsModelsFactor, new MultivariateElasticModel(strengthOfSelectionMatrixParamFactor));
        System.out.println("\nTest gradient precision.");
        testGradient(diffusionModelFactor, diffusionProcessDelegateFactor, dataModelFactor, rootPriorFactor, rootMeanFactor, precisionMatrixFactor, false);

        // Repeated Measures Model
        System.out.println("\nTest gradient precision repeated measures.");
        testGradient(diffusionModel, diffusionProcessDelegate, dataModelRepeatedMeasures, rootPrior, meanRoot, precisionMatrix, false, null, null, samplingPrecision);
        testGradient(diffusionModel, diffusionProcessDelegate, dataModelRepeatedMeasuresInv, rootPrior, meanRoot, precisionMatrix, false, null, null, samplingPrecisionInv);

        //************//
        // Single opt
        List<BranchRateModel> optimalTraitsModelsSingle = new ArrayList<BranchRateModel>();
        CompoundParameter optParamSingle = new CompoundParameter("opt");
        for (int i = 0; i < dimTrait; i++) {
            Parameter rate = new Parameter.Default("opt." + (i + 1), new double[]{2.2 * i + 3.1});
            optParamSingle.addParameter(rate);
            optimalTraitsModelsSingle.add(new StrictClockBranchRates(rate));
        }

        // Wrt Precision
        DiffusionProcessDelegate diffusionProcessDelegateSingle
                = new OUDiffusionModelDelegate(treeModel, diffusionModel,
                optimalTraitsModelsSingle, new MultivariateElasticModel(strengthOfSelectionMatrixParam));
        System.out.println("\nTest OU single opt gradient precision.");
        testGradient(diffusionModel, diffusionProcessDelegateSingle, dataModel, precisionMatrix, optParamSingle);
        System.out.println("\nTest OU single opt gradient precision with missing.");
        testGradient(diffusionModel, diffusionProcessDelegateSingle, dataModelMissing, precisionMatrix, optParamSingle);

        // Wrt Variance
        DiffusionProcessDelegate diffusionProcessDelegateVarianceSingle
                = new OUDiffusionModelDelegate(treeModel, diffusionModelVar,
                optimalTraitsModelsSingle, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        System.out.println("\nTest OU single opt gradient variance.");
        testGradient(diffusionModelVar, diffusionProcessDelegateVarianceSingle, dataModel, precisionMatrixInv, optParamSingle);
        System.out.println("\nTest OU single opt gradient variance with missing.");
        testGradient(diffusionModelVar, diffusionProcessDelegateVarianceSingle, dataModelMissing, precisionMatrixInv, optParamSingle);

        // Repeated Measures Model
        System.out.println("\nTest gradient precision repeated measures.");
        testGradient(diffusionModel, diffusionProcessDelegateSingle, dataModelRepeatedMeasures, rootPrior, meanRoot, precisionMatrix, false, null, optParamSingle, samplingPrecision);
        testGradient(diffusionModel, diffusionProcessDelegateSingle, dataModelRepeatedMeasuresInv, rootPrior, meanRoot, precisionMatrix, false, null, optParamSingle, samplingPrecisionInv);

        //**************//
        // Same mean
        // Diffusion
        List<BranchRateModel> optimalTraitsModelsSame = new ArrayList<BranchRateModel>();
        for (int i = 0; i < dimTrait; i++) {
            optimalTraitsModelsSame.add(new StrictClockBranchRates(meanRoot.getParameter(i)));
        }

        // Wrt Precision
        DiffusionProcessDelegate diffusionProcessDelegateSame
                = new OUDiffusionModelDelegate(treeModel, diffusionModel,
                optimalTraitsModelsSame, new MultivariateElasticModel(strengthOfSelectionMatrixParam));
        System.out.println("\nTest OU Same opt gradient precision.");
        testGradient(diffusionModel, diffusionProcessDelegateSame, dataModel, precisionMatrix, meanRoot);
        System.out.println("\nTest OU Same opt gradient precision with missing.");
        testGradient(diffusionModel, diffusionProcessDelegateSame, dataModelMissing, precisionMatrix, meanRoot);

        // Wrt Variance
        DiffusionProcessDelegate diffusionProcessDelegateVarianceSame
                = new OUDiffusionModelDelegate(treeModel, diffusionModelVar,
                optimalTraitsModelsSame, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        System.out.println("\nTest OU Same opt gradient variance.");
        testGradient(diffusionModelVar, diffusionProcessDelegateVarianceSame, dataModel, precisionMatrixInv, meanRoot);
        System.out.println("\nTest OU Same opt gradient variance with missing.");
        testGradient(diffusionModelVar, diffusionProcessDelegateVarianceSame, dataModelMissing, precisionMatrixInv, meanRoot);

        // Repeated Measures Model
        System.out.println("\nTest gradient precision repeated measures.");
        testGradient(diffusionModel, diffusionProcessDelegateSame, dataModelRepeatedMeasures, rootPrior, meanRoot, precisionMatrix, false, null, meanRoot, samplingPrecision);
        testGradient(diffusionModel, diffusionProcessDelegateSame, dataModelRepeatedMeasuresInv, rootPrior, meanRoot, precisionMatrix, false, null, meanRoot, samplingPrecisionInv);
    }

    public void testGradientDiagonalOUWithMissing() {

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{1.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{2.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.3", new double[]{-2.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.4", new double[]{1.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.5", new double[]{2.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.6", new double[]{-2.0})));

        DiagonalMatrix strengthOfSelectionMatrixParam
                = new DiagonalMatrix(new Parameter.Default(new double[]{0.0, 0.000005, 1.0, 5.0, 10.0, 50.0}));

        // Wrt Precision
        DiffusionProcessDelegate diffusionProcessDelegate
                = new OUDiffusionModelDelegate(treeModel, diffusionModel,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParam));
        System.out.println("\nTest Diagonal OU gradient precision.");
        testGradient(diffusionModel, diffusionProcessDelegate, dataModel, precisionMatrix, strengthOfSelectionMatrixParam);
        System.out.println("\nTest Diagonal OU gradient precision with missing.");
        testGradient(diffusionModel, diffusionProcessDelegate, dataModelMissing, precisionMatrix, strengthOfSelectionMatrixParam);

        // Wrt Variance
        DiffusionProcessDelegate diffusionProcessDelegateVariance
                = new OUDiffusionModelDelegate(treeModel, diffusionModelVar,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        System.out.println("\nTest Diagonal OU gradient variance.");
        testGradient(diffusionModelVar, diffusionProcessDelegateVariance, dataModel, precisionMatrixInv, strengthOfSelectionMatrixParam);
        System.out.println("\nTest Diagonal OU gradient variance with missing.");
        testGradient(diffusionModelVar, diffusionProcessDelegateVariance, dataModelMissing, precisionMatrixInv, strengthOfSelectionMatrixParam);

        // Factor Model
        List<BranchRateModel> optimalTraitsModelsFactor = new ArrayList<BranchRateModel>();
        optimalTraitsModelsFactor.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{0.0})));
        optimalTraitsModelsFactor.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{2.0})));

        DiagonalMatrix strengthOfSelectionMatrixParamFactor
                = new DiagonalMatrix(new Parameter.Default(new double[]{0.0, 50.0}));

        DiffusionProcessDelegate diffusionProcessDelegateFactor
                = new OUDiffusionModelDelegate(treeModel, diffusionModelFactor,
                optimalTraitsModelsFactor, new MultivariateElasticModel(strengthOfSelectionMatrixParamFactor));
        System.out.println("\nTest gradient precision.");
        testGradient(diffusionModelFactor, diffusionProcessDelegateFactor, dataModelFactor, rootPriorFactor, rootMeanFactor, precisionMatrixFactor, false, strengthOfSelectionMatrixParamFactor, null, null);

        // Repeated Measures Model
        System.out.println("\nTest gradient precision repeated measures.");
        testGradient(diffusionModel, diffusionProcessDelegate, dataModelRepeatedMeasures, rootPrior, meanRoot, precisionMatrix, false, strengthOfSelectionMatrixParam, null, samplingPrecision);
        testGradient(diffusionModel, diffusionProcessDelegate, dataModelRepeatedMeasuresInv, rootPrior, meanRoot, precisionMatrix, false, strengthOfSelectionMatrixParam, null, samplingPrecisionInv);
    }

    private void testGradient(MultivariateDiffusionModel diffusionModel,
                              DiffusionProcessDelegate diffusionProcessDelegate,
                              ContinuousTraitPartialsProvider dataModel,
                              MatrixParameterInterface precision,
                              Boolean wishart) {
        testGradient(diffusionModel, diffusionProcessDelegate, dataModel, rootPrior, meanRoot, precision, wishart, null, null, null);
    }

    private void testGradient(MultivariateDiffusionModel diffusionModel,
                              DiffusionProcessDelegate diffusionProcessDelegate,
                              ContinuousTraitPartialsProvider dataModel,
                              ConjugateRootTraitPrior rootPrior,
                              Parameter meanRoot,
                              MatrixParameterInterface precision,
                              Boolean wishart) {
        testGradient(diffusionModel, diffusionProcessDelegate, dataModel, rootPrior, meanRoot, precision, wishart, null, null, null);
    }

    private void testGradient(MultivariateDiffusionModel diffusionModel,
                              DiffusionProcessDelegate diffusionProcessDelegate,
                              ContinuousTraitPartialsProvider dataModel,
                              MatrixParameterInterface precision) {
        testGradient(diffusionModel, diffusionProcessDelegate, dataModel, rootPrior, meanRoot, precision, false, null, null, null);
    }

    private void testGradient(MultivariateDiffusionModel diffusionModel,
                              DiffusionProcessDelegate diffusionProcessDelegate,
                              ContinuousTraitPartialsProvider dataModel,
                              MatrixParameterInterface precision,
                              Parameter driftParam) {
        testGradient(diffusionModel, diffusionProcessDelegate, dataModel, rootPrior, meanRoot, precision, false, null, driftParam, null);
    }

    private void testGradient(MultivariateDiffusionModel diffusionModel,
                              DiffusionProcessDelegate diffusionProcessDelegate,
                              ContinuousTraitPartialsProvider dataModel,
                              MatrixParameterInterface precision,
                              MatrixParameterInterface attenuation,
                              Boolean wishart) {
        testGradient(diffusionModel, diffusionProcessDelegate, dataModel, rootPrior, meanRoot, precision, wishart, attenuation, null, null);
    }

    private void testGradient(MultivariateDiffusionModel diffusionModel,
                              DiffusionProcessDelegate diffusionProcessDelegate,
                              ContinuousTraitPartialsProvider dataModel,
                              MatrixParameterInterface precision,
                              MatrixParameterInterface attenuation) {
        testGradient(diffusionModel, diffusionProcessDelegate, dataModel, rootPrior, meanRoot, precision, false, attenuation, null, null);
    }

    private void testGradient(MultivariateDiffusionModel diffusionModel,
                              DiffusionProcessDelegate diffusionProcessDelegate,
                              ContinuousTraitPartialsProvider dataModel,
                              MatrixParameterInterface precision,
                              MatrixParameterInterface attenuation,
                              Parameter drift) {
        testGradient(diffusionModel, diffusionProcessDelegate, dataModel, rootPrior, meanRoot, precision, false, attenuation, drift, null);
    }

    private void testGradient(MultivariateDiffusionModel diffusionModel,
                              DiffusionProcessDelegate diffusionProcessDelegate,
                              ContinuousTraitPartialsProvider dataModel,
                              ConjugateRootTraitPrior rootPrior,
                              Parameter meanRoot,
                              MatrixParameterInterface precision, Boolean wishart,
                              MatrixParameterInterface attenuation,
                              Parameter drift,
                              MatrixParameterInterface samplingPrecision) {
        int dimLocal = rootPrior.getMean().length;

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

        // Variance
        ContinuousDataLikelihoodDelegate cdld = (ContinuousDataLikelihoodDelegate) dataLikelihood.getDataLikelihoodDelegate();

        if (precision != null) {
            // Branch Specific
            ContinuousProcessParameterGradient traitGradient =
                    new ContinuousProcessParameterGradient(
                            rootPrior.getMean().length, treeModel, cdld,
                            new ArrayList<>(
                                    Arrays.asList(DerivationParameter.WRT_VARIANCE)
                            ));

            BranchSpecificGradient branchSpecificGradient =
                    new BranchSpecificGradient("trait", dataLikelihood, cdld, traitGradient, precision);

            GradientWrtPrecisionProvider gPPBranchSpecific = new GradientWrtPrecisionProvider.BranchSpecificGradientWrtPrecisionProvider(branchSpecificGradient);

            // Correlation Gradient Branch Specific
            CorrelationPrecisionGradient gradientProviderBranchSpecific = new CorrelationPrecisionGradient(gPPBranchSpecific, dataLikelihood, precision);

            double[] gradientAnalyticalBS = testOneGradient(gradientProviderBranchSpecific);

            // Diagonal Gradient Branch Specific
            DiagonalPrecisionGradient gradientDiagonalProviderBS = new DiagonalPrecisionGradient(gPPBranchSpecific, dataLikelihood, precision);

            double[] gradientDiagonalAnalyticalBS = testOneGradient(gradientDiagonalProviderBS);

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
                assertEquals("Sizes", gradientAnalyticalW.length, gradientAnalyticalBS.length);

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

                assertEquals("Sizes", gradientDiagonalAnalyticalW.length, gradientDiagonalAnalyticalBS.length);

                for (int k = 0; k < gradientDiagonalAnalyticalW.length; k++) {
                    assertEquals("gradient diagonal k=" + k,
                            gradientDiagonalAnalyticalW[k],
                            gradientDiagonalAnalyticalBS[k],
                            delta);
                }
            }
        }

        // Diagonal Attenuation Gradient Branch Specific
        if (attenuation != null) {
            ContinuousProcessParameterGradient traitGradientAtt =
                    new ContinuousProcessParameterGradient(
                            dimLocal, treeModel, cdld,
                            new ArrayList<>(
                                    Arrays.asList(DerivationParameter.WRT_DIAGONAL_SELECTION_STRENGTH)
                            ));

            BranchSpecificGradient branchSpecificGradientAtt =
                    new BranchSpecificGradient("trait", dataLikelihood, cdld, traitGradientAtt, attenuation);

            AbstractDiffusionGradient.ParameterDiffusionGradient gABranchSpecific = createDiagonalAttenuationGradient(branchSpecificGradientAtt, dataLikelihood, attenuation);
            testOneGradient(gABranchSpecific);

        }

        // WRT root mean
        boolean sameRoot = (drift == meanRoot);
        ContinuousProcessParameterGradient traitGradientRoot =
                new ContinuousProcessParameterGradient(
                        dimLocal, treeModel, cdld,
                        new ArrayList<>(
                                Arrays.asList(sameRoot ? DerivationParameter.WRT_CONSTANT_DRIFT_AND_ROOT_MEAN : DerivationParameter.WRT_ROOT_MEAN)
                        ));

        BranchSpecificGradient branchSpecificGradientRoot =
                new BranchSpecificGradient("trait", dataLikelihood, cdld, traitGradientRoot, meanRoot);

        AbstractDiffusionGradient.ParameterDiffusionGradient gRootBranchSpecific = createDriftGradient(branchSpecificGradientRoot, dataLikelihood, meanRoot);
        testOneGradient(gRootBranchSpecific);

        // Drift Gradient Branch Specific
        if (drift != null && !sameRoot) {
            ContinuousProcessParameterGradient traitGradientDrift =
                    new ContinuousProcessParameterGradient(
                            dimLocal, treeModel, cdld,
                            new ArrayList<>(
                                    Arrays.asList(DerivationParameter.WRT_CONSTANT_DRIFT)
                            ));

            BranchSpecificGradient branchSpecificGradientDrift =
                    new BranchSpecificGradient("trait", dataLikelihood, cdld, traitGradientDrift, drift);

            AbstractDiffusionGradient.ParameterDiffusionGradient gDriftBranchSpecific = createDriftGradient(branchSpecificGradientDrift, dataLikelihood, drift);
            testOneGradient(gDriftBranchSpecific);
        }

        // Sampling Precision
        if (samplingPrecision != null) {
            ContinuousTraitGradientForBranch.SamplingVarianceGradient traitGradientSampling =
                    new ContinuousTraitGradientForBranch.SamplingVarianceGradient(
                            dimLocal, treeModel, likelihoodDelegate,
                            (ModelExtensionProvider.NormalExtensionProvider) dataModel);

            BranchSpecificGradient branchSpecificGradientSampling =
                    new BranchSpecificGradient("trait", dataLikelihood, cdld, traitGradientSampling, samplingPrecision);

            GradientWrtPrecisionProvider gPPBranchSpecificSampling = new GradientWrtPrecisionProvider.BranchSpecificGradientWrtPrecisionProvider(branchSpecificGradientSampling);

            // Correlation Gradient Branch Specific
//            CorrelationPrecisionGradient gradientProviderBranchSpecificSampling = new CorrelationPrecisionGradient(gPPBranchSpecificSampling, dataLikelihood, samplingPrecision);
//
//            testOneGradient(gradientProviderBranchSpecificSampling);

            // Diagonal Gradient Branch Specific
            DiagonalPrecisionGradient gradientDiagonalProviderBSSampling = new DiagonalPrecisionGradient(gPPBranchSpecificSampling, dataLikelihood, samplingPrecision);

            testOneGradient(gradientDiagonalProviderBSSampling);
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

    private double[] testOneGradient(AbstractDiffusionGradient gradientProviderBranchSpecific) {
        String sBS = gradientProviderBranchSpecific.getReport();
        System.err.println(sBS);
        double[] gradientAnalyticalBS = parseGradient(sBS, "analytic");
        double[] gradientNumeric = parseGradient(sBS, "numeric :");

        assertEquals("Sizes", gradientAnalyticalBS.length, gradientNumeric.length);

        for (int k = 0; k < gradientAnalyticalBS.length; k++) {
            assertEquals("gradient k=" + k,
                    gradientAnalyticalBS[k],
                    gradientNumeric[k],
                    delta);
        }

        return gradientAnalyticalBS;
    }
}

