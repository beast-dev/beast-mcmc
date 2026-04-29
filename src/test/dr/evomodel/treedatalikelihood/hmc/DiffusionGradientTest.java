/*
 * DiffusionGradientTest.java
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

package test.dr.evomodel.treedatalikelihood.hmc;

import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.TreeTrait;
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
import dr.evomodel.treedatalikelihood.preorder.BranchConditionalDistributionDelegate;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.evomodel.treedatalikelihood.preorder.*;
import dr.inference.model.*;
import dr.math.MultivariateFunction;
import test.dr.evomodel.treedatalikelihood.continuous.ContinuousTraitTest;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static dr.evomodel.treedatalikelihood.continuous.ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient;
import static dr.evomodel.treedatalikelihood.continuous.ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter;
import static dr.evomodel.treedatalikelihood.hmc.AbstractDiffusionGradient.ParameterDiffusionGradient.createDiagonalAttenuationGradient;
import static dr.evomodel.treedatalikelihood.hmc.AbstractDiffusionGradient.ParameterDiffusionGradient.createDriftGradient;

public class DiffusionGradientTest extends ContinuousTraitTest {

    private static final double ANALYTIC_REGRESSION_TOLERANCE = 1e-8;

    /**
     * Fixed regression baseline for OU (dense selection matrix) precision gradient
     * in the non-missing setup.
     */
    private static final double[] EXPECTED_OU_NON_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            -2012.127519876533, -3787.7530373745367, -3426.214318032913,
            -7036.583460241445, -4018.096502174195, -4183.1564069705055,
            -5373.7858514049585, -8907.56103481203, -5465.185580098292,
            -8065.756829069363, -10510.668358646974, -10767.513324782352,
            -14240.92412646215, -15551.55947436676, -16941.61261285212
    };

    /**
     * Fixed regression baseline for diagonal-OU precision gradient
     * in the non-missing setup.
     */
    private static final double[] EXPECTED_DIAGONAL_OU_NON_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            -1850.328847999107, -2217.311570386829, -2992.92814085657,
            -4620.7839353297895, -2273.4729070927874, -2532.855314214248,
            -4513.7210813907695, -5834.929363875373, -4176.784171535893,
            -4586.555695275532, -6188.717478630532, -5041.310280965145,
            -8830.676376595093, -10267.697715325847, -10024.035433493409
    };

    /**
     * Fixed regression baseline for OU (dense selection matrix) precision gradient
     * under missing observations.
     */
    private static final double[] EXPECTED_OU_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            74.81412830364764, -258.8947558689557, 344.45582590843685,
            583.9296371742928, -136.52873168781966, -157.62714026311528,
            546.8452930351989, -211.5120616379187, -393.840964864784,
            -958.7566153189082, 1344.79636855238, 165.77240807171225,
            1106.0518906036552, 1542.6609607380626, 1495.4826838084582
    };

    /**
     * Fixed regression baseline for diagonal-OU precision gradient
     * under missing observations.
     */
    private static final double[] EXPECTED_DIAGONAL_OU_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            -47.4651521746416, 128.03418950082556, 56.24155717568088,
            409.4326733059245, 33.29452653571957, 179.79517412355858,
            114.92483761711222, 71.4085872025047, -147.28687327044906,
            -99.56318995900003, 155.69951217027022, 496.14703204294847,
            393.2262692791022, 971.3033044828563, 665.974550525364
    };

    /**
     * Fixed regression baseline for OU variance-parameterization precision gradient
     * in the non-missing setup.
     */
    private static final double[] EXPECTED_OU_VARIANCE_NON_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            530.5315456101581, 664.4705355229135, 482.8847773523566,
            723.4550261857936, 433.4566617734267, 559.8427844975721,
            544.2785158339794, 849.7834027385901, 213.17927963333398,
            578.5315342058839, 921.0664900153004, 318.5698519227596,
            778.5006123467441, 52.57294931842925, 252.52779037725682
    };

    /**
     * Fixed regression baseline for OU variance-parameterization precision gradient
     * under missing observations.
     */
    private static final double[] EXPECTED_OU_VARIANCE_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            -10.482437350865844, 27.0313040567109, -31.092711918647616,
            -102.38301687315386, 7.756758662822043, -19.244990977505136,
            -12.755678651396357, 6.465320565983597, 18.798919905625663,
            100.92325596924675, -10.980981100998513, -59.557876563478985,
            41.72305184685205, -112.92353191311265, -89.89097676178369
    };

    /**
     * Fixed regression baseline for diagonal-OU factor-model precision gradient.
     */
    private static final double[] EXPECTED_DIAGONAL_OU_FACTOR_PRECISION_CORRELATION_ANALYTIC = new double[]{
            2.474488004265763
    };

    /**
     * Fixed regression baseline for dense-OU factor-model precision gradient.
     */
    private static final double[] EXPECTED_OU_FACTOR_PRECISION_CORRELATION_ANALYTIC = new double[]{
            6.402724355632401
    };

    /**
     * Fixed regression baseline for OU single-opt precision gradient
     * in the non-missing setup.
     */
    private static final double[] EXPECTED_OU_SINGLE_OPT_NON_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            -1622.0894878807821, -3774.3995639703307, -2915.325514354333,
            -13243.604795626403, -242.0704089771807, -4528.018371460863,
            -4987.054769383916, -14784.68858069885, -1505.6872170230727,
            -12263.952702767445, -24565.416894482838, -14632.184456870398,
            -32444.84630412644, -26809.398304189403, -45821.2053023945
    };

    /**
     * Fixed regression baseline for OU single-opt precision gradient
     * under missing observations.
     */
    private static final double[] EXPECTED_OU_SINGLE_OPT_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            669.8653636679958, -1100.8297232129419, 784.0846529384476,
            -6976.889189236126, 2199.144046762984, -1726.5394150875409,
            376.45247000850304, -7868.141764035929, 5362.483381185795,
            -7500.650822752845, -17021.270930979, -6221.813664529498,
            -19291.08017009647, -9980.321177007034, -29538.032844942984
    };

    /**
     * Fixed regression baseline for OU single-opt variance-parameterization precision gradient
     * in the non-missing setup.
     */
    private static final double[] EXPECTED_OU_SINGLE_OPT_VARIANCE_NON_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            488.2009536008699, 529.1293899515387, 410.2140957373376,
            627.8428481181519, 256.24312490560914, 612.7994631877528,
            619.7566646566453, 980.5640659515681, 234.6436311070708,
            911.2907118296798, 1444.1979582617273, 518.018046192173,
            1490.3465751718616, 559.2972316195027, 994.870435968866
    };

    /**
     * Fixed regression baseline for OU single-opt variance-parameterization precision gradient
     * under missing observations.
     */
    private static final double[] EXPECTED_OU_SINGLE_OPT_VARIANCE_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            9.822423272921199, -51.20056560477846, -66.35071484503358,
            -94.73391738415705, 149.65039618532455, 26.218053989152192,
            30.579426006920194, 89.6728764395342, 62.66917648158248,
            592.0227731044948, 805.3361146791153, 367.1359231341601,
            877.196233440105, 347.6997808807747, 838.5736717555991
    };

    /**
     * Fixed regression baseline for OU same-opt precision gradient
     * in the non-missing setup.
     */
    private static final double[] EXPECTED_OU_SAME_OPT_NON_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            -2515.367588644446, -2740.863460433436, -3726.683908887034,
            -7176.048550100972, -3524.903065898143, -2962.212896489451,
            -6858.331156814085, -9253.318832137942, -5190.383325217147,
            -6420.707623229729, -10354.503767013442, -8143.210857861226,
            -15128.521778856662, -13223.841962065268, -16075.644017221393
    };

    /**
     * Fixed regression baseline for OU same-opt precision gradient
     * under missing observations.
     */
    private static final double[] EXPECTED_OU_SAME_OPT_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            229.5800641478644, 285.6877643650421, 68.75722262615182,
            725.8897536524848, 546.3369789564646, 533.9272750517844,
            -95.18639618951181, -176.7622309667853, 1941.019313517878,
            422.7234861644795, 588.0840075174815, -6.6301491810298785,
            1762.8801461771718, 3887.2080573436397, 1354.5970318455275
    };

    /**
     * Fixed regression baseline for OU same-opt variance-parameterization precision gradient
     * in the non-missing setup.
     */
    private static final double[] EXPECTED_OU_SAME_OPT_VARIANCE_NON_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            699.9141698447906, 587.4869488149476, 623.4527374939938,
            788.1233630982264, 259.23457959852533, 534.0293098237914,
            757.509150004868, 953.6993392389652, 70.62113257885778,
            581.8095930880203, 922.0797979875407, 177.25769593312324,
            918.405629245733, -58.58745042633665, 146.47550396233487
    };

    /**
     * Fixed regression baseline for OU same-opt variance-parameterization precision gradient
     * under missing observations.
     */
    private static final double[] EXPECTED_OU_SAME_OPT_VARIANCE_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            14.686058959410179, -65.89222309931523, 53.49007572178578,
            -93.72616405459894, 0.4720513295325593, -98.1877566010237,
            127.61094579774434, 54.46464167983924, -165.4105484010552,
            17.710039371790348, 14.518756321438346, -39.28967146696643,
            44.26157222992314, -255.45864804302022, -153.98661974845177
    };

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
        assertCorrelationPrecisionAnalyticRegressionWithMissing(
                diffusionModel,
                diffusionProcessDelegate,
                dataModel,
                rootPrior,
                precisionMatrix,
                EXPECTED_OU_NON_MISSING_PRECISION_CORRELATION_ANALYTIC,
                "OU precision correlation analytic baseline");
        System.out.println("\nTest OU gradient precision with missing.");
        assertCorrelationPrecisionAnalyticRegressionWithMissing(
                diffusionModel,
                diffusionProcessDelegate,
                dataModelMissing,
                rootPrior,
                precisionMatrix,
                EXPECTED_OU_MISSING_PRECISION_CORRELATION_ANALYTIC,
                "OU-with-missing precision correlation analytic baseline");

        // Wrt Variance
        DiffusionProcessDelegate diffusionProcessDelegateVariance
                = new OUDiffusionModelDelegate(treeModel, diffusionModelVar,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        System.out.println("\nTest OU gradient variance.");
        assertCorrelationPrecisionAnalyticRegressionWithMissing(
                diffusionModelVar,
                diffusionProcessDelegateVariance,
                dataModel,
                rootPrior,
                precisionMatrixInv,
                EXPECTED_OU_VARIANCE_NON_MISSING_PRECISION_CORRELATION_ANALYTIC,
                "OU-variance precision correlation analytic baseline");
        System.out.println("\nTest OU gradient variance with missing.");
        assertCorrelationPrecisionAnalyticRegressionWithMissing(
                diffusionModelVar,
                diffusionProcessDelegateVariance,
                dataModelMissing,
                rootPrior,
                precisionMatrixInv,
                EXPECTED_OU_VARIANCE_MISSING_PRECISION_CORRELATION_ANALYTIC,
                "OU-variance-with-missing precision correlation analytic baseline");

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
        assertCorrelationPrecisionAnalyticRegressionWithMissing(
                diffusionModelFactor,
                diffusionProcessDelegateFactor,
                dataModelFactor,
                rootPriorFactor,
                precisionMatrixFactor,
                EXPECTED_OU_FACTOR_PRECISION_CORRELATION_ANALYTIC,
                "OU factor precision correlation analytic baseline");

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
        assertCorrelationPrecisionAnalyticRegressionWithMissing(
                diffusionModel,
                diffusionProcessDelegateSingle,
                dataModel,
                rootPrior,
                precisionMatrix,
                EXPECTED_OU_SINGLE_OPT_NON_MISSING_PRECISION_CORRELATION_ANALYTIC,
                "OU single-opt precision correlation analytic baseline");
        System.out.println("\nTest OU single opt gradient precision with missing.");
        assertCorrelationPrecisionAnalyticRegressionWithMissing(
                diffusionModel,
                diffusionProcessDelegateSingle,
                dataModelMissing,
                rootPrior,
                precisionMatrix,
                EXPECTED_OU_SINGLE_OPT_MISSING_PRECISION_CORRELATION_ANALYTIC,
                "OU single-opt-with-missing precision correlation analytic baseline");

        // Wrt Variance
        DiffusionProcessDelegate diffusionProcessDelegateVarianceSingle
                = new OUDiffusionModelDelegate(treeModel, diffusionModelVar,
                optimalTraitsModelsSingle, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        System.out.println("\nTest OU single opt gradient variance.");
        assertCorrelationPrecisionAnalyticRegressionWithMissing(
                diffusionModelVar,
                diffusionProcessDelegateVarianceSingle,
                dataModel,
                rootPrior,
                precisionMatrixInv,
                EXPECTED_OU_SINGLE_OPT_VARIANCE_NON_MISSING_PRECISION_CORRELATION_ANALYTIC,
                "OU single-opt variance precision correlation analytic baseline");
        System.out.println("\nTest OU single opt gradient variance with missing.");
        assertCorrelationPrecisionAnalyticRegressionWithMissing(
                diffusionModelVar,
                diffusionProcessDelegateVarianceSingle,
                dataModelMissing,
                rootPrior,
                precisionMatrixInv,
                EXPECTED_OU_SINGLE_OPT_VARIANCE_MISSING_PRECISION_CORRELATION_ANALYTIC,
                "OU single-opt variance-with-missing precision correlation analytic baseline");

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
        assertCorrelationPrecisionAnalyticRegressionWithMissing(
                diffusionModel,
                diffusionProcessDelegateSame,
                dataModel,
                rootPrior,
                precisionMatrix,
                EXPECTED_OU_SAME_OPT_NON_MISSING_PRECISION_CORRELATION_ANALYTIC,
                "OU same-opt precision correlation analytic baseline");
        System.out.println("\nTest OU Same opt gradient precision with missing.");
        assertCorrelationPrecisionAnalyticRegressionWithMissing(
                diffusionModel,
                diffusionProcessDelegateSame,
                dataModelMissing,
                rootPrior,
                precisionMatrix,
                EXPECTED_OU_SAME_OPT_MISSING_PRECISION_CORRELATION_ANALYTIC,
                "OU same-opt-with-missing precision correlation analytic baseline");

        // Wrt Variance
        DiffusionProcessDelegate diffusionProcessDelegateVarianceSame
                = new OUDiffusionModelDelegate(treeModel, diffusionModelVar,
                optimalTraitsModelsSame, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        System.out.println("\nTest OU Same opt gradient variance.");
        assertCorrelationPrecisionAnalyticRegressionWithMissing(
                diffusionModelVar,
                diffusionProcessDelegateVarianceSame,
                dataModel,
                rootPrior,
                precisionMatrixInv,
                EXPECTED_OU_SAME_OPT_VARIANCE_NON_MISSING_PRECISION_CORRELATION_ANALYTIC,
                "OU same-opt variance precision correlation analytic baseline");
        System.out.println("\nTest OU Same opt gradient variance with missing.");
        assertCorrelationPrecisionAnalyticRegressionWithMissing(
                diffusionModelVar,
                diffusionProcessDelegateVarianceSame,
                dataModelMissing,
                rootPrior,
                precisionMatrixInv,
                EXPECTED_OU_SAME_OPT_VARIANCE_MISSING_PRECISION_CORRELATION_ANALYTIC,
                "OU same-opt variance-with-missing precision correlation analytic baseline");

        // Repeated Measures Model
        System.out.println("\nTest gradient precision repeated measures.");
        testGradient(diffusionModel, diffusionProcessDelegateSame, dataModelRepeatedMeasures, rootPrior, meanRoot, precisionMatrix, false, null, meanRoot, samplingPrecision);
        testGradient(diffusionModel, diffusionProcessDelegateSame, dataModelRepeatedMeasuresInv, rootPrior, meanRoot, precisionMatrix, false, null, meanRoot, samplingPrecisionInv);
    }

    public void testGradientOUSingleOptWithMissing_BranchBreakdownMatchesLegacyAccumulator() {

        final List<BranchRateModel> optimalTraitsModelsSingle = new ArrayList<BranchRateModel>();
        final CompoundParameter optParamSingle = new CompoundParameter("opt.breakdown");
        for (int i = 0; i < dimTrait; i++) {
            final Parameter rate = new Parameter.Default("opt.breakdown." + (i + 1), new double[]{2.2 * i + 3.1});
            optParamSingle.addParameter(rate);
            optimalTraitsModelsSingle.add(new StrictClockBranchRates(rate));
        }

        final MatrixParameter selectionStrength = buildDefaultOuSelectionStrengthMatrix();
        final OUDiffusionModelDelegate diffusionProcessDelegate =
                new OUDiffusionModelDelegate(treeModel, diffusionModel, optimalTraitsModelsSingle,
                        new MultivariateElasticModel(selectionStrength));

        final ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelMissing, rootPrior, rateTransformation, rateModel, true);
        final TreeDataLikelihood dataLikelihood = new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel);

        final ContinuousProcessParameterGradient traitGradient =
                new ContinuousProcessParameterGradient(
                        dimTrait,
                        treeModel,
                        likelihoodDelegate,
                        new ArrayList<>(Arrays.asList(DerivationParameter.WRT_CONSTANT_DRIFT)));

        final BranchSpecificGradient branchSpecificGradient =
                new BranchSpecificGradient("trait", dataLikelihood, likelihoodDelegate, traitGradient, optParamSingle);

        final double[] analyticGlobal = branchSpecificGradient.getGradientLogDensity();

        @SuppressWarnings("unchecked")
        final TreeTrait<List<BranchSufficientStatistics>> traitProvider =
                dataLikelihood.getTreeTrait(BranchConditionalDistributionDelegate.getName("trait"));
        assertNotNull("Branch conditional trait must be available", traitProvider);

        final double[] branchSummed = new double[optParamSingle.getDimension()];
        for (int i = 0; i < treeModel.getNodeCount(); ++i) {
            final NodeRef node = treeModel.getNode(i);
            final List<BranchSufficientStatistics> statisticsForNode = traitProvider.getTrait(treeModel, node);
            assertEquals("single-trait setup expected", 1, statisticsForNode.size());

            final BranchSufficientStatistics statistics = statisticsForNode.get(0);
            final double[] observedBranchGradient = traitGradient.getGradientForBranch(statistics, node);
            for (int d = 0; d < branchSummed.length; ++d) {
                final double observed = observedBranchGradient[d];
                assertTrue("Observed branch gradient must be finite at node " + node.getNumber() + ", dim " + d,
                        Double.isFinite(observed));
                branchSummed[d] += observed;
            }
        }

        assertEquals("global gradient length", analyticGlobal.length, branchSummed.length);
        for (int d = 0; d < analyticGlobal.length; ++d) {
            assertEquals("global branch sum dim " + d, branchSummed[d], analyticGlobal[d], 1e-10);
        }

        final double[] numericGlobal = numericalGradient(dataLikelihood, optParamSingle, 1e-6);
        assertEquals("global numeric length", analyticGlobal.length, numericGlobal.length);
        for (int d = 0; d < analyticGlobal.length; ++d) {
            assertEquals("global numeric dim " + d, analyticGlobal[d], numericGlobal[d], 2e-2);
        }
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
        assertCorrelationPrecisionAnalyticRegressionWithMissing(
                diffusionModel,
                diffusionProcessDelegate,
                dataModel,
                rootPrior,
                precisionMatrix,
                EXPECTED_DIAGONAL_OU_NON_MISSING_PRECISION_CORRELATION_ANALYTIC,
                "Diagonal OU precision correlation analytic baseline");
        System.out.println("\nTest Diagonal OU gradient precision with missing.");
        assertCorrelationPrecisionAnalyticRegressionWithMissing(
                diffusionModel,
                diffusionProcessDelegate,
                dataModelMissing,
                rootPrior,
                precisionMatrix,
                EXPECTED_DIAGONAL_OU_MISSING_PRECISION_CORRELATION_ANALYTIC,
                "Diagonal-OU-with-missing precision correlation analytic baseline");

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
        assertCorrelationPrecisionAnalyticRegressionWithMissing(
                diffusionModelFactor,
                diffusionProcessDelegateFactor,
                dataModelFactor,
                rootPriorFactor,
                precisionMatrixFactor,
                EXPECTED_DIAGONAL_OU_FACTOR_PRECISION_CORRELATION_ANALYTIC,
                "Diagonal-OU factor precision correlation analytic baseline");

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

    private void assertCorrelationPrecisionAnalyticRegressionWithMissing(
            final MultivariateDiffusionModel diffusionModel,
            final DiffusionProcessDelegate diffusionProcessDelegate,
            final ContinuousTraitPartialsProvider dataModel,
            final ConjugateRootTraitPrior rootPrior,
            final MatrixParameterInterface precision,
            final double[] expectedAnalytic,
            final String label) {

        final int dimLocal = rootPrior.getMean().length;
        final ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(
                treeModel, diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel, true);
        final TreeDataLikelihood dataLikelihood = new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel);

        final ProcessSimulationDelegate simulationDelegate =
                likelihoodDelegate.getPrecisionType() == PrecisionType.SCALAR
                        ? new ConditionalOnTipsRealizedDelegate(
                        "trait", treeModel, diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate)
                        : new MultivariateConditionalOnTipsRealizedDelegate(
                        "trait", treeModel, diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
        final TreeTraitProvider traitProvider = new ProcessSimulation(dataLikelihood, simulationDelegate);
        dataLikelihood.addTraits(traitProvider.getTreeTraits());

        final ProcessSimulationDelegate fullConditionalDelegate = new TipRealizedValuesViaFullConditionalDelegate(
                "trait", treeModel, diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
        dataLikelihood.addTraits(new ProcessSimulation(dataLikelihood, fullConditionalDelegate).getTreeTraits());

        final ContinuousProcessParameterGradient traitGradient = new ContinuousProcessParameterGradient(
                dimLocal, treeModel, likelihoodDelegate,
                new ArrayList<>(Arrays.asList(DerivationParameter.WRT_VARIANCE)));

        final BranchSpecificGradient branchSpecificGradient =
                new BranchSpecificGradient("trait", dataLikelihood, likelihoodDelegate, traitGradient, precision);

        final GradientWrtPrecisionProvider gPP =
                new GradientWrtPrecisionProvider.BranchSpecificGradientWrtPrecisionProvider(branchSpecificGradient);
        final CorrelationPrecisionGradient gradientProvider =
                new CorrelationPrecisionGradient(gPP, dataLikelihood, precision);

        final String report = gradientProvider.getReport();
        final double[] analytic = parseGradient(report, "analytic");

        assertEquals(label + " length", expectedAnalytic.length, analytic.length);
        for (int i = 0; i < analytic.length; i++) {
            assertEquals(label + " idx=" + i, expectedAnalytic[i], analytic[i], ANALYTIC_REGRESSION_TOLERANCE);
        }
    }

    private double[] parseVector(String s, String sep) {
        String[] vectorString = s.split(sep);
        double[] gradient = new double[vectorString.length];
        for (int i = 0; i < gradient.length; i++) {
            gradient[i] = Double.parseDouble(vectorString[i]);
        }
        return gradient;
    }

    private MatrixParameter buildDefaultOuSelectionStrengthMatrix() {
        final Parameter[] strengthOfSelectionParameters = new Parameter[6];
        strengthOfSelectionParameters[0] = new Parameter.Default(new double[]{10, 1.0, 0.0, 0.0, 0.0, 2.0});
        strengthOfSelectionParameters[1] = new Parameter.Default(new double[]{1.0, 20, 0.0, 0.0, 0.0, 0.0});
        strengthOfSelectionParameters[2] = new Parameter.Default(new double[]{0.0, 0.0, 30, 0.0, 0.0, 0.0});
        strengthOfSelectionParameters[3] = new Parameter.Default(new double[]{0.0, 0.0, 0.0, 40, 3.0, 0.0});
        strengthOfSelectionParameters[4] = new Parameter.Default(new double[]{0.0, 0.0, 0.0, 3.0, 50, 0.0});
        strengthOfSelectionParameters[5] = new Parameter.Default(new double[]{2.0, 0.0, 0.0, 0.0, 0.0, 60});
        return new MatrixParameter("strengthOfSelectionMatrix.breakdown", strengthOfSelectionParameters);
    }

    private double[] numericalGradient(final TreeDataLikelihood likelihood,
                                       final Parameter parameter,
                                       final double step) {
        final double[] values = parameter.getParameterValues().clone();
        final MultivariateFunction numeric = new MultivariateFunction() {
            @Override
            public double evaluate(final double[] argument) {
                for (int i = 0; i < argument.length; ++i) {
                    parameter.setParameterValue(i, argument[i]);
                }
                likelihood.makeDirty();
                return likelihood.getLogLikelihood();
            }

            @Override
            public int getNumArguments() {
                return parameter.getDimension();
            }

            @Override
            public double getLowerBound(final int n) {
                return Double.NEGATIVE_INFINITY;
            }

            @Override
            public double getUpperBound(final int n) {
                return Double.POSITIVE_INFINITY;
            }
        };

        final double[] gradient = new double[values.length];
        for (int i = 0; i < values.length; ++i) {
            final double x = values[i];

            values[i] = x + step;
            final double fxPlus = numeric.evaluate(values);

            values[i] = x - step;
            final double fxMinus = numeric.evaluate(values);

            gradient[i] = (fxPlus - fxMinus) / (2.0 * step);
            values[i] = x;
        }

        for (int i = 0; i < values.length; ++i) {
            parameter.setParameterValue(i, values[i]);
        }
        likelihood.makeDirty();
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
