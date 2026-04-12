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
            -2012.1795817606735, -3787.784464817003, -3426.2410228063436,
            -7036.539636254487, -4018.2642532715577, -4183.164826311454,
            -5373.785052935776, -8907.373338327994, -5465.523619531839,
            -8065.809493616765, -10510.678288896994, -10767.593579233519,
            -14240.864060526541, -15551.692834024976, -16943.44740858466
    };

    /**
     * Fixed regression baseline for diagonal-OU precision gradient
     * in the non-missing setup.
     */
    private static final double[] EXPECTED_DIAGONAL_OU_NON_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            -1850.3288931038878, -2217.311621496853, -2992.9281549462553,
            -4620.783881716824, -2273.53683495885, -2532.8553783004713,
            -4513.721112107637, -5834.929328420172, -4176.843988270135,
            -4586.555713289437, -6188.71743675792, -5041.366183701355,
            -8830.676498694364, -10267.748337353914, -10024.08844420044
    };

    /**
     * Fixed regression baseline for OU (dense selection matrix) precision gradient
     * under missing observations.
     */
    private static final double[] EXPECTED_OU_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            74.77521730783411, -258.91629314069223, 344.4388601867604,
            583.9969383383511, -136.64730444531838, -157.63545585619872,
            546.8456126015939, -211.32570666481473, -394.0709230467711,
            -958.8069616087868, 1344.7908760966486, 165.7747365073456,
            1106.1182677540924, 1542.6046364767465, 1493.8996215111927
    };

    /**
     * Fixed regression baseline for diagonal-OU precision gradient
     * under missing observations.
     */
    private static final double[] EXPECTED_DIAGONAL_OU_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            -47.465188418725745, 128.03414643106794, 56.24154886603152,
            409.432725801331, 33.24895215706398, 179.7951140273837,
            114.92480858803546, 71.40861849850492, -147.3294583787569,
            -99.56320968699453, 155.69954742029697, 496.1072344404216,
            393.226156115573, 971.2675165197544, 665.9357688239339
    };

    /**
     * Fixed regression baseline for OU variance-parameterization precision gradient
     * in the non-missing setup.
     */
    private static final double[] EXPECTED_OU_VARIANCE_NON_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            530.5315456200768, 664.4704537337831, 482.8848423836275,
            723.4557301694928, 433.4557447139943, 559.8427774356218,
            544.2785218427101, 849.7834895737302, 213.1792852230237,
            578.5320756376082, 921.071265862378, 318.56485017663255,
            778.4974862893661, 52.576295233008835, 252.5680116116773
    };

    /**
     * Fixed regression baseline for OU variance-parameterization precision gradient
     * under missing observations.
     */
    private static final double[] EXPECTED_OU_VARIANCE_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            -10.482439891425514, 27.03121931652482, -31.092675471394024,
            -102.38243478839736, 7.7563223690237075, -19.244998091681246,
            -12.755661038445322, 6.465455959714896, 18.798764922407226,
            100.9238191913746, -10.976111237354552, -59.563041846252204,
            41.72026586613858, -112.92185362620714, -89.8579337605535
    };

    /**
     * Fixed regression baseline for diagonal-OU factor-model precision gradient.
     */
    private static final double[] EXPECTED_DIAGONAL_OU_FACTOR_PRECISION_CORRELATION_ANALYTIC = new double[]{
            0.6616157521699866
    };

    /**
     * Fixed regression baseline for dense-OU factor-model precision gradient.
     */
    private static final double[] EXPECTED_OU_FACTOR_PRECISION_CORRELATION_ANALYTIC = new double[]{
            3.03205106645377
    };

    /**
     * Fixed regression baseline for OU single-opt precision gradient
     * in the non-missing setup.
     */
    private static final double[] EXPECTED_OU_SINGLE_OPT_NON_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            -1622.3240880164374, -3774.5277535610044, -2915.1408633857354,
            -13241.288270813326, -243.45392361467202, -4528.01676306718,
            -4987.120755206228, -14783.757790987926, -1509.2680361975654,
            -12263.948100269296, -24564.155783096656, -14634.689280179566,
            -32442.105470057824, -26808.571271195233, -45809.428209962054
    };

    /**
     * Fixed regression baseline for OU single-opt precision gradient
     * under missing observations.
     */
    private static final double[] EXPECTED_OU_SINGLE_OPT_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            669.6364711366855, -1100.9569618131447, 784.2527293000892,
            -6974.625417352369, 2197.7998348390156, -1726.537441136292,
            376.3858864730021, -7867.2055697257365, 5359.011219299139,
            -7500.647424009136, -17020.01899837163, -6224.266791825359,
            -19288.34311338752, -9979.586636761518, -29526.81807680111
    };

    /**
     * Fixed regression baseline for OU single-opt variance-parameterization precision gradient
     * in the non-missing setup.
     */
    private static final double[] EXPECTED_OU_SINGLE_OPT_VARIANCE_NON_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            488.2008903260898, 529.1291344991621, 410.21253345447553,
            627.8373745202729, 256.2380932905658, 612.7994110362752,
            619.755784475312, 980.5608268180988, 234.6398068969133,
            911.2866196820083, 1444.1831942094825, 518.0026953963004,
            1490.2491663234334, 559.2109732348184, 994.5710817509905
    };

    /**
     * Fixed regression baseline for OU single-opt variance-parameterization precision gradient
     * under missing observations.
     */
    private static final double[] EXPECTED_OU_SINGLE_OPT_VARIANCE_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            9.82236191608851, -51.20080948510681, -66.35222791637457,
            -94.73924710659901, 149.64534580109594, 26.218008531001846,
            30.578631837219266, 89.6699281605342, 62.66547285195611,
            592.0186762042639, 805.3214176602125, 367.1212614217344,
            877.0982991962223, 347.616213347652, 838.2822366924668
    };

    /**
     * Fixed regression baseline for OU same-opt precision gradient
     * in the non-missing setup.
     */
    private static final double[] EXPECTED_OU_SAME_OPT_NON_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            -2515.3676145387863, -2740.863726705725, -3726.7013964646403,
            -7175.8313979072955, -3524.939767143063, -2962.213315239569,
            -6858.330602862178, -9253.201452095285, -5190.487892096534,
            -6420.707205401952, -10354.389395445512, -8143.310388596071,
            -15128.608286328985, -13224.061106818237, -16075.09649609822
    };

    /**
     * Fixed regression baseline for OU same-opt precision gradient
     * under missing observations.
     */
    private static final double[] EXPECTED_OU_SAME_OPT_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            229.58234283163756, 285.6897988362582, 68.74003449597637,
            726.112468053775, 546.3063068891853, 533.9271564348537,
            -95.18558461567034, -176.6449985885731, 1940.926846797614,
            422.7232226688895, 588.1975964697981, -6.71402971211225,
            1762.796394523283, 3886.9855397607043, 1355.2141488444602
    };

    /**
     * Fixed regression baseline for OU same-opt variance-parameterization precision gradient
     * in the non-missing setup.
     */
    private static final double[] EXPECTED_OU_SAME_OPT_VARIANCE_NON_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            699.9142175541317, 587.4869902140994, 623.4529542409521,
            788.1231292242356, 259.2341229583513, 534.0292836673011,
            757.5089020057927, 953.7001834779971, 70.62411125461149,
            581.8092553090733, 922.0805252622079, 177.26035681553705,
            918.4112027160407, -58.57277656238341, 146.45400570621774
    };

    /**
     * Fixed regression baseline for OU same-opt variance-parameterization precision gradient
     * under missing observations.
     */
    private static final double[] EXPECTED_OU_SAME_OPT_VARIANCE_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            14.686107223558674, -65.89218395241254, 53.490294273557765,
            -93.72641111598142, 0.47164996580445623, -98.18778338529107,
            127.61069765719785, 54.46547446798413, -165.40753795720485,
            17.709705476030226, 14.51949769566256, -39.28714667096465,
            44.26716890668358, -255.44382051043192, -154.00885692825605
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
