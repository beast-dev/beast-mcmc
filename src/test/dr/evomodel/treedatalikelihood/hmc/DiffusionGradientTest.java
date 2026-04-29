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
            -2012.127768501196, -3787.7531988258534, -3426.2145058439846,
            -7036.5833457890185, -4018.097365417933, -4183.156426300409,
            -5373.785842470059, -8907.560205306832, -5465.187327209527,
            -8065.7570358879875, -10510.668199256786, -10767.514083273374,
            -14240.923489004213, -15551.560709452022, -16941.62119344948
    };

    /**
     * Fixed regression baseline for diagonal-OU precision gradient
     * in the non-missing setup.
     */
    private static final double[] EXPECTED_DIAGONAL_OU_NON_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            -1850.3288480697158, -2217.31157046912, -2992.928140870981,
            -4620.7839352104775, -2273.473136596309, -2532.8553143195177,
            -4513.721081434111, -5834.929363789856, -4176.784386493142,
            -4586.555695302157, -6188.717478538401, -5041.310481699882,
            -8830.676376803682, -10267.697899372732, -10024.03561546394
    };

    /**
     * Fixed regression baseline for OU (dense selection matrix) precision gradient
     * under missing observations.
     */
    private static final double[] EXPECTED_OU_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            74.81391862667019, -258.89488804396933, 344.4556704568904,
            583.929816421075, -136.52944184586673, -157.6271603462804,
            546.8452998103903, -211.5112380085236, -393.8423768503005,
            -958.7568155827505, 1344.7965381938636, 165.7719104987697,
            1106.0525381603493, 1542.6599949579872, 1495.4748103407157
    };

    /**
     * Fixed regression baseline for diagonal-OU precision gradient
     * under missing observations.
     */
    private static final double[] EXPECTED_DIAGONAL_OU_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            -47.46515223451834, 128.03418942876442, 56.24155716780447,
            409.43267342188796, 33.294340437731805, 179.79517402329742,
            114.92483757465985, 71.408587280647, -147.2870475020756,
            -99.56318998813293, 155.69951225274968, 496.1468693477042,
            393.2262690816413, 971.3031555915368, 665.9744017787768
    };

    /**
     * Fixed regression baseline for OU variance-parameterization precision gradient
     * in the non-missing setup.
     */
    private static final double[] EXPECTED_OU_VARIANCE_NON_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            530.5315456066558, 664.4705352396178, 482.884777661865,
            723.4550291853658, 433.4566581584673, 559.8427844806681,
            544.2785158674445, 849.7834030285774, 213.17927942793426,
            578.5315361196583, 921.0665038304086, 318.56983457911605,
            778.5005956944276, 52.57296725816232, 252.52796946879525
    };

    /**
     * Fixed regression baseline for OU variance-parameterization precision gradient
     * under missing observations.
     */
    private static final double[] EXPECTED_OU_VARIANCE_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            -10.482437358766541, 27.031303767529526, -31.092711685298088,
            -102.3830141970221, 7.756756272319861, -19.244990994556108,
            -12.755678599828784, 6.4653209313050635, 18.798919422942145,
            100.92325792379417, -10.980967111504567, -59.557894247718195,
            41.7230360138148, -112.92351852995164, -89.89081710783219
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
            -1622.0909953894848, -3774.4005042351605, -2915.324844558352,
            -13243.59486160547, -242.07861628795945, -4528.018322046966,
            -4987.054963341072, -14784.685041197823, -1505.708166022031,
            -12263.95270019819, -24565.412042249, -14632.200023242582,
            -32444.83566618452, -26809.397180483746, -45821.15008383
    };

    /**
     * Fixed regression baseline for OU single-opt precision gradient
     * under missing observations.
     */
    private static final double[] EXPECTED_OU_SINGLE_OPT_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            669.8638770197365, -1100.8306552773595, 784.0852816907243,
            -6976.879391568891, 2199.1359654675703, -1726.5393654131212,
            376.45227292644404, -7868.138207433156, 5362.462775707946,
            -7500.650824437087, -17021.266097223364, -6221.829029892011,
            -19291.069526308136, -9980.320261176941, -29537.979120247386
    };

    /**
     * Fixed regression baseline for OU single-opt variance-parameterization precision gradient
     * in the non-missing setup.
     */
    private static final double[] EXPECTED_OU_SINGLE_OPT_VARIANCE_NON_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            488.2009534145884, 529.1293890995275, 410.21408873969324,
            627.8428230716063, 256.24309920864675, 612.7994630874826,
            619.7566624852499, 980.5640578352967, 234.6436198004557,
            911.2907006321641, 1444.197916921717, 518.0179944509362,
            1490.3461872751895, 559.2968237432422, 994.8689815178469
    };

    /**
     * Fixed regression baseline for OU single-opt variance-parameterization precision gradient
     * under missing observations.
     */
    private static final double[] EXPECTED_OU_SINGLE_OPT_VARIANCE_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            9.822423090185577, -51.20056643291538, -66.35072171925201,
            -94.7339420541207, 149.65037050740114, 26.218053897919216,
            30.57942397605362, 89.67286880786354, 62.6691653967964,
            592.0227619029495, 805.3360734619064, 367.13587282948293,
            877.1958444174901, 347.69938004682206, 838.572238843495
    };

    /**
     * Fixed regression baseline for OU same-opt precision gradient
     * in the non-missing setup.
     */
    private static final double[] EXPECTED_OU_SAME_OPT_NON_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            -2515.367588122637, -2740.8634588469395, -3726.683959137922,
            -7176.047606005556, -3524.9032360205124, -2962.2128950558067,
            -6858.331150176706, -9253.318318179077, -5190.3838110969355,
            -6420.7076222197265, -10354.503260215482, -8143.211301398615,
            -15128.521771454009, -13223.842815471444, -16075.641229624873
    };

    /**
     * Fixed regression baseline for OU same-opt precision gradient
     * under missing observations.
     */
    private static final double[] EXPECTED_OU_SAME_OPT_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            229.58007124282503, 285.6877720624481, 68.75717471918074,
            725.8907127868663, 546.3368286858915, 533.9272768045805,
            -95.1863893076461, -176.7617190114118, 1941.0188679697405,
            422.72348587504905, 588.0845126120216, -6.6305481538464335,
            1762.880157891705, 3887.207211415569, 1354.6000105760058
    };

    /**
     * Fixed regression baseline for OU same-opt variance-parameterization precision gradient
     * in the non-missing setup.
     */
    private static final double[] EXPECTED_OU_SAME_OPT_VARIANCE_NON_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            699.9141699827549, 587.4869489563561, 623.4527383707058,
            788.1233616717072, 259.23457611356565, 534.0293097822279,
            757.5091496145517, 953.6993411861988, 70.62114109657973,
            581.8095923750493, 922.079800124396, 177.25770481531646,
            918.4056462176139, -58.58739460589427, 146.47540170262934
    };

    /**
     * Fixed regression baseline for OU same-opt variance-parameterization precision gradient
     * under missing observations.
     */
    private static final double[] EXPECTED_OU_SAME_OPT_VARIANCE_MISSING_PRECISION_CORRELATION_ANALYTIC = new double[]{
            14.686059098424375, -65.89222296256858, 53.49007660243335,
            -93.72616551502153, 0.47204795730720983, -98.18775664346707,
            127.61094540671614, 54.46464360790316, -165.4105398230568,
            17.710038665902243, 14.518758484366751, -39.28966286823404,
            44.261589256420386, -255.45859192204554, -153.98672397503432
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
