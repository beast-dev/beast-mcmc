package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodel.treedatalikelihood.preorder.MultivariateConditionalOnTipsRealizedDelegate;
import dr.evomodel.treedatalikelihood.preorder.ProcessSimulationDelegate;
import dr.evomodel.treelikelihood.utilities.TreeTraitLogger;
import dr.inference.model.DiagonalMatrix;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.Vector;

import java.util.ArrayList;
import java.util.List;

import static dr.evomodel.branchratemodel.ArbitraryBranchRates.make;


/**
 * @author Paul Bastide
 */

public class ContinuousDataLikelihoodDelegateTest extends ContinuousTraitTest {

    public ContinuousDataLikelihoodDelegateTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();
    }

    public void testLikelihoodBM() {
        System.out.println("\nTest Likelihood using vanilla BM:");

        // Diffusion
        DiffusionProcessDelegate diffusionProcessDelegate
                = new HomogeneousDiffusionModelDelegate(treeModel, diffusionModel);

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel, true);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihood = new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel);

        testLikelihood("likelihoodBM", dataLikelihood);

        // Conditional moments (preorder)
        testConditionalMoments(dataLikelihood, likelihoodDelegate);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        double[] expectedTraits = new double[]{-1.0, 2.0, 0.0, 0.45807521679597646, 2.6505355982097605, 3.4693334367360538, 0.5, 2.64206285585883, 5.5, 2.0, 5.0, -8.0, 11.0, 1.0, -1.5, 1.0, 2.5, 4.0};
        testConditionalSimulations(dataLikelihood, likelihoodDelegate, diffusionModel, dataModel, rootPrior, expectedTraits);

        // Fixed Root
        ContinuousDataLikelihoodDelegate likelihoodDelegateInf = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPriorInf, rateTransformation, rateModel, true);
        TreeDataLikelihood dataLikelihoodInf = new TreeDataLikelihood(likelihoodDelegateInf, treeModel, rateModel);
        testLikelihood("likelihoodBMInf", dataLikelihoodInf);
        testConditionalMoments(dataLikelihoodInf, likelihoodDelegateInf);

    }

    public void testLikelihoodDrift() {
        System.out.println("\nTest Likelihood using Drifted BM:");

        // Diffusion
        List<BranchRateModel> driftModels = new ArrayList<BranchRateModel>();
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{100.0})));
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{200.0})));
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.3", new double[]{-200.0})));
        DiffusionProcessDelegate diffusionProcessDelegate
                = new DriftDiffusionModelDelegate(treeModel, diffusionModel, driftModels);

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel, false);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihood = new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel);

        testLikelihood("likelihoodDrift", dataLikelihood);

        // Conditional moments (preorder)
        testConditionalMoments(dataLikelihood, likelihoodDelegate);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        double[] expectedTraits = new double[]{-1.0, 2.0, 0.0, 0.5457621072639138, 3.28662834718796, 3.2939596558001845, 0.5, 1.0742799493604265, 5.5, 2.0, 5.0, -8.0, 11.0, 1.0, -1.5, 1.0, 2.5, 4.0};
        testConditionalSimulations(dataLikelihood, likelihoodDelegate, diffusionModel, dataModel, rootPrior, expectedTraits);


        // Fixed Root
        ContinuousDataLikelihoodDelegate likelihoodDelegateInf = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPriorInf, rateTransformation, rateModel, true);
        TreeDataLikelihood dataLikelihoodInf = new TreeDataLikelihood(likelihoodDelegateInf, treeModel, rateModel);
        testLikelihood("likelihoodDriftInf", dataLikelihoodInf);
        testConditionalMoments(dataLikelihoodInf, likelihoodDelegateInf);

    }

    public void testLikelihoodDriftRelaxed() {
        System.out.println("\nTest Likelihood using Drifted relaxed BM:");

        // Diffusion
        List<BranchRateModel> driftModels = new ArrayList<BranchRateModel>();
        ArbitraryBranchRates.BranchRateTransform transform = make(false, false, false);
        driftModels.add(new ArbitraryBranchRates(treeModel,
                new Parameter.Default("rate.1", new double[]{0, 100, 200, 300, 400, 500, 600, 700, 800, 900}),
                transform, false));
        driftModels.add(new ArbitraryBranchRates(treeModel,
                new Parameter.Default("rate.2", new double[]{0, -100, 200, -300, 400, -500, 600, -700, 800, -900}),
                transform, false));
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.3", new double[]{-2.0})));
        DiffusionProcessDelegate diffusionProcessDelegate
                = new DriftDiffusionModelDelegate(treeModel, diffusionModel, driftModels);

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel, false);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihood = new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel);

        testLikelihood("likelihoodDriftRelaxed", dataLikelihood);

        // Conditional moments (preorder)
        testConditionalMoments(dataLikelihood, likelihoodDelegate);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        double[] expectedTraits = new double[]{-1.0, 2.0, 0.0, 2.843948876154644, 10.866053719140933, 3.467579698926694, 0.5, 12.000214659757933, 5.5, 2.0, 5.0, -8.0, 11.0, 1.0, -1.5, 1.0, 2.5, 4.0};
        testConditionalSimulations(dataLikelihood, likelihoodDelegate, diffusionModel, dataModel, rootPrior, expectedTraits);


        // Fixed Root
        ContinuousDataLikelihoodDelegate likelihoodDelegateInf = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPriorInf, rateTransformation, rateModel, true);
        TreeDataLikelihood dataLikelihoodInf = new TreeDataLikelihood(likelihoodDelegateInf, treeModel, rateModel);
        testLikelihood("likelihoodDriftRelaxedInf", dataLikelihoodInf);
        testConditionalMoments(dataLikelihoodInf, likelihoodDelegateInf);
    }

    public void testLikelihoodDiagonalOU() {
        System.out.println("\nTest Likelihood using Diagonal OU:");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{1.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{2.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.3", new double[]{-2.0})));

        DiagonalMatrix strengthOfSelectionMatrixParam
                = new DiagonalMatrix(new Parameter.Default(new double[]{0.1, 100.0, 50.0}));

        DiffusionProcessDelegate diffusionProcessDelegate
                = new OUDiffusionModelDelegate(treeModel, diffusionModel,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel, false);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihood = new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel);

        testLikelihood("likelihoodDiagonalOU", dataLikelihood);

        // Conditional moments (preorder)
        testConditionalMoments(dataLikelihood, likelihoodDelegate);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        double[] expectedTraits = new double[]{-1.0, 2.0, 0.0, 1.0369756434332398, 2.0462155265038557, 0.6054236464363175, 0.5, 2.0829935706195615, 5.5, 2.0, 5.0, -8.0, 11.0, 1.0, -1.5, 1.0, 2.5, 4.0};
        testConditionalSimulations(dataLikelihood, likelihoodDelegate, diffusionModel, dataModel, rootPrior, expectedTraits);


        // Fixed Root
        ContinuousDataLikelihoodDelegate likelihoodDelegateInf = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPriorInf, rateTransformation, rateModel, true);
        TreeDataLikelihood dataLikelihoodInf = new TreeDataLikelihood(likelihoodDelegateInf, treeModel, rateModel);
        testLikelihood("likelihoodDiagonalOUInf", dataLikelihoodInf);
        testConditionalMoments(dataLikelihoodInf, likelihoodDelegateInf);
    }

    public void testLikelihoodDiagonalOURelaxed() {
        System.out.println("\nTest Likelihood using Diagonal OU Relaxed:");

        // Diffusion

        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        ArbitraryBranchRates.BranchRateTransform transform = make(false, false, false);
        optimalTraitsModels.add(new ArbitraryBranchRates(treeModel,
                new Parameter.Default("rate.1", new double[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}),
                transform, false));
        optimalTraitsModels.add(new ArbitraryBranchRates(treeModel,
                new Parameter.Default("rate.2", new double[]{0, -1, 2, -3, 4, -5, 6, -7, 8, -9}),
                transform, false));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.3", new double[]{-2.0})));

        DiagonalMatrix strengthOfSelectionMatrixParam
                = new DiagonalMatrix(new Parameter.Default(new double[]{1.0, 100.0, 100.0}));

        DiffusionProcessDelegate diffusionProcessDelegate
                = new OUDiffusionModelDelegate(treeModel, diffusionModel,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel, false);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihood = new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel);

        testLikelihood("likelihoodDiagonalOURelaxed", dataLikelihood);

        // Conditional moments (preorder)
        testConditionalMoments(dataLikelihood, likelihoodDelegate);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        double[] expectedTraits = new double[]{-1.0, 2.0, 0.0, 1.811936920575975, 0.6645253568041535, -1.0807136572611238, 0.5, 3.8623525502275142, 5.5, 2.0, 5.0, -8.0, 11.0, 1.0, -1.5, 1.0, 2.5, 4.0};
        testConditionalSimulations(dataLikelihood, likelihoodDelegate, diffusionModel, dataModel, rootPrior, expectedTraits);

        // Fixed Root
        ContinuousDataLikelihoodDelegate likelihoodDelegateInf = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPriorInf, rateTransformation, rateModel, true);
        TreeDataLikelihood dataLikelihoodInf = new TreeDataLikelihood(likelihoodDelegateInf, treeModel, rateModel);
        testLikelihood("likelihoodDiagonalOURelaxedInf", dataLikelihoodInf);
        testConditionalMoments(dataLikelihoodInf, likelihoodDelegateInf);
    }

    public void testLikelihoodDiagonalOUBM() {
        System.out.println("\nTest Likelihood using Diagonal OU / BM:");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{1.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{2.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.3", new double[]{-2.0})));

        DiagonalMatrix strengthOfSelectionMatrixParam
                = new DiagonalMatrix(new Parameter.Default(new double[]{0.0, 0.000001, 50.0}));

        DiffusionProcessDelegate diffusionProcessDelegate
                = new OUDiffusionModelDelegate(treeModel, diffusionModel,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel, false);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihood = new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel);

        testLikelihood("likelihoodDiagonalOUBM", dataLikelihood);

        // Conditional moments (preorder)
        testConditionalMoments(dataLikelihood, likelihoodDelegate);

        // Fixed Root
        ContinuousDataLikelihoodDelegate likelihoodDelegateInf = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPriorInf, rateTransformation, rateModel, true);
        TreeDataLikelihood dataLikelihoodInf = new TreeDataLikelihood(likelihoodDelegateInf, treeModel, rateModel);
        testLikelihood("likelihoodDiagonalOUBMInf", dataLikelihoodInf);
        testConditionalMoments(dataLikelihoodInf, likelihoodDelegateInf);
    }

    public void testLikelihoodDiagonalOUBMInd() {
        System.out.println("\nTest Likelihood using Diagonal OU / BM:");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{1.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{-3.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.3", new double[]{-2.0})));

        DiagonalMatrix strengthOfSelectionMatrixParamOUBM
                = new DiagonalMatrix(new Parameter.Default(new double[]{0.0, 0.0, 50.0}));
        DiagonalMatrix strengthOfSelectionMatrixParamOU
                = new DiagonalMatrix(new Parameter.Default(new double[]{10.0, 20.0, 50.0}));

        DiagonalMatrix diffusionPrecisionMatrixParameter
                = new DiagonalMatrix(new Parameter.Default(new double[]{1.0, 2.0, 3.0}));
        MultivariateDiffusionModel diffusionModel = new MultivariateDiffusionModel(diffusionPrecisionMatrixParameter);

        DiffusionProcessDelegate diffusionProcessDelegateOUBM
                = new OUDiffusionModelDelegate(treeModel, diffusionModel, optimalTraitsModels,
                new MultivariateElasticModel(strengthOfSelectionMatrixParamOUBM));

        DiffusionProcessDelegate diffusionProcessDelegateOU
                = new OUDiffusionModelDelegate(treeModel, diffusionModel, optimalTraitsModels,
                new MultivariateElasticModel(strengthOfSelectionMatrixParamOU));

        DiffusionProcessDelegate diffusionProcessDelegateBM
                = new HomogeneousDiffusionModelDelegate(treeModel, diffusionModel);

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateOUBM = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegateOUBM, dataModel, rootPriorInf, rateTransformation, rateModel, false);
        ContinuousDataLikelihoodDelegate likelihoodDelegateOU = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegateOU, dataModel, rootPriorInf, rateTransformation, rateModel, false);
        ContinuousDataLikelihoodDelegate likelihoodDelegateBM = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegateBM, dataModel, rootPriorInf, rateTransformation, rateModel, false);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodOUBM = new TreeDataLikelihood(likelihoodDelegateOUBM, treeModel, rateModel);
        TreeDataLikelihood dataLikelihoodOU = new TreeDataLikelihood(likelihoodDelegateOU, treeModel, rateModel);
        TreeDataLikelihood dataLikelihoodBM = new TreeDataLikelihood(likelihoodDelegateBM, treeModel, rateModel);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        double[] traitsOUBM = getConditionalSimulations(dataLikelihoodOUBM, likelihoodDelegateOUBM, diffusionModel, dataModel, rootPriorInf, treeModel, rateTransformation);
        System.err.println(new Vector(traitsOUBM));
        MathUtils.setSeed(17890826);
        double[] traitsOU = getConditionalSimulations(dataLikelihoodOU, likelihoodDelegateOU, diffusionModel, dataModel, rootPriorInf, treeModel, rateTransformation);
        System.err.println(new Vector(traitsOU));
        MathUtils.setSeed(17890826);
        double[] traitsBM = getConditionalSimulations(dataLikelihoodBM, likelihoodDelegateBM, diffusionModel, dataModel, rootPriorInf, treeModel, rateTransformation);
        System.err.println(new Vector(traitsBM));

        // Check that missing dimensions with the same process have the same values
        assertEquals(format.format(traitsBM[3]), format.format(traitsOUBM[3]));
        assertEquals(format.format(traitsBM[4]), format.format(traitsOUBM[4]));
        assertEquals(format.format(traitsBM[7]), format.format(traitsOUBM[7]));
        assertEquals(format.format(traitsOU[5]), format.format(traitsOUBM[5]));

    }


    public void testLikelihoodFullOU() {
        System.out.println("\nTest Likelihood using Full OU:");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{1.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{2.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.3", new double[]{-2.0})));

        Parameter[] strengthOfSelectionParameters = new Parameter[3];
        strengthOfSelectionParameters[0] = new Parameter.Default(new double[]{0.5, 0.2, 0.0});
        strengthOfSelectionParameters[1] = new Parameter.Default(new double[]{0.2, 100.0, 0.1});
        strengthOfSelectionParameters[2] = new Parameter.Default(new double[]{0.0, 0.1, 50.5});
        MatrixParameter strengthOfSelectionMatrixParam
                = new MatrixParameter("strengthOfSelectionMatrix", strengthOfSelectionParameters);

        DiffusionProcessDelegate diffusionProcessDelegate
                = new OUDiffusionModelDelegate(treeModel, diffusionModel,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel, false);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihood = new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel);

        testLikelihood("likelihoodFullOU", dataLikelihood);

        // Conditional moments (preorder)
        testConditionalMoments(dataLikelihood, likelihoodDelegate);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        double[] expectedTraits = new double[]{-1.0, 2.0, 0.0, 1.0428617934711057, 2.0411050451484463, 0.5794607846668219, 0.5, 2.07249828895442, 5.5, 2.0, 5.0, -8.0, 11.0, 1.0, -1.5, 1.0, 2.5, 4.0};
        testConditionalSimulations(dataLikelihood, likelihoodDelegate, diffusionModel, dataModel, rootPrior, expectedTraits);

        // Fixed Root
        ContinuousDataLikelihoodDelegate likelihoodDelegateInf = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPriorInf, rateTransformation, rateModel, true);
        TreeDataLikelihood dataLikelihoodInf = new TreeDataLikelihood(likelihoodDelegateInf, treeModel, rateModel);
        testLikelihood("likelihoodFullOUInf", dataLikelihoodInf);
        testConditionalMoments(dataLikelihoodInf, likelihoodDelegateInf);
    }

    public void testLikelihoodFullOURelaxed() {
        System.out.println("\nTest Likelihood using Full OU Relaxed:");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        ArbitraryBranchRates.BranchRateTransform transform = make(false, false, false);
        optimalTraitsModels.add(new ArbitraryBranchRates(treeModel,
                new Parameter.Default("rate.1", new double[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}),
                transform, false));
        optimalTraitsModels.add(new ArbitraryBranchRates(treeModel,
                new Parameter.Default("rate.2", new double[]{0, -1, 2, -3, 4, -5, 6, -7, 8, -9}),
                transform, false));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.3", new double[]{-2.0})));

        Parameter[] strengthOfSelectionParameters = new Parameter[3];
        strengthOfSelectionParameters[0] = new Parameter.Default(new double[]{0.5, 0.2, 0.0});
        strengthOfSelectionParameters[1] = new Parameter.Default(new double[]{0.2, 10.5, 0.1});
        strengthOfSelectionParameters[2] = new Parameter.Default(new double[]{0.0, 0.1, 100.0});
        MatrixParameter strengthOfSelectionMatrixParam
                = new MatrixParameter("strengthOfSelectionMatrix", strengthOfSelectionParameters);

        DiffusionProcessDelegate diffusionProcessDelegate
                = new OUDiffusionModelDelegate(treeModel, diffusionModel,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel, false);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihood = new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel);

        testLikelihood("likelihoodFullOURelaxed", dataLikelihood);

        // Conditional moments (preorder)
        testConditionalMoments(dataLikelihood, likelihoodDelegate);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        double[] expectedTraits = new double[]{-1.0, 2.0, 0.0, 1.635010498015842, 2.864821373865659, -1.0853111798010757, 0.5, 3.3661883786009166, 5.5, 2.0, 5.0, -8.0, 11.0, 1.0, -1.5, 1.0, 2.5, 4.0};
        testConditionalSimulations(dataLikelihood, likelihoodDelegate, diffusionModel, dataModel, rootPrior, expectedTraits);

        // Fixed Root
        ContinuousDataLikelihoodDelegate likelihoodDelegateInf = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPriorInf, rateTransformation, rateModel, true);
        TreeDataLikelihood dataLikelihoodInf = new TreeDataLikelihood(likelihoodDelegateInf, treeModel, rateModel);
        testLikelihood("likelihoodFullOURelaxedInf", dataLikelihoodInf);
        testConditionalMoments(dataLikelihoodInf, likelihoodDelegateInf);
    }

    public void testLikelihoodFullAndDiagonalOU() {
        System.out.println("\nTest Likelihood comparing Full and Diagonal OU:");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        ArbitraryBranchRates.BranchRateTransform transform = make(false, false, false);
        optimalTraitsModels.add(new ArbitraryBranchRates(treeModel,
                new Parameter.Default("rate.1", new double[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}),
                transform, false));
        optimalTraitsModels.add(new ArbitraryBranchRates(treeModel,
                new Parameter.Default("rate.2", new double[]{0, -1, 2, -3, 4, -5, 6, -7, 8, -9}),
                transform, false));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.3", new double[]{-2.0})));

        Parameter[] strengthOfSelectionParameters = new Parameter[3];
        strengthOfSelectionParameters[0] = new Parameter.Default(new double[]{0.5, 0.0, 0.0});
        strengthOfSelectionParameters[1] = new Parameter.Default(new double[]{0.0, 10.5, 0.0});
        strengthOfSelectionParameters[2] = new Parameter.Default(new double[]{0.0, 0.0, 100.0});
        MatrixParameter strengthOfSelectionMatrixParam
                = new MatrixParameter("strengthOfSelectionMatrix", strengthOfSelectionParameters);

        DiagonalMatrix strengthOfSelectionMatrixParamDiagonal
                = new DiagonalMatrix(new Parameter.Default(new double[]{0.5, 10.5, 100.0}));

        DiffusionProcessDelegate diffusionProcessDelegate
                = new OUDiffusionModelDelegate(treeModel, diffusionModel,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        DiffusionProcessDelegate diffusionProcessDelegateDiagonal
                = new OUDiffusionModelDelegate(treeModel, diffusionModel,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParamDiagonal));

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel, false);

        ContinuousDataLikelihoodDelegate likelihoodDelegateDiagonal = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegateDiagonal, dataModel, rootPrior, rateTransformation, rateModel, false);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihood = new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel);

        TreeDataLikelihood dataLikelihoodDiagonal = new TreeDataLikelihood(likelihoodDelegateDiagonal, treeModel, rateModel);

        assertEquals("likelihoodFullDiagonalOU",
                format.format(dataLikelihood.getLogLikelihood()),
                format.format(dataLikelihoodDiagonal.getLogLikelihood()));
    }

    public void testLikelihoodIBM() {
        System.out.println("\nTest Likelihood using IBM:");

        // Diffusion
        List<BranchRateModel> driftModels = new ArrayList<BranchRateModel>();
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{100.0})));
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{200.0})));
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.3", new double[]{-200.0})));

        DiffusionProcessDelegate diffusionProcessDelegate
                = new IntegratedBMDiffusionModelDelegate(treeModel, diffusionModel,
                driftModels);

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, true, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelIntegrated, rootPriorIntegrated, rateTransformation, rateModel, false);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihood = new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel);

        String s = dataLikelihood.getReport();
        int indLikBeg = s.indexOf("logDatumLikelihood:") + 20;
        int indLikEnd = s.indexOf("\n", indLikBeg);
        char[] logDatumLikelihoodChar = new char[indLikEnd - indLikBeg + 1];
        s.getChars(indLikBeg, indLikEnd, logDatumLikelihoodChar, 0);
        double logDatumLikelihood = Double.parseDouble(String.valueOf(logDatumLikelihoodChar));

        assertEquals("likelihoodIBM",
                format.format(logDatumLikelihood),
                format.format(dataLikelihood.getLogLikelihood()));

        System.out.println("likelihoodIBM: " + format.format(logDatumLikelihood));

        // Conditional moments (preorder)
//        new TreeTipGradient("" +
//                "trait", dataLikelihood, likelihoodDelegate, null);
//        TreeTraitLogger treeTraitLogger = new TreeTraitLogger(treeModel,
//                new TreeTrait[]{dataLikelihood.getTreeTrait("fcd.trait")},
//                TreeTraitLogger.NodeRestriction.EXTERNAL);
//
//        String moments = treeTraitLogger.getReport();
//        double[] partials = parseVector(moments, "\t");
//        testCMeans(s, "cMean ", partials);
//        testCVariances(s, "cVar ", partials);
    }

    public void testLikelihoodFullIOU() {
        System.out.println("\nTest Likelihood using Full IOU:");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{1.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{2.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.3", new double[]{-2.0})));

        Parameter[] strengthOfSelectionParameters = new Parameter[3];
        strengthOfSelectionParameters[0] = new Parameter.Default(new double[]{0.5, 0.5, 0.0});
        strengthOfSelectionParameters[1] = new Parameter.Default(new double[]{0.2, 5, 0.1});
        strengthOfSelectionParameters[2] = new Parameter.Default(new double[]{0.0, 1.0, 10.0});
        MatrixParameter strengthOfSelectionMatrixParam
                = new MatrixParameter("strengthOfSelectionMatrix", strengthOfSelectionParameters);

        DiffusionProcessDelegate diffusionProcessDelegate
                = new IntegratedOUDiffusionModelDelegate(treeModel, diffusionModel,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, true, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelIntegrated, rootPriorIntegrated, rateTransformation, rateModel, false);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihood = new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel);

        String s = dataLikelihood.getReport();
        int indLikBeg = s.indexOf("logDatumLikelihood:") + 20;
        int indLikEnd = s.indexOf("\n", indLikBeg);
        char[] logDatumLikelihoodChar = new char[indLikEnd - indLikBeg + 1];
        s.getChars(indLikBeg, indLikEnd, logDatumLikelihoodChar, 0);
        double logDatumLikelihood = Double.parseDouble(String.valueOf(logDatumLikelihoodChar));

        assertEquals("likelihoodIOU",
                format.format(logDatumLikelihood),
                format.format(dataLikelihood.getLogLikelihood()));

        System.out.println("likelihoodIOU: " + format.format(logDatumLikelihood));

        // Conditional moments (preorder)
//        testConditionalMoments(dataLikelihood, likelihoodDelegate);
    }

    public void testLikelihoodFullNonSymmetricOU() {
        System.out.println("\nTest Likelihood using Full Non symmetric OU:");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{1.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{2.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.3", new double[]{-2.0})));

        Parameter[] strengthOfSelectionParameters = new Parameter[3];
        strengthOfSelectionParameters[0] = new Parameter.Default(new double[]{0.5, 0.0, 0.0});
        strengthOfSelectionParameters[1] = new Parameter.Default(new double[]{0.2, 100.0, 0.1});
        strengthOfSelectionParameters[2] = new Parameter.Default(new double[]{10.0, 0.1, 50.5});
        MatrixParameter strengthOfSelectionMatrixParam
                = new MatrixParameter("strengthOfSelectionMatrix", strengthOfSelectionParameters);

        DiffusionProcessDelegate diffusionProcessDelegate
                = new OUDiffusionModelDelegate(treeModel, diffusionModel,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel, false);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihood = new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel);

        testLikelihood("likelihoodFullNonSymmetricOU", dataLikelihood);

        // Conditional moments (preorder)
        testConditionalMoments(dataLikelihood, likelihoodDelegate);

        // Fixed Root
        ContinuousDataLikelihoodDelegate likelihoodDelegateInf = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPriorInf, rateTransformation, rateModel, true);
        TreeDataLikelihood dataLikelihoodInf = new TreeDataLikelihood(likelihoodDelegateInf, treeModel, rateModel);
        testLikelihood("likelihoodFullNonSymmetricOUInf", dataLikelihoodInf);
        testConditionalMoments(dataLikelihoodInf, likelihoodDelegateInf);
    }

    public void testLikelihoodFullOUNonSymmetricRelaxed() {
        System.out.println("\nTest Likelihood using Full Non symmetric OU Relaxed:");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        ArbitraryBranchRates.BranchRateTransform transform = make(false, false, false);
        optimalTraitsModels.add(new ArbitraryBranchRates(treeModel,
                new Parameter.Default("rate.1", new double[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}),
                transform, false));
        optimalTraitsModels.add(new ArbitraryBranchRates(treeModel,
                new Parameter.Default("rate.2", new double[]{0, -1, 2, -3, 4, -5, 6, -7, 8, -9}),
                transform, false));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.3", new double[]{-2.0})));

        Parameter[] strengthOfSelectionParameters = new Parameter[3];
        strengthOfSelectionParameters[0] = new Parameter.Default(new double[]{0.5, 0.0, 0.0});
        strengthOfSelectionParameters[1] = new Parameter.Default(new double[]{0.2, 100.0, 0.1});
        strengthOfSelectionParameters[2] = new Parameter.Default(new double[]{10.0, 0.1, 50.5});
        MatrixParameter strengthOfSelectionMatrixParam
                = new MatrixParameter("strengthOfSelectionMatrix", strengthOfSelectionParameters);

        DiffusionProcessDelegate diffusionProcessDelegate
                = new OUDiffusionModelDelegate(treeModel, diffusionModel,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel, false);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihood = new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel);

        testLikelihood("likelihoodFullNonSymmetricOURelaxed", dataLikelihood);

        // Conditional moments (preorder)
        testConditionalMoments(dataLikelihood, likelihoodDelegate);

        // Fixed Root
        ContinuousDataLikelihoodDelegate likelihoodDelegateInf = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPriorInf, rateTransformation, rateModel, true);
        TreeDataLikelihood dataLikelihoodInf = new TreeDataLikelihood(likelihoodDelegateInf, treeModel, rateModel);
        testLikelihood("likelihoodFullNonSymmetricOURelaxedInf", dataLikelihoodInf);
        testConditionalMoments(dataLikelihoodInf, likelihoodDelegateInf);
    }

    //// Factor Model //// *********************************************************************************************

    public void testLikelihoodBMFactor() {
        System.out.println("\nTest Likelihood using vanilla BM and factor:");

        // Diffusion
        DiffusionProcessDelegate diffusionProcessDelegate
                = new HomogeneousDiffusionModelDelegate(treeModel, diffusionModelFactor);

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateFactors = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelFactor, rootPriorFactor,
                rateTransformation, rateModel, true);

        dataModelFactor.setLikelihoodDelegate(likelihoodDelegateFactors);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodFactors
                = new TreeDataLikelihood(likelihoodDelegateFactors, treeModel, rateModel);

        testLikelihood("likelihoodBMFactor", dataModelFactor, dataLikelihoodFactors);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        double[] expectedTraits = new double[]{0.6002879987080073, 1.3630884580519484, 0.5250449300511655, 1.4853676908300644, 0.6673202215955497, 1.399820047380221, 1.0853554355129353, 1.6054879123935393, 0.4495494080256063, 1.4427296475118248, 0.8750789069500045, 1.8099596179292183};
        testConditionalSimulations(dataLikelihoodFactors, likelihoodDelegateFactors, diffusionModelFactor, dataModelFactor, rootPrior, expectedTraits);

    }

    public void testLikelihoodDriftFactor() {
        System.out.println("\nTest Likelihood using drifted BM and factor:");

        // Diffusion
        List<BranchRateModel> driftModels = new ArrayList<BranchRateModel>();
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{0.0})));
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{-40.0})));
        DiffusionProcessDelegate diffusionProcessDelegate
                = new DriftDiffusionModelDelegate(treeModel, diffusionModelFactor, driftModels);

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateFactors = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelFactor, rootPriorFactor,
                rateTransformation, rateModel, false);

        dataModelFactor.setLikelihoodDelegate(likelihoodDelegateFactors);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodFactors
                = new TreeDataLikelihood(likelihoodDelegateFactors, treeModel, rateModel);

        testLikelihood("likelihoodDriftFactor", dataModelFactor, dataLikelihoodFactors);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        double[] expectedTraits = new double[]{1.5058510863259034, -2.344107747791032, 1.415239714927795, -2.225937980916329, 1.5639840062954773, -2.3082612693286513, 1.9875205911751028, -2.1049011248405525, 1.3355460225282372, -2.2848471441564056, 1.742347318026791, -1.940903337116235};
        testConditionalSimulations(dataLikelihoodFactors, likelihoodDelegateFactors, diffusionModelFactor, dataModelFactor, rootPriorFactor, expectedTraits);

    }

    public void testLikelihoodDriftRelaxedFactor() {
        System.out.println("\nTest Likelihood using drifted Relaxed BM and factor:");

        // Diffusion
        List<BranchRateModel> driftModels = new ArrayList<BranchRateModel>();
        ArbitraryBranchRates.BranchRateTransform transform = make(false, false, false);
        driftModels.add(new ArbitraryBranchRates(treeModel,
                new Parameter.Default("rate.1", new double[]{0, 10, 20, 30, 40, 40, 30, 20, 10, 0}),
                transform, false));
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{0.0})));
        DiffusionProcessDelegate diffusionProcessDelegate
                = new DriftDiffusionModelDelegate(treeModel, diffusionModelFactor, driftModels);

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateFactors = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelFactor, rootPriorFactor,
                rateTransformation, rateModel, false);

        dataModelFactor.setLikelihoodDelegate(likelihoodDelegateFactors);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodFactors
                = new TreeDataLikelihood(likelihoodDelegateFactors, treeModel, rateModel);

        testLikelihood("likelihoodDriftRelaxedFactor", dataModelFactor, dataLikelihoodFactors);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        double[] expectedTraits = new double[]{0.21992781609528125, 1.271388273711557, 0.40761548539751596, 1.3682648770877144, 0.6599021787120436, 1.2830636141108613, 1.1488658943588324, 1.472103688153391, 0.8971632986744889, 1.20748933414854, 1.603739823726808, 1.4761482401796842};
        testConditionalSimulations(dataLikelihoodFactors, likelihoodDelegateFactors, diffusionModelFactor, dataModelFactor, rootPriorFactor, expectedTraits);

    }

    public void testLikelihoodDiagonalOUFactor() {
        System.out.println("\nTest Likelihood using diagonal OU and factor:");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{1.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{-1.5})));

        DiagonalMatrix strengthOfSelectionMatrixParam = new DiagonalMatrix(new Parameter.Default(new double[]{0.5, 50.0}));

        DiffusionProcessDelegate diffusionProcessDelegate
                = new OUDiffusionModelDelegate(treeModel, diffusionModelFactor,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateFactors = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelFactor, rootPriorFactor,
                rateTransformation, rateModel, false);

        dataModelFactor.setLikelihoodDelegate(likelihoodDelegateFactors);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodFactors
                = new TreeDataLikelihood(likelihoodDelegateFactors, treeModel, rateModel);

        testLikelihood("likelihoodDiagonalOUFactor", dataModelFactor, dataLikelihoodFactors);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        double[] expectedTraits = new double[]{1.3270345780684147, -1.5589839744772112, 1.2414736635261745, -1.4667294402006905, 1.3880171920550715, -1.533399261284528, 1.8040948422714405, -1.4189758124698626, 1.1408165196730307, -1.4607180452574147, 1.6048925584768932, -1.433392241781935};
        testConditionalSimulations(dataLikelihoodFactors, likelihoodDelegateFactors, diffusionModelFactor, dataModelFactor, rootPriorFactor, expectedTraits);

    }

    public void testLikelihoodDiagonalOURelaxedFactor() {
        System.out.println("\nTest Likelihood using diagonal Relaxed OU and factor:");

        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        ArbitraryBranchRates.BranchRateTransform transform = make(false, false, false);
        optimalTraitsModels.add(new ArbitraryBranchRates(treeModel,
                new Parameter.Default("rate.1", new double[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}),
                transform, false));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{-1.5})));

        DiagonalMatrix strengthOfSelectionMatrixParam = new DiagonalMatrix(new Parameter.Default(new double[]{1.5, 20.0}));

        DiffusionProcessDelegate diffusionProcessDelegate
                = new OUDiffusionModelDelegate(treeModel, diffusionModelFactor,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateFactors = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelFactor, rootPriorFactor,
                rateTransformation, rateModel, false);

        dataModelFactor.setLikelihoodDelegate(likelihoodDelegateFactors);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodFactors
                = new TreeDataLikelihood(likelihoodDelegateFactors, treeModel, rateModel);

        testLikelihood("likelihoodDiagonalOURelaxedFactor", dataModelFactor, dataLikelihoodFactors);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        double[] expectedTraits = new double[]{1.254609711504247, -1.1761389610061128, 1.305808317174877, -1.0709005664064868, 1.4571577865578949, -1.147788545433862, 1.7495515066909766, -0.9890375865560959, 1.0763987352911082, -1.0671848962791533, 1.5276137552159164, -0.9822950803483552};
        testConditionalSimulations(dataLikelihoodFactors, likelihoodDelegateFactors, diffusionModelFactor, dataModelFactor, rootPriorFactor, expectedTraits);

    }

    public void testLikelihoodFullOUFactor() {
        System.out.println("\nTest Likelihood using full OU and factor:");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{1.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{2.0})));

        Parameter[] strengthOfSelectionParameters = new Parameter[2];
        strengthOfSelectionParameters[0] = new Parameter.Default(new double[]{0.5, 0.05});
        strengthOfSelectionParameters[1] = new Parameter.Default(new double[]{0.05, 25.5});
        MatrixParameter strengthOfSelectionMatrixParam
                = new MatrixParameter("strengthOfSelectionMatrix", strengthOfSelectionParameters);

        DiffusionProcessDelegate diffusionProcessDelegate
                = new OUDiffusionModelDelegate(treeModel, diffusionModelFactor,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateFactors = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelFactor, rootPriorFactor,
                rateTransformation, rateModel, false);

        dataModelFactor.setLikelihoodDelegate(likelihoodDelegateFactors);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodFactors
                = new TreeDataLikelihood(likelihoodDelegateFactors, treeModel, rateModel);

        testLikelihood("likelihoodFullOUFactor", dataModelFactor, dataLikelihoodFactors);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        double[] expectedTraits = new double[]{0.4889263055380016, 1.866143125332805, 0.418517497134906, 1.9704713910092426, 0.5589398189745178, 1.8942177988462992, 0.9699471558714021, 2.0423474264069226, 0.32888191115035126, 1.9759942580128798, 0.8081782261607062, 2.0382998491256417};
        testConditionalSimulations(dataLikelihoodFactors, likelihoodDelegateFactors, diffusionModelFactor, dataModelFactor, rootPriorFactor, expectedTraits);

    }

    public void testLikelihoodFullOURelaxedFactor() {
        System.out.println("\nTest Likelihood using full Relaxed OU and factor:");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        ArbitraryBranchRates.BranchRateTransform transform = make(false, false, false);
        optimalTraitsModels.add(new ArbitraryBranchRates(treeModel,
                new Parameter.Default("rate.1", new double[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}),
                transform, false));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{1.5})));

        Parameter[] strengthOfSelectionParameters = new Parameter[2];
        strengthOfSelectionParameters[0] = new Parameter.Default(new double[]{0.5, 0.15});
        strengthOfSelectionParameters[1] = new Parameter.Default(new double[]{0.15, 25.5});
        MatrixParameter strengthOfSelectionMatrixParam
                = new MatrixParameter("strengthOfSelectionMatrix", strengthOfSelectionParameters);

        DiffusionProcessDelegate diffusionProcessDelegate
                = new OUDiffusionModelDelegate(treeModel, diffusionModelFactor,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParam));


        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateFactors = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelFactor, rootPriorFactor,
                rateTransformation, rateModel, false);

        dataModelFactor.setLikelihoodDelegate(likelihoodDelegateFactors);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodFactors
                = new TreeDataLikelihood(likelihoodDelegateFactors, treeModel, rateModel);

        testLikelihood("likelihoodFullRelaxedOUFactor", dataModelFactor, dataLikelihoodFactors);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        double[] expectedTraits = new double[]{0.6074917697347216, 1.4240248940048557, 0.5819299147586043, 1.5372685585813965, 0.7248840309526001, 1.4623057817648863, 1.0961030599134718, 1.6036947173685623, 0.442809377799788, 1.5374906895611142, 0.9206989848884695, 1.6011019729425795};
        testConditionalSimulations(dataLikelihoodFactors, likelihoodDelegateFactors, diffusionModelFactor, dataModelFactor, rootPriorFactor, expectedTraits);

    }

    public void testLikelihoodFullDiagonalOUFactor() {
        System.out.println("\nTest Likelihood comparing full and diagonal OU and factor:");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        ArbitraryBranchRates.BranchRateTransform transform = make(false, false, false);
        optimalTraitsModels.add(new ArbitraryBranchRates(treeModel,
                new Parameter.Default("rate.1", new double[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}),
                transform, false));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{1.5})));

        Parameter[] strengthOfSelectionParameters = new Parameter[2];
        strengthOfSelectionParameters[0] = new Parameter.Default(new double[]{0.5, 0.0});
        strengthOfSelectionParameters[1] = new Parameter.Default(new double[]{0.0, 1.5});
        MatrixParameter strengthOfSelectionMatrixParam
                = new MatrixParameter("strengthOfSelectionMatrix", strengthOfSelectionParameters);

        DiagonalMatrix strengthOfSelectionMatrixParamDiagonal
                = new DiagonalMatrix(new Parameter.Default(new double[]{0.5, 1.5}));

        DiffusionProcessDelegate diffusionProcessDelegate
                = new OUDiffusionModelDelegate(treeModel, diffusionModelFactor,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        DiffusionProcessDelegate diffusionProcessDelegateDiagonal
                = new OUDiffusionModelDelegate(treeModel, diffusionModelFactor,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParamDiagonal));

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateFactors = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelFactor, rootPriorFactor,
                rateTransformation, rateModel, false);

        ContinuousDataLikelihoodDelegate likelihoodDelegateFactorsDiagonal = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegateDiagonal, dataModelFactor, rootPriorFactor,
                rateTransformation, rateModel, false);


        dataModelFactor.setLikelihoodDelegate(likelihoodDelegateFactors);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodFactors
                = new TreeDataLikelihood(likelihoodDelegateFactors, treeModel, rateModel);

        TreeDataLikelihood dataLikelihoodFactorsDiagonal
                = new TreeDataLikelihood(likelihoodDelegateFactorsDiagonal, treeModel, rateModel);

        double likelihoodFactorData = dataLikelihoodFactors.getLogLikelihood();
        double likelihoodFactorDiffusion = dataModelFactor.getLogLikelihood();

        double likelihoodFactorDataDiagonal = dataLikelihoodFactorsDiagonal.getLogLikelihood();
        double likelihoodFactorDiffusionDiagonal = dataModelFactor.getLogLikelihood();


        assertEquals("likelihoodFullDiagonalOUFactor",
                format.format(likelihoodFactorData + likelihoodFactorDiffusion),
                format.format(likelihoodFactorDataDiagonal + likelihoodFactorDiffusionDiagonal));
    }

    private static double[] parseVectorLine(String s, String sep) {
        String[] vectorString = s.split(sep);
        double[] vec = new double[vectorString.length];
        vec[0] = Double.parseDouble(vectorString[0].substring(1));
        for (int i = 1; i < vec.length - 1; i++) {
            vec[i] = Double.parseDouble(vectorString[i]);
        }
        vec[vec.length - 1] = Double.parseDouble(vectorString[vec.length - 1].substring(0, vectorString[vec.length - 1].length() - 1));
        return vec;
    }

    private double[] parseVector(String s, String sep) {
        String[] vectorString = s.split(sep);
        double[] gradient = new double[vectorString.length];
        for (int i = 0; i < gradient.length; i++) {
            gradient[i] = Double.parseDouble(vectorString[i]);
        }
        return gradient;
    }

    private void testCMeans(TreeDataLikelihood dataLikelihood, String name, double[] partials) {
        String s = dataLikelihood.getReport();
        int offset = 0;
        int indBeg = 0;
        for (int tip = 0; tip < nTips; tip++) {
            indBeg = s.indexOf(name, indBeg + 1) + name.length() + 4;
            int indEnd = s.indexOf("]", indBeg);
            double[] vector = parseVector(s.substring(indBeg, indEnd), ",");
            for (int i = 0; i < vector.length; i++) {
//                System.out.println("cMean Mat: " + vector[i]);

                System.out.println("cMean preorder: " + partials[offset + i]);
                assertEquals("cMean " + tip + "; " + i,
                        format.format(partials[offset + i]),
                        format.format(vector[i]));
            }
            offset += dimTrait + PrecisionType.FULL.getMatrixLength(dimTrait);
        }
    }

    private void testCVariances(TreeDataLikelihood dataLikelihood, String name, double[] partials) {
        String s = dataLikelihood.getReport();
        int offset = 0;
        int indBeg = 0;
        for (int tip = 0; tip < nTips; tip++) {
            indBeg = s.indexOf(name, indBeg + 1) + name.length() + 3;
            int indEnd = s.indexOf("]", indBeg);
            double[] vector = parseVector(s.substring(indBeg, indEnd - 2), "\\s+|\\}\\n\\{ ");
            for (int i = 0; i < vector.length; i++) {
//                System.out.println("cMean Mat: " + vector[i]);
//                System.out.println("cMean preorder: " + partials[offset + dimTrait + dimTrait * dimTrait + i]);
                assertEquals("cVar " + tip + "; " + i,
                        format.format(partials[offset + dimTrait + dimTrait * dimTrait + i]),
                        format.format(vector[i]));
            }
            offset += dimTrait + PrecisionType.FULL.getMatrixLength(dimTrait);
        }
    }

    static double getLogDatumLikelihood(TreeDataLikelihood dataLikelihood) {
        String s = dataLikelihood.getReport();
        int indLikBeg = s.indexOf("logDatumLikelihood:") + 20;
        int indLikEnd = s.indexOf("\n", indLikBeg);
        char[] logDatumLikelihoodChar = new char[indLikEnd - indLikBeg + 1];
        s.getChars(indLikBeg, indLikEnd, logDatumLikelihoodChar, 0);
        return Double.parseDouble(String.valueOf(logDatumLikelihoodChar));
    }

    private void testLikelihood(String message, TreeDataLikelihood dataLikelihood) {
        double logDatumLikelihood = getLogDatumLikelihood(dataLikelihood);
        double integratedLikelihood = dataLikelihood.getLogLikelihood();
        assertEquals(message, format.format(logDatumLikelihood), format.format(integratedLikelihood));
        System.out.println(message + format.format(logDatumLikelihood));
    }

    static double getLogDatumLikelihood(IntegratedFactorAnalysisLikelihood dataModelFactor) {
        String sf = dataModelFactor.getReport();
        int indLikBegF = sf.indexOf("logMultiVariateNormalDensity = ") + 31;
        int indLikEndF = sf.indexOf("\n", indLikBegF);
        char[] logDatumLikelihoodCharF = new char[indLikEndF - indLikBegF + 1];
        sf.getChars(indLikBegF, indLikEndF, logDatumLikelihoodCharF, 0);
        return Double.parseDouble(String.valueOf(logDatumLikelihoodCharF));
    }

    private void testLikelihood(String message, IntegratedFactorAnalysisLikelihood dataModelFactor, TreeDataLikelihood dataLikelihoodFactors) {
        double logDatumLikelihoodFactor = getLogDatumLikelihood(dataModelFactor);

        double likelihoodFactorData = dataLikelihoodFactors.getLogLikelihood();
        double likelihoodFactorDiffusion = dataModelFactor.getLogLikelihood();

        assertEquals(message,
                format.format(logDatumLikelihoodFactor),
                format.format(likelihoodFactorData + likelihoodFactorDiffusion));
    }

    private void testConditionalMoments(TreeDataLikelihood dataLikelihood, ContinuousDataLikelihoodDelegate likelihoodDelegate) {
        new TreeTipGradient("" +
                "trait", dataLikelihood, likelihoodDelegate, null);
        TreeTraitLogger treeTraitLogger = new TreeTraitLogger(treeModel,
                new TreeTrait[]{dataLikelihood.getTreeTrait("fcd.trait")},
                TreeTraitLogger.NodeRestriction.EXTERNAL, false);

        String moments = treeTraitLogger.getReport();
        double[] partials = parseVector(moments, "\t");
        testCMeans(dataLikelihood, "cMean ", partials);
        testCVariances(dataLikelihood, "cVar ", partials);
    }

    private void testConditionalSimulations(TreeDataLikelihood dataLikelihood,
                                            ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                            MultivariateDiffusionModel diffusionModel,
                                            ContinuousTraitPartialsProvider dataModel,
                                            ConjugateRootTraitPrior rootPrior,
                                            double[] expectedTraits) {
        double[] traits = getConditionalSimulations(dataLikelihood, likelihoodDelegate, diffusionModel, dataModel, rootPrior, treeModel, rateTransformation);

        for (int i = 0; i < traits.length; i++) {
            assertEquals(format.format(expectedTraits[i]), format.format(traits[i]));
        }
    }

    static double[] getConditionalSimulations(TreeDataLikelihood dataLikelihood,
                                              ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                              MultivariateDiffusionModel diffusionModel,
                                              ContinuousTraitPartialsProvider dataModel,
                                              ConjugateRootTraitPrior rootPrior,
                                              TreeModel treeModel,
                                              ContinuousRateTransformation rateTransformation) {
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihood, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        return parseVectorLine(treeTrait[0].getTraitString(treeModel, null), ",");
    }
}