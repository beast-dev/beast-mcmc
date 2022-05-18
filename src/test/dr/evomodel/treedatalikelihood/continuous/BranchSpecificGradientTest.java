package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.List;

import static dr.evomodel.branchratemodel.ArbitraryBranchRates.make;


/**
 * @author Paul Bastide
 */

public class BranchSpecificGradientTest extends ContinuousTraitTest {

    public BranchSpecificGradientTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();
    }

    public void testRateGradient() {
        System.out.println("\nTest Likelihood using vanilla BM:");

        // Diffusion
        List<BranchRateModel> driftModels = new ArrayList<BranchRateModel>();
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{1.0})));
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{2.0})));
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.3", new double[]{-2.0})));
        DiffusionProcessDelegate diffusionProcessDelegate
                = new DriftDiffusionModelDelegate(treeModel, diffusionModel, driftModels);

        // Rates
        ArbitraryBranchRates.BranchRateTransform transform = make(false, true, false, null, null);
        Parameter branchRates = new Parameter.Default(new double[]{10, 9, 8, 7, 6, 5, 4, 3, 2, 1});
        ArbitraryBranchRates rateModel = new ArbitraryBranchRates(treeModel, branchRates, transform, false);

        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel, false);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihood = new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel);

        // Gradient (Rates)
        BranchRateGradient branchGradient1
                = new BranchRateGradient("trait", dataLikelihood, likelihoodDelegate, branchRates);
        double[] gradient1 = branchGradient1.getGradientLogDensity();

        // Gradient (Specific)
        ContinuousTraitGradientForBranch.RateGradient traitGradient =
                new ContinuousTraitGradientForBranch.RateGradient(dimTrait, treeModel, rateModel, diffusionProcessDelegate);
        BranchSpecificGradient branchGradient2
                = new BranchSpecificGradient("trait", dataLikelihood, likelihoodDelegate,
                traitGradient, branchRates);
        double[] gradient2 = branchGradient2.getGradientLogDensity();
        double[] numericalGradient = branchGradient1.getNumericalGradient();

        System.err.println("\tGradient with rate method    = " + new dr.math.matrixAlgebra.Vector(gradient1));
        System.err.println("\tGradient with general method = " + new dr.math.matrixAlgebra.Vector(gradient2));
        System.err.println("\tNumerical gradient           = " + new dr.math.matrixAlgebra.Vector(numericalGradient));

        assertEquals("length", gradient1.length, gradient2.length);
        for (int i = 0; i < gradient1.length; i++) {
            assertEquals("numeric " + i,
                    gradient1[i],
                    numericalGradient[i], 1E-4);
        }
        for (int i = 0; i < gradient1.length; i++) {
            assertEquals("gradient " + i,
                    format.format(gradient1[i]),
                    format.format(gradient2[i]));
        }
    }

    public void testRateGradientDriftScaled() {
        System.out.println("\nTest Likelihood using vanilla BM:");

        // Diffusion
        List<BranchRateModel> driftModels = new ArrayList<BranchRateModel>();
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{1.0})));
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{2.0})));
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.3", new double[]{-2.0})));
        DiffusionProcessDelegate diffusionProcessDelegate
                = new DriftDiffusionModelDelegate(treeModel, diffusionModel, driftModels, true);

        // Rates
        ArbitraryBranchRates.BranchRateTransform transform = make(false, true, false, null, null);
        Parameter branchRates = new Parameter.Default(new double[]{10, 9, 8, 7, 6, 5, 4, 3, 2, 1});
        ArbitraryBranchRates rateModel = new ArbitraryBranchRates(treeModel, branchRates, transform, false);

        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel, false);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihood = new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel);

        // Gradient (Rates)
        BranchRateGradient branchGradient1
                = new BranchRateGradient("trait", dataLikelihood, likelihoodDelegate, branchRates);
        double[] gradient1 = branchGradient1.getGradientLogDensity();

        // Gradient (Specific)
        ContinuousTraitGradientForBranch.RateGradient traitGradient =
                new ContinuousTraitGradientForBranch.RateGradient(dimTrait, treeModel, rateModel, diffusionProcessDelegate);
        BranchSpecificGradient branchGradient2
                = new BranchSpecificGradient("trait", dataLikelihood, likelihoodDelegate,
                traitGradient, branchRates);
        double[] gradient2 = branchGradient2.getGradientLogDensity();
        double[] numericalGradient = branchGradient1.getNumericalGradient();

        System.err.println("\tGradient with rate method    = " + new dr.math.matrixAlgebra.Vector(gradient1));
        System.err.println("\tGradient with general method = " + new dr.math.matrixAlgebra.Vector(gradient2));
        System.err.println("\tNumerical gradient           = " + new dr.math.matrixAlgebra.Vector(numericalGradient));

        assertEquals("length", gradient1.length, gradient2.length);
        for (int i = 0; i < gradient1.length; i++) {
            assertEquals("numeric " + i,
                    gradient1[i],
                    numericalGradient[i], 1E-4);
        }
        for (int i = 0; i < gradient1.length; i++) {
            assertEquals("gradient " + i,
                    format.format(gradient1[i]),
                    format.format(gradient2[i]));
        }
    }
}