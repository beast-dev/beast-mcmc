package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.datatype.Nucleotides;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.inference.model.CompoundParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import static dr.evomodel.branchratemodel.ArbitraryBranchRates.make;


/**
 * @author Paul Bastide
 */

public class BranchSpecificGradientTest extends TraceCorrelationAssert {

    private int dimTrait;

    private MultivariateDiffusionModel diffusionModel;
    private ContinuousTraitPartialsProvider dataModel;
    private ConjugateRootTraitPrior rootPrior;

    private NumberFormat format = NumberFormat.getNumberInstance(Locale.ENGLISH);

    public BranchSpecificGradientTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        dimTrait = 3;

        format.setMaximumFractionDigits(5);

        // Tree
        createAlignment(PRIMATES_TAXON_SEQUENCE, Nucleotides.INSTANCE);
        treeModel = createPrimateTreeModel();

        // Data
        Parameter[] dataTraits = new Parameter[6];
        dataTraits[0] = new Parameter.Default("human", new double[]{-1.0, 2.0, 3.0});
        dataTraits[1] = new Parameter.Default("chimp", new double[]{10.0, 12.0, 14.0});
        dataTraits[2] = new Parameter.Default("bonobo", new double[]{0.5, -2.0, 5.5});
        dataTraits[3] = new Parameter.Default("gorilla", new double[]{2.0, 5.0, -8.0});
        dataTraits[4] = new Parameter.Default("orangutan", new double[]{11.0, 1.0, -1.5});
        dataTraits[5] = new Parameter.Default("siamang", new double[]{1.0, 2.5, 4.0});
        CompoundParameter traitParameter = new CompoundParameter("trait", dataTraits);

        List<Integer> missingIndices = new ArrayList<Integer>();
        traitParameter.setParameterValue(2, 0);
        missingIndices.add(3);
        missingIndices.add(4);
        missingIndices.add(5);
        missingIndices.add(7);

        //// Standard Model //// ***************************************************************************************

        // Diffusion
        Parameter[] precisionParameters = new Parameter[dimTrait];
        precisionParameters[0] = new Parameter.Default(new double[]{1.0, 0.1, 0.2});
        precisionParameters[1] = new Parameter.Default(new double[]{0.1, 2.0, 0.0});
        precisionParameters[2] = new Parameter.Default(new double[]{0.2, 0.0, 3.0});
        MatrixParameterInterface diffusionPrecisionMatrixParameter
                = new MatrixParameter("precisionMatrix", precisionParameters);
        diffusionModel = new MultivariateDiffusionModel(diffusionPrecisionMatrixParameter);

        PrecisionType precisionType = PrecisionType.FULL;

        // Root prior
        Parameter rootMean = new Parameter.Default(new double[]{-1.0, -3.0, 2.5});
        Parameter rootSampleSize = new Parameter.Default(10.0);
        rootPrior = new ConjugateRootTraitPrior(rootMean, rootSampleSize);

        // Data Model
        dataModel = new ContinuousTraitDataModel("dataModel",
                traitParameter,
                missingIndices, true,
                dimTrait, precisionType);
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
        ArbitraryBranchRates.BranchRateTransform transform = make(false, true, null, null);
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
                new ContinuousTraitGradientForBranch.RateGradient(dimTrait, treeModel, rateModel);
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