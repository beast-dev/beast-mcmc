package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.datatype.Nucleotides;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.inference.model.*;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;


/**
 * @author Paul Bastide
 */

public class ContinuousDataLikelihoodDelegateTest extends TraceCorrelationAssert {

    private MultivariateDiffusionModel diffusionModel;
    private ContinuousTraitPartialsProvider dataModel;
    private ConjugateRootTraitPrior rootPrior;

    private MultivariateDiffusionModel diffusionModelFactor;
    private IntegratedFactorAnalysisLikelihood dataModelFactor;
    private ConjugateRootTraitPrior rootPriorFactor;

    private NumberFormat format = NumberFormat.getNumberInstance(Locale.ENGLISH);

    public ContinuousDataLikelihoodDelegateTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        format.setMaximumFractionDigits(5);

        // Tree
        createAlignment(PRIMATES_TAXON_SEQUENCE, Nucleotides.INSTANCE);
        treeModel = createPrimateTreeModel();

        // Data
        Parameter[] dataTraits = new Parameter[6];
        dataTraits[0] = new Parameter.Default("human", new double[]{1.0, 2.0});
        dataTraits[1] = new Parameter.Default("chimp", new double[]{10.0, 12.0});
        dataTraits[2] = new Parameter.Default("bonobo", new double[]{0.5, 2.0});
        dataTraits[3] = new Parameter.Default("gorilla", new double[]{2.0, 5.0});
        dataTraits[4] = new Parameter.Default("orangutan", new double[]{11.0, 1.0});
        dataTraits[5] = new Parameter.Default("siamang", new double[]{1.0, 2.5});
        CompoundParameter traitParameter = new CompoundParameter("trait", dataTraits);

        List<Integer> missingIndices = new ArrayList<Integer>();
        traitParameter.setParameterValue(2, 0);
        missingIndices.add(2);
        missingIndices.add(3);
        missingIndices.add(7);

        //// Standard Model //// ***************************************************************************************

        // Diffusion
        Parameter[] precisionParameters = new Parameter[2];
        precisionParameters[0] = new Parameter.Default(new double[]{1.0, 0.1});
        precisionParameters[1] = new Parameter.Default(new double[]{0.1, 2.0});
        MatrixParameterInterface diffusionPrecisionMatrixParameter = new MatrixParameter("precisionMatrix", precisionParameters);
        diffusionModel = new MultivariateDiffusionModel(diffusionPrecisionMatrixParameter);

        PrecisionType precisionType = PrecisionType.FULL;

        // Root prior
        rootPrior = new ConjugateRootTraitPrior(new double[]{-1.0, -3.0}, 10.0, true);

        // Data Model
        dataModel = new ContinuousTraitDataModel("dataModel",
                traitParameter,
                missingIndices, true,
                2, precisionType);

        //// Factor Model //// *****************************************************************************************
        // Diffusion
        Parameter[] precisionParametersFactor = new Parameter[1];
        precisionParametersFactor[0] = new Parameter.Default(new double[]{1.0});
        MatrixParameterInterface diffusionPrecisionMatrixParameterFactor = new MatrixParameter("precisionMatrixFactor", precisionParametersFactor);
        diffusionModelFactor = new MultivariateDiffusionModel(diffusionPrecisionMatrixParameterFactor);

        // Root prior
        rootPriorFactor = new ConjugateRootTraitPrior(new double[]{-1.0}, 10.0, true);

        // Error model
        Parameter factorPrecisionParameters = new Parameter.Default("factorPrecision", new double[]{1.0, 5.0});

        // Loadings
        Parameter[] loadingsParameters = new Parameter[1];
        loadingsParameters[0] = new Parameter.Default(new double[]{1.0, 2.0});
        MatrixParameterInterface loadingsMatrixParameters = new MatrixParameter("loadings", loadingsParameters);

        dataModelFactor = new IntegratedFactorAnalysisLikelihood("dataModelFactors",
                traitParameter,
                missingIndices,
                loadingsMatrixParameters,
                factorPrecisionParameters, 0.0);
    }

    public void testLikelihoodBM() {
        System.out.println("\nTest Likelihood using vanilla BM:");

        // Diffusion
        DiffusionProcessDelegate diffusionProcessDelegate = new HomogeneousDiffusionModelDelegate(treeModel, diffusionModel);

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel, true);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihood = new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel);

        String s = dataLikelihood.getReport();
        int indLikBeg = s.indexOf("logDatumLikelihood:") + 20;
        int indLikEnd = s.indexOf("\n", indLikBeg);
        char[] logDatumLikelihoodChar = new char[indLikEnd - indLikBeg + 1];
        s.getChars(indLikBeg, indLikEnd, logDatumLikelihoodChar, 0);
        double logDatumLikelihood = Double.parseDouble(String.valueOf(logDatumLikelihoodChar));

        double integratedLikelihood = dataLikelihood.getLogLikelihood();

        assertEquals("likelihoodBM", format.format(logDatumLikelihood), format.format(integratedLikelihood));

        System.out.println("likelihoodBM: " + format.format(logDatumLikelihood));
    }

    public void testLikelihoodDrift() {
        System.out.println("\nTest Likelihood using Drifted BM:");

        // Diffusion
        List<BranchRateModel> driftModels = new ArrayList<BranchRateModel>();
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{1.0})));
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{2.0})));
        DiffusionProcessDelegate diffusionProcessDelegate = new DriftDiffusionModelDelegate(treeModel, diffusionModel, driftModels);

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel, false);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihood = new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel);

        String s = dataLikelihood.getReport();
        int indLikBeg = s.indexOf("logDatumLikelihood:") + 20;
        int indLikEnd = s.indexOf("\n", indLikBeg);
        char[] logDatumLikelihoodChar = new char[indLikEnd - indLikBeg + 1];
        s.getChars(indLikBeg, indLikEnd, logDatumLikelihoodChar, 0);
        double logDatumLikelihood = Double.parseDouble(String.valueOf(logDatumLikelihoodChar));

        assertEquals("likelihoodDrift", format.format(logDatumLikelihood), format.format(dataLikelihood.getLogLikelihood()));

        System.out.println("likelihoodDrift: " + format.format(logDatumLikelihood));
    }

    public void testLikelihoodDiagonalOU() {
        System.out.println("\nTest Likelihood using Diagonal OU:");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{1.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{2.0})));

        DiagonalMatrix strengthOfSelectionMatrixParam = new DiagonalMatrix(new Parameter.Default(new double[]{1.0, 2.0}));

        DiffusionProcessDelegate diffusionProcessDelegate = new DiagonalOrnsteinUhlenbeckDiffusionModelDelegate(treeModel, diffusionModel, optimalTraitsModels, strengthOfSelectionMatrixParam);

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel, false);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihood = new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel);

        String s = dataLikelihood.getReport();
        int indLikBeg = s.indexOf("logDatumLikelihood:") + 20;
        int indLikEnd = s.indexOf("\n", indLikBeg);
        char[] logDatumLikelihoodChar = new char[indLikEnd - indLikBeg + 1];
        s.getChars(indLikBeg, indLikEnd, logDatumLikelihoodChar, 0);
        double logDatumLikelihood = Double.parseDouble(String.valueOf(logDatumLikelihoodChar));

        assertEquals("likelihoodDiagOU", format.format(logDatumLikelihood), format.format(dataLikelihood.getLogLikelihood()));

        System.out.println("likelihoodDiagOU: " + format.format(logDatumLikelihood));
    }

    public void testLikelihoodFullOU() {
        System.out.println("\nTest Likelihood using Full OU:");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{1.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{2.0})));

        Parameter[] strengthOfSelectionParameters = new Parameter[2];
        strengthOfSelectionParameters[0] = new Parameter.Default(new double[]{0.5, 0.2});
        strengthOfSelectionParameters[1] = new Parameter.Default(new double[]{0.2, 2.0});
        MatrixParameter strengthOfSelectionMatrixParam = new MatrixParameter("strengthOfSelectionMatrix", strengthOfSelectionParameters);

        DiffusionProcessDelegate diffusionProcessDelegate = new OrnsteinUhlenbeckDiffusionModelDelegate(treeModel, diffusionModel, optimalTraitsModels, strengthOfSelectionMatrixParam);

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel, false);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihood = new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel);

        String s = dataLikelihood.getReport();
        int indLikBeg = s.indexOf("logDatumLikelihood:") + 20;
        int indLikEnd = s.indexOf("\n", indLikBeg);
        char[] logDatumLikelihoodChar = new char[indLikEnd - indLikBeg + 1];
        s.getChars(indLikBeg, indLikEnd, logDatumLikelihoodChar, 0);
        double logDatumLikelihood = Double.parseDouble(String.valueOf(logDatumLikelihoodChar));

        assertEquals("likelihoodFullOU", format.format(logDatumLikelihood), format.format(dataLikelihood.getLogLikelihood()));

        System.out.println("likelihoodFullOU: " + format.format(logDatumLikelihood));
    }

    //// Factor Model //// *********************************************************************************************

    public void testLikelihoodBMFactor() {
        System.out.println("\nTest Likelihood using vanilla BM and factor:");

        // Diffusion
        DiffusionProcessDelegate diffusionProcessDelegate = new HomogeneousDiffusionModelDelegate(treeModel, diffusionModelFactor);

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateFactors = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelFactor, rootPriorFactor, rateTransformation, rateModel, true);

        dataModelFactor.setLikelihoodDelegate(likelihoodDelegateFactors);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodFactors = new TreeDataLikelihood(likelihoodDelegateFactors, treeModel, rateModel);

        String sf = dataModelFactor.getReport();
        int indLikBegF = sf.indexOf("logMultiVariateNormalDensity = ") + 31;
        int indLikEndF = sf.indexOf("\n", indLikBegF);
        char[] logDatumLikelihoodCharF = new char[indLikEndF - indLikBegF + 1];
        sf.getChars(indLikBegF, indLikEndF, logDatumLikelihoodCharF, 0);
        double logDatumLikelihoodFactor = Double.parseDouble(String.valueOf(logDatumLikelihoodCharF));

        double likelihoodFactorData = dataLikelihoodFactors.getLogLikelihood();
        double likelihoodFactorDiffusion = dataModelFactor.getLogLikelihood();


        assertEquals("likelihoodBMFactor", format.format(logDatumLikelihoodFactor), format.format(likelihoodFactorData + likelihoodFactorDiffusion));
    }

    public void testLikelihoodDriftFactor() {
        System.out.println("\nTest Likelihood using drifted BM and factor:");

        // Diffusion
        List<BranchRateModel> driftModels = new ArrayList<BranchRateModel>();
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{0.0})));
        DiffusionProcessDelegate diffusionProcessDelegate = new DriftDiffusionModelDelegate(treeModel, diffusionModelFactor, driftModels);

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateFactors = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelFactor, rootPriorFactor, rateTransformation, rateModel, false);

        dataModelFactor.setLikelihoodDelegate(likelihoodDelegateFactors);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodFactors = new TreeDataLikelihood(likelihoodDelegateFactors, treeModel, rateModel);

        String sf = dataModelFactor.getReport();
        int indLikBegF = sf.indexOf("logMultiVariateNormalDensity = ") + 31;
        int indLikEndF = sf.indexOf("\n", indLikBegF);
        char[] logDatumLikelihoodCharF = new char[indLikEndF - indLikBegF + 1];
        sf.getChars(indLikBegF, indLikEndF, logDatumLikelihoodCharF, 0);
        double logDatumLikelihoodFactor = Double.parseDouble(String.valueOf(logDatumLikelihoodCharF));

        double likelihoodFactorData = dataLikelihoodFactors.getLogLikelihood();
        double likelihoodFactorDiffusion = dataModelFactor.getLogLikelihood();


        assertEquals("likelihoodDriftFactor", format.format(logDatumLikelihoodFactor), format.format(likelihoodFactorData + likelihoodFactorDiffusion));
    }

    public void testLikelihoodDiagonalOUFactor() {
        System.out.println("\nTest Likelihood using diagonal OU and factor:");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{1.0})));

        DiagonalMatrix strengthOfSelectionMatrixParam = new DiagonalMatrix(new Parameter.Default(new double[]{1.5}));

        DiffusionProcessDelegate diffusionProcessDelegate = new DiagonalOrnsteinUhlenbeckDiffusionModelDelegate(treeModel, diffusionModelFactor, optimalTraitsModels, strengthOfSelectionMatrixParam);

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateFactors = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelFactor, rootPriorFactor, rateTransformation, rateModel, false);

        dataModelFactor.setLikelihoodDelegate(likelihoodDelegateFactors);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodFactors = new TreeDataLikelihood(likelihoodDelegateFactors, treeModel, rateModel);

        String sf = dataModelFactor.getReport();
        int indLikBegF = sf.indexOf("logMultiVariateNormalDensity = ") + 31;
        int indLikEndF = sf.indexOf("\n", indLikBegF);
        char[] logDatumLikelihoodCharF = new char[indLikEndF - indLikBegF + 1];
        sf.getChars(indLikBegF, indLikEndF, logDatumLikelihoodCharF, 0);
        double logDatumLikelihoodFactor = Double.parseDouble(String.valueOf(logDatumLikelihoodCharF));

        double likelihoodFactorData = dataLikelihoodFactors.getLogLikelihood();
        double likelihoodFactorDiffusion = dataModelFactor.getLogLikelihood();


        assertEquals("likelihoodDiagonalOUFactor", format.format(logDatumLikelihoodFactor), format.format(likelihoodFactorData + likelihoodFactorDiffusion));
    }

    public void testLikelihoodFullOUFactor() {
        System.out.println("\nTest Likelihood using full OU and factor:");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{1.0})));

        Parameter[] strengthOfSelectionParameters = new Parameter[1];
        strengthOfSelectionParameters[0] = new Parameter.Default(new double[]{0.5});
        MatrixParameter strengthOfSelectionMatrixParam = new MatrixParameter("strengthOfSelectionMatrix", strengthOfSelectionParameters);

        DiffusionProcessDelegate diffusionProcessDelegate = new OrnsteinUhlenbeckDiffusionModelDelegate(treeModel, diffusionModelFactor, optimalTraitsModels, strengthOfSelectionMatrixParam);

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateFactors = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelFactor, rootPriorFactor, rateTransformation, rateModel, false);

        dataModelFactor.setLikelihoodDelegate(likelihoodDelegateFactors);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodFactors = new TreeDataLikelihood(likelihoodDelegateFactors, treeModel, rateModel);

        String sf = dataModelFactor.getReport();
        int indLikBegF = sf.indexOf("logMultiVariateNormalDensity = ") + 31;
        int indLikEndF = sf.indexOf("\n", indLikBegF);
        char[] logDatumLikelihoodCharF = new char[indLikEndF - indLikBegF + 1];
        sf.getChars(indLikBegF, indLikEndF, logDatumLikelihoodCharF, 0);
        double logDatumLikelihoodFactor = Double.parseDouble(String.valueOf(logDatumLikelihoodCharF));

        double likelihoodFactorData = dataLikelihoodFactors.getLogLikelihood();
        double likelihoodFactorDiffusion = dataModelFactor.getLogLikelihood();


        assertEquals("likelihoodFullOUFactor", format.format(logDatumLikelihoodFactor), format.format(likelihoodFactorData + likelihoodFactorDiffusion));
    }

}