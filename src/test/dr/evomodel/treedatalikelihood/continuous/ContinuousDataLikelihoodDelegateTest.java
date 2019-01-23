package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodel.treedatalikelihood.preorder.MultivariateConditionalOnTipsRealizedDelegate;
import dr.evomodel.treedatalikelihood.preorder.ProcessSimulationDelegate;
import dr.evomodel.treelikelihood.utilities.TreeTraitLogger;
import dr.inference.model.*;
import dr.math.MathUtils;
import dr.xml.AttributeParser;
import dr.xml.XMLParser;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.io.StringReader;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import static dr.evomodel.branchratemodel.ArbitraryBranchRates.make;


/**
 * @author Paul Bastide
 */

public class ContinuousDataLikelihoodDelegateTest extends TraceCorrelationAssert {

    private int dimTrait;
    private int nTips;

    private MultivariateDiffusionModel diffusionModel;
    private ContinuousTraitPartialsProvider dataModel;
    private ContinuousTraitPartialsProvider dataModelIntegrated;
    private ConjugateRootTraitPrior rootPrior;
    private ConjugateRootTraitPrior rootPriorIntegrated;

    private MultivariateDiffusionModel diffusionModelFactor;
    private IntegratedFactorAnalysisLikelihood dataModelFactor;
    private ConjugateRootTraitPrior rootPriorFactor;

    private NumberFormat format = NumberFormat.getNumberInstance(Locale.ENGLISH);

    public ContinuousDataLikelihoodDelegateTest(String name) {
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
        nTips = 6;
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
        String s = "<beast>\n" +
                "    <conjugateRootPrior>\n" +
                "        <meanParameter>\n" +
                "            <parameter id=\"meanRoot\"  value=\"-1.0 -3.0 2.5\"/>\n" +
                "        </meanParameter>\n" +
                "        <priorSampleSize>\n" +
                "            <parameter id=\"sampleSizeRoot\" value=\"10.0\"/>\n" +
                "        </priorSampleSize>\n" +
                "    </conjugateRootPrior>\n" +
                "</beast>";
        XMLParser parser = new XMLParser(true, true, true, null);
        parser.addXMLObjectParser(new AttributeParser());
        parser.addXMLObjectParser(new ParameterParser());
        parser.parse(new StringReader(s), true);
        rootPrior = ConjugateRootTraitPrior.parseConjugateRootTraitPrior(parser.getRoot(), dimTrait);

        // Data Model
        dataModel = new ContinuousTraitDataModel("dataModel",
                traitParameter,
                missingIndices, true,
                dimTrait, precisionType);

        //// Factor Model //// *****************************************************************************************
        // Diffusion
        Parameter[] precisionParametersFactor = new Parameter[2];
        precisionParametersFactor[0] = new Parameter.Default(new double[]{1.0, 0.1});
        precisionParametersFactor[1] = new Parameter.Default(new double[]{0.1, 1.5});
        MatrixParameterInterface diffusionPrecisionMatrixParameterFactor
                = new MatrixParameter("precisionMatrixFactor", precisionParametersFactor);
        diffusionModelFactor = new MultivariateDiffusionModel(diffusionPrecisionMatrixParameterFactor);

        // Root prior
        String sFactor = "<beast>\n" +
                "    <conjugateRootPrior>\n" +
                "        <meanParameter>\n" +
                "            <parameter id=\"meanRoot\"  value=\"-1.0 2.0\"/>\n" +
                "        </meanParameter>\n" +
                "        <priorSampleSize>\n" +
                "            <parameter id=\"sampleSizeRoot\" value=\"10.0\"/>\n" +
                "        </priorSampleSize>\n" +
                "    </conjugateRootPrior>\n" +
                "</beast>";
        XMLParser parserFactor = new XMLParser(true, true, true, null);
        parserFactor.addXMLObjectParser(new AttributeParser());
        parserFactor.addXMLObjectParser(new ParameterParser());
        parserFactor.parse(new StringReader(sFactor), true);
        rootPriorFactor = ConjugateRootTraitPrior.parseConjugateRootTraitPrior(parserFactor.getRoot(), 2);
//        rootPriorFactor = new ConjugateRootTraitPrior(new double[]{-1.0, 2.0}, 10.0, true);

        // Error model
        Parameter[] factorPrecision = new Parameter[3];
        factorPrecision[0] = new Parameter.Default(new double[]{1.0, 0.5, 0.2});
        factorPrecision[1] = new Parameter.Default(new double[]{0.5, 5.0, 0.0});
        factorPrecision[2] = new Parameter.Default(new double[]{0.2, 0.0, 0.5});
        MatrixParameterInterface factorPrecisionParameters
                = new MatrixParameter("factorPrecision", factorPrecision);

//        MatrixParameterInterface factorPrecisionParameters = new DiagonalMatrix(new Parameter.Default("factorPrecision", new double[]{1.0, 5.0, 0.5}));

        // Loadings
        Parameter[] loadingsParameters = new Parameter[2];
        loadingsParameters[0] = new Parameter.Default(new double[]{1.0, 2.0, 3.0});
        loadingsParameters[1] = new Parameter.Default(new double[]{0.0, 0.5, 1.0});
        MatrixParameterInterface loadingsMatrixParameters = new MatrixParameter("loadings", loadingsParameters);

        dataModelFactor = new IntegratedFactorAnalysisLikelihood("dataModelFactors",
                traitParameter,
                missingIndices,
                loadingsMatrixParameters,
                factorPrecisionParameters, 0.0, null);

        //// Integrated Process //// ***********************************************************************************
        // Data Model
        dataModelIntegrated = new IntegratedProcessTraitDataModel("dataModelIntegrated",
                traitParameter,
                missingIndices, true,
                dimTrait, precisionType);

        // Root prior
        String sI = "<beast>\n" +
                "    <conjugateRootPrior>\n" +
                "        <meanParameter>\n" +
                "            <parameter id=\"meanRoot\"  value=\"0.0 1.2 -0.5 -1.0 -3.0 2.5\"/>\n" +
                "        </meanParameter>\n" +
                "        <priorSampleSize>\n" +
                "            <parameter id=\"sampleSizeRoot\" value=\"10.0\"/>\n" +
                "        </priorSampleSize>\n" +
                "    </conjugateRootPrior>\n" +
                "</beast>";
        XMLParser parserI = new XMLParser(true, true, true, null);
        parserI.addXMLObjectParser(new AttributeParser());
        parserI.addXMLObjectParser(new ParameterParser());
        parserI.parse(new StringReader(sI), true);
        rootPriorIntegrated = ConjugateRootTraitPrior.parseConjugateRootTraitPrior(parserI.getRoot(), 2 * dimTrait);
    }

    public void testLikelihoodBM() {
        System.out.println("\nTest Likelihood using vanilla BM:");

        // Diffusion
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

        String s = dataLikelihood.getReport();
        int indLikBeg = s.indexOf("logDatumLikelihood:") + 20;
        int indLikEnd = s.indexOf("\n", indLikBeg);
        char[] logDatumLikelihoodChar = new char[indLikEnd - indLikBeg + 1];
        s.getChars(indLikBeg, indLikEnd, logDatumLikelihoodChar, 0);
        double logDatumLikelihood = Double.parseDouble(String.valueOf(logDatumLikelihoodChar));

        double integratedLikelihood = dataLikelihood.getLogLikelihood();

        assertEquals("likelihoodBM", format.format(logDatumLikelihood), format.format(integratedLikelihood));

        System.out.println("likelihoodBM: " + format.format(logDatumLikelihood));

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihood, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{-1.0,2.0,0.0,0.45807521679597646,2.6505355982097605,3.4693334367360538,0.5,2.64206285585883,5.5,2.0,5.0,-8.0,11.0,1.0,-1.5,1.0,2.5,4.0};
        double[] traits = parseVectorLine(treeTrait[0].getTraitString(treeModel, null), ",");
        for (int i = 0; i < traits.length; i++) {
            assertEquals(format.format(expectedTraits[i]), format.format(traits[i]));
        }

        // Conditional moments (preorder)
        new TreeTipGradient("" +
                "trait", dataLikelihood, likelihoodDelegate, null);
        TreeTraitLogger treeTraitLogger = new TreeTraitLogger(treeModel,
                new TreeTrait[]{dataLikelihood.getTreeTrait("fcd.trait")},
                TreeTraitLogger.NodeRestriction.EXTERNAL, false);

        String moments = treeTraitLogger.getReport();
        double[] partials = parseVector(moments, "\t");
        testCMeans(s, "cMean ", partials);
        testCVariances(s, "cVar ", partials);

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

        assertEquals("likelihoodDrift",
                format.format(logDatumLikelihood),
                format.format(dataLikelihood.getLogLikelihood()));

        System.out.println("likelihoodDrift: " + format.format(logDatumLikelihood));

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihood, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{-1.0,2.0,0.0,0.5457621072639138,3.28662834718796,3.2939596558001845,0.5,1.0742799493604265,5.5,2.0,5.0,-8.0,11.0,1.0,-1.5,1.0,2.5,4.0};
        double[] traits = parseVectorLine(treeTrait[0].getTraitString(treeModel, null), ",");
        for (int i = 0; i < traits.length; i++) {
            assertEquals(format.format(expectedTraits[i]), format.format(traits[i]));
        }

        // Conditional moments (preorder)
        new TreeTipGradient("" +
                "trait", dataLikelihood, likelihoodDelegate, null);
        TreeTraitLogger treeTraitLogger = new TreeTraitLogger(treeModel,
                new TreeTrait[]{dataLikelihood.getTreeTrait("fcd.trait")},
                TreeTraitLogger.NodeRestriction.EXTERNAL, false);

        String moments = treeTraitLogger.getReport();
        double[] partials = parseVector(moments, "\t");
        testCMeans(s, "cMean ", partials);
        testCVariances(s, "cVar ", partials);
    }

    public void testLikelihoodDriftRelaxed() {
        System.out.println("\nTest Likelihood using Drifted relaxed BM:");

        // Diffusion
        List<BranchRateModel> driftModels = new ArrayList<BranchRateModel>();
        ArbitraryBranchRates.BranchRateTransform transform = make(false, false);
        driftModels.add(new ArbitraryBranchRates(treeModel,
                new Parameter.Default("rate.1", new double[]{0, 100, 200, 300, 400, 500, 600, 700, 800, 900}),
                transform, false));
        driftModels.add(new ArbitraryBranchRates(treeModel,
                new Parameter.Default("rate.2", new double[]{0, -100, 200, -300, 400, -500, 600, -700, 800, -900}),
                transform, false));
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.3", new double[]{-2.0})));
        DiffusionProcessDelegate diffusionProcessDelegate
                = new DriftDiffusionModelDelegate(treeModel, diffusionModel, driftModels);

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

        assertEquals("likelihoodDriftRelaxed",
                format.format(logDatumLikelihood),
                format.format(dataLikelihood.getLogLikelihood()));

        System.out.println("likelihoodDriftRelaxed: " + format.format(logDatumLikelihood));

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihood, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{-1.0,2.0,0.0,2.843948876154644,10.866053719140933,3.467579698926694,0.5,12.000214659757933,5.5,2.0,5.0,-8.0,11.0,1.0,-1.5,1.0,2.5,4.0};
        double[] traits = parseVectorLine(treeTrait[0].getTraitString(treeModel, null), ",");
        for (int i = 0; i < traits.length; i++) {
            assertEquals(format.format(expectedTraits[i]), format.format(traits[i]));
        }
        // Conditional moments (preorder)
        new TreeTipGradient("" +
                "trait", dataLikelihood, likelihoodDelegate, null);
        TreeTraitLogger treeTraitLogger = new TreeTraitLogger(treeModel,
                new TreeTrait[]{dataLikelihood.getTreeTrait("fcd.trait")},
                TreeTraitLogger.NodeRestriction.EXTERNAL, false);

        String moments = treeTraitLogger.getReport();
        double[] partials = parseVector(moments, "\t");
        testCMeans(s, "cMean ", partials);
        testCVariances(s, "cVar ", partials);
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

        assertEquals("likelihoodDiagonalOU",
                format.format(logDatumLikelihood),
                format.format(dataLikelihood.getLogLikelihood()));

        System.out.println("likelihoodDiagonalOU: " + format.format(logDatumLikelihood));

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihood, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{-1.0,2.0,0.0,1.0369622398437415,2.065450266793184,0.6174755164694558,0.5,2.0829935706195615,5.5,2.0,5.0,-8.0,11.0,1.0,-1.5,1.0,2.5,4.0};
        double[] traits = parseVectorLine(treeTrait[0].getTraitString(treeModel, null), ",");
        for (int i = 0; i < traits.length; i++) {
            assertEquals(format.format(expectedTraits[i]), format.format(traits[i]));
        }
        // Conditional moments (preorder)
        new TreeTipGradient("" +
                "trait", dataLikelihood, likelihoodDelegate, null);
        TreeTraitLogger treeTraitLogger = new TreeTraitLogger(treeModel,
                new TreeTrait[]{dataLikelihood.getTreeTrait("fcd.trait")},
                TreeTraitLogger.NodeRestriction.EXTERNAL, false);

        String moments = treeTraitLogger.getReport();
        double[] partials = parseVector(moments, "\t");
        testCMeans(s, "cMean ", partials);
        testCVariances(s, "cVar ", partials);
    }

    public void testLikelihoodDiagonalOURelaxed() {
        System.out.println("\nTest Likelihood using Diagonal OU Relaxed:");

        // Diffusion

        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        ArbitraryBranchRates.BranchRateTransform transform = make(false, false);
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

        assertEquals("likelihoodDiagonalOURelaxed",
                format.format(logDatumLikelihood),
                format.format(dataLikelihood.getLogLikelihood()));

        System.out.println("likelihoodDiagonalOURelaxed: " + format.format(logDatumLikelihood));

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihood, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{-1.0,2.0,0.0,1.811803424441062,0.6837595819961084,-1.0607909328094163,0.5,3.8623525502275142,5.5,2.0,5.0,-8.0,11.0,1.0,-1.5,1.0,2.5,4.0};
        double[] traits = parseVectorLine(treeTrait[0].getTraitString(treeModel, null), ",");
        for (int i = 0; i < traits.length; i++) {
            assertEquals(format.format(expectedTraits[i]), format.format(traits[i]));
        }
        // Conditional moments (preorder)
        new TreeTipGradient("" +
                "trait", dataLikelihood, likelihoodDelegate, null);
        TreeTraitLogger treeTraitLogger = new TreeTraitLogger(treeModel,
                new TreeTrait[]{dataLikelihood.getTreeTrait("fcd.trait")},
                TreeTraitLogger.NodeRestriction.EXTERNAL, false);

        String moments = treeTraitLogger.getReport();
        double[] partials = parseVector(moments, "\t");
        testCMeans(s, "cMean ", partials);
        testCVariances(s, "cVar ", partials);
    }

    public void testLikelihoodDiagonalOUBM() {
        System.out.println("\nTest Likelihood using Diagonal OU / BM:");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{1.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{2.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.3", new double[]{-2.0})));

        DiagonalMatrix strengthOfSelectionMatrixParam
                = new DiagonalMatrix(new Parameter.Default(new double[]{0.0, 0.0, 50.0}));

        DiffusionProcessDelegate diffusionProcessDelegate
                = new OUDiffusionModelDelegate(treeModel, diffusionModel,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

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

        assertEquals("likelihoodDiagonalOU",
                format.format(logDatumLikelihood),
                format.format(dataLikelihood.getLogLikelihood()));

        System.out.println("likelihoodDiagonalOU: " + format.format(logDatumLikelihood));

        // Conditional moments (preorder)
        new TreeTipGradient("" +
                "trait", dataLikelihood, likelihoodDelegate, null);
        TreeTraitLogger treeTraitLogger = new TreeTraitLogger(treeModel,
                new TreeTrait[]{dataLikelihood.getTreeTrait("fcd.trait")},
                TreeTraitLogger.NodeRestriction.EXTERNAL,
                false);

        String moments = treeTraitLogger.getReport();
        double[] partials = parseVector(moments, "\t");
        testCMeans(s, "cMean ", partials);
        testCVariances(s, "cVar ", partials);
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

        assertEquals("likelihoodFullOU",
                format.format(logDatumLikelihood),
                format.format(dataLikelihood.getLogLikelihood()));

        System.out.println("likelihoodFullOU: " + format.format(logDatumLikelihood));

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihood, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{-1.0,2.0,0.0,1.0427958776637916,2.060317467842193,0.5916377446549433,0.5,2.07249828895442,5.5,2.0,5.0,-8.0,11.0,1.0,-1.5,1.0,2.5,4.0};
        double[] traits = parseVectorLine(treeTrait[0].getTraitString(treeModel, null), ",");
        for (int i = 0; i < traits.length; i++) {
            assertEquals(format.format(expectedTraits[i]), format.format(traits[i]));
        }
        // Conditional moments (preorder)
        new TreeTipGradient("" +
                "trait", dataLikelihood, likelihoodDelegate, null);
        TreeTraitLogger treeTraitLogger = new TreeTraitLogger(treeModel,
                new TreeTrait[]{dataLikelihood.getTreeTrait("fcd.trait")},
                TreeTraitLogger.NodeRestriction.EXTERNAL, false);

        String moments = treeTraitLogger.getReport();
        double[] partials = parseVector(moments, "\t");
        testCMeans(s, "cMean ", partials);
        testCVariances(s, "cVar ", partials);
    }

    public void testLikelihoodFullOURelaxed() {
        System.out.println("\nTest Likelihood using Full OU Relaxed:");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        ArbitraryBranchRates.BranchRateTransform transform = make(false, false);
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

        assertEquals("likelihoodFullOURelaxed",
                format.format(logDatumLikelihood),
                format.format(dataLikelihood.getLogLikelihood()));

        System.out.println("likelihoodFullOURelaxed: " + format.format(logDatumLikelihood));

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihood, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{-1.0,2.0,0.0,1.6349449153945943,2.8676718538313635,-1.0653412418514505,0.5,3.3661883786009166,5.5,2.0,5.0,-8.0,11.0,1.0,-1.5,1.0,2.5,4.0};
        double[] traits = parseVectorLine(treeTrait[0].getTraitString(treeModel, null), ",");
        for (int i = 0; i < traits.length; i++) {
            assertEquals(format.format(expectedTraits[i]), format.format(traits[i]));
        }
        // Conditional moments (preorder)
        new TreeTipGradient("" +
                "trait", dataLikelihood, likelihoodDelegate, null);
        TreeTraitLogger treeTraitLogger = new TreeTraitLogger(treeModel,
                new TreeTrait[]{dataLikelihood.getTreeTrait("fcd.trait")},
                TreeTraitLogger.NodeRestriction.EXTERNAL, false);

        String moments = treeTraitLogger.getReport();
        double[] partials = parseVector(moments, "\t");
        testCMeans(s, "cMean ", partials);
        testCVariances(s, "cVar ", partials);
    }

    public void testLikelihoodFullAndDiagonalOU() {
        System.out.println("\nTest Likelihood comparing Full and Diagonal OU:");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        ArbitraryBranchRates.BranchRateTransform transform = make(false, false);
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

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

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

        assertEquals("likelihoodFullOURelaxed",
                format.format(logDatumLikelihood),
                format.format(dataLikelihood.getLogLikelihood()));

        System.out.println("likelihoodFullOURelaxed: " + format.format(logDatumLikelihood));

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

        assertEquals("likelihoodFullNonSymmetricOU",
                format.format(logDatumLikelihood),
                format.format(dataLikelihood.getLogLikelihood()));

        System.out.println("likelihoodFullNonSymmetricOU: " + format.format(logDatumLikelihood));

        // Conditional moments (preorder)
        new TreeTipGradient("" +
                "trait", dataLikelihood, likelihoodDelegate, null);
        TreeTraitLogger treeTraitLogger = new TreeTraitLogger(treeModel,
                new TreeTrait[]{dataLikelihood.getTreeTrait("fcd.trait")},
                TreeTraitLogger.NodeRestriction.EXTERNAL, false);

        String moments = treeTraitLogger.getReport();
        double[] partials = parseVector(moments, "\t");
        testCMeans(s, "cMean ", partials);
        testCVariances(s, "cVar ", partials);
    }

    public void testLikelihoodFullOUNonSymmetricRelaxed() {
        System.out.println("\nTest Likelihood using Full Non symmetric OU Relaxed:");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        ArbitraryBranchRates.BranchRateTransform transform = make(false, false);
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

        assertEquals("likelihoodFullNonSymmetricOURelaxed",
                format.format(logDatumLikelihood),
                format.format(dataLikelihood.getLogLikelihood()));

        System.out.println("likelihoodFullNonSymmetricOURelaxed: " + format.format(logDatumLikelihood));

        // Conditional moments (preorder)
        new TreeTipGradient("" +
                "trait", dataLikelihood, likelihoodDelegate, null);
        TreeTraitLogger treeTraitLogger = new TreeTraitLogger(treeModel,
                new TreeTrait[]{dataLikelihood.getTreeTrait("fcd.trait")},
                TreeTraitLogger.NodeRestriction.EXTERNAL, false);

        String moments = treeTraitLogger.getReport();
        double[] partials = parseVector(moments, "\t");
        testCMeans(s, "cMean ", partials);
        testCVariances(s, "cVar ", partials);
    }

    //// Factor Model //// *********************************************************************************************

    public void testLikelihoodBMFactor() {
        System.out.println("\nTest Likelihood using vanilla BM and factor:");

        // Diffusion
        DiffusionProcessDelegate diffusionProcessDelegate
                = new HomogeneousDiffusionModelDelegate(treeModel, diffusionModelFactor);

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateFactors = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelFactor, rootPriorFactor,
                rateTransformation, rateModel, true);

        dataModelFactor.setLikelihoodDelegate(likelihoodDelegateFactors);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodFactors
                = new TreeDataLikelihood(likelihoodDelegateFactors, treeModel, rateModel);

        String sf = dataModelFactor.getReport();
        int indLikBegF = sf.indexOf("logMultiVariateNormalDensity = ") + 31;
        int indLikEndF = sf.indexOf("\n", indLikBegF);
        char[] logDatumLikelihoodCharF = new char[indLikEndF - indLikBegF + 1];
        sf.getChars(indLikBegF, indLikEndF, logDatumLikelihoodCharF, 0);
        double logDatumLikelihoodFactor = Double.parseDouble(String.valueOf(logDatumLikelihoodCharF));

        double likelihoodFactorData = dataLikelihoodFactors.getLogLikelihood();
        double likelihoodFactorDiffusion = dataModelFactor.getLogLikelihood();

        assertEquals("likelihoodBMFactor",
                format.format(logDatumLikelihoodFactor),
                format.format(likelihoodFactorData + likelihoodFactorDiffusion));

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModelFactor, dataModelFactor, rootPrior, rateTransformation, likelihoodDelegateFactors);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihoodFactors, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{0.5465660348743779,1.4630588821162218,0.5269128635581201,1.5893136820798857,0.6710852429187679,1.5029828619880496,1.086143018033472,1.7199651678375532,0.8259954834667216,1.6198597117499574,0.8943759297423082,1.883525907495704};
        double[] traits = parseVectorLine(treeTrait[0].getTraitString(treeModel, null), ",");
        for (int i = 0; i < traits.length; i++) {
            assertEquals(format.format(expectedTraits[i]), format.format(traits[i]));
        }
    }

    public void testLikelihoodDriftFactor() {
        System.out.println("\nTest Likelihood using drifted BM and factor:");

        // Diffusion
        List<BranchRateModel> driftModels = new ArrayList<BranchRateModel>();
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{0.0})));
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{-40.0})));
        DiffusionProcessDelegate diffusionProcessDelegate
                = new DriftDiffusionModelDelegate(treeModel, diffusionModelFactor, driftModels);

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateFactors = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelFactor, rootPriorFactor,
                rateTransformation, rateModel, false);

        dataModelFactor.setLikelihoodDelegate(likelihoodDelegateFactors);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodFactors
                = new TreeDataLikelihood(likelihoodDelegateFactors, treeModel, rateModel);

        String sf = dataModelFactor.getReport();
        int indLikBegF = sf.indexOf("logMultiVariateNormalDensity = ") + 31;
        int indLikEndF = sf.indexOf("\n", indLikBegF);
        char[] logDatumLikelihoodCharF = new char[indLikEndF - indLikBegF + 1];
        sf.getChars(indLikBegF, indLikEndF, logDatumLikelihoodCharF, 0);
        double logDatumLikelihoodFactor = Double.parseDouble(String.valueOf(logDatumLikelihoodCharF));

        double likelihoodFactorData = dataLikelihoodFactors.getLogLikelihood();
        double likelihoodFactorDiffusion = dataModelFactor.getLogLikelihood();


        assertEquals("likelihoodDriftFactor",
                format.format(logDatumLikelihoodFactor),
                format.format(likelihoodFactorData + likelihoodFactorDiffusion));

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModelFactor, dataModelFactor, rootPrior, rateTransformation, likelihoodDelegateFactors);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihoodFactors, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{1.4253173030875967,-2.303529760054276,1.3911367775584844,-2.1804943569626927,1.5416862752420926,-2.2643516055441855,1.9616668293129644,-2.049007739844355,1.6865650403835397,-2.1618813046449183,1.7376526612221093,-1.9154773481169824};
        double[] traits = parseVectorLine(treeTrait[0].getTraitString(treeModel, null), ",");
        for (int i = 0; i < traits.length; i++) {
            assertEquals(format.format(expectedTraits[i]), format.format(traits[i]));
        }
    }

    public void testLikelihoodDriftRelaxedFactor() {
        System.out.println("\nTest Likelihood using drifted Relaxed BM and factor:");

        // Diffusion
        List<BranchRateModel> driftModels = new ArrayList<BranchRateModel>();
        ArbitraryBranchRates.BranchRateTransform transform = make(false, false);
        driftModels.add(new ArbitraryBranchRates(treeModel,
                new Parameter.Default("rate.1", new double[]{0, 10, 20, 30, 40, 40, 30, 20, 10, 0}),
                transform, false));
        driftModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{0.0})));
        DiffusionProcessDelegate diffusionProcessDelegate
                = new DriftDiffusionModelDelegate(treeModel, diffusionModelFactor, driftModels);

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateFactors = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelFactor, rootPriorFactor,
                rateTransformation, rateModel, false);

        dataModelFactor.setLikelihoodDelegate(likelihoodDelegateFactors);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodFactors
                = new TreeDataLikelihood(likelihoodDelegateFactors, treeModel, rateModel);

        String sf = dataModelFactor.getReport();
        int indLikBegF = sf.indexOf("logMultiVariateNormalDensity = ") + 31;
        int indLikEndF = sf.indexOf("\n", indLikBegF);
        char[] logDatumLikelihoodCharF = new char[indLikEndF - indLikBegF + 1];
        sf.getChars(indLikBegF, indLikEndF, logDatumLikelihoodCharF, 0);
        double logDatumLikelihoodFactor = Double.parseDouble(String.valueOf(logDatumLikelihoodCharF));

        double likelihoodFactorData = dataLikelihoodFactors.getLogLikelihood();
        double likelihoodFactorDiffusion = dataModelFactor.getLogLikelihood();


        assertEquals("likelihoodDriftRelaxedFactor",
                format.format(logDatumLikelihoodFactor),
                format.format(likelihoodFactorData + likelihoodFactorDiffusion));

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModelFactor, dataModelFactor, rootPrior, rateTransformation, likelihoodDelegateFactors);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihoodFactors, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{0.18706585597925554,1.3789549552109548,0.4154391209373806,1.4804852322163773,0.6694121384379099,1.3944259204124303,1.1456184528469746,1.5952491597077192,1.2346606015535597,1.3984078447275292,1.556929180473793,1.5691060274867998};
        double[] traits = parseVectorLine(treeTrait[0].getTraitString(treeModel, null), ",");
        for (int i = 0; i < traits.length; i++) {
            assertEquals(format.format(expectedTraits[i]), format.format(traits[i]));
        }
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

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateFactors = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelFactor, rootPriorFactor,
                rateTransformation, rateModel, false);

        dataModelFactor.setLikelihoodDelegate(likelihoodDelegateFactors);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodFactors
                = new TreeDataLikelihood(likelihoodDelegateFactors, treeModel, rateModel);

        String sf = dataModelFactor.getReport();
        int indLikBegF = sf.indexOf("logMultiVariateNormalDensity = ") + 31;
        int indLikEndF = sf.indexOf("\n", indLikBegF);
        char[] logDatumLikelihoodCharF = new char[indLikEndF - indLikBegF + 1];
        sf.getChars(indLikBegF, indLikEndF, logDatumLikelihoodCharF, 0);
        double logDatumLikelihoodFactor = Double.parseDouble(String.valueOf(logDatumLikelihoodCharF));

        double likelihoodFactorData = dataLikelihoodFactors.getLogLikelihood();
        double likelihoodFactorDiffusion = dataModelFactor.getLogLikelihood();


        assertEquals("likelihoodDiagonalOUFactor",
                format.format(logDatumLikelihoodFactor),
                format.format(likelihoodFactorData + likelihoodFactorDiffusion));

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModelFactor, dataModelFactor, rootPrior, rateTransformation, likelihoodDelegateFactors);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihoodFactors, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{1.2661812048403958,-1.5601987678298495,1.237446720479603,-1.4524772824103267,1.38571414295445,-1.534159429561456,1.7996426251091364,-1.4170923915969982,1.5205462963225242,-1.4540114665751946,1.6138006594873264,-1.4347696894136706};
        double[] traits = parseVectorLine(treeTrait[0].getTraitString(treeModel, null), ",");
        for (int i = 0; i < traits.length; i++) {
            assertEquals(format.format(expectedTraits[i]), format.format(traits[i]));
        }
    }

    public void testLikelihoodDiagonalOURelaxedFactor() {
        System.out.println("\nTest Likelihood using diagonal Relaxed OU and factor:");

        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        ArbitraryBranchRates.BranchRateTransform transform = make(false, false);
        optimalTraitsModels.add(new ArbitraryBranchRates(treeModel,
                new Parameter.Default("rate.1", new double[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}),
                transform, false));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{-1.5})));

        DiagonalMatrix strengthOfSelectionMatrixParam = new DiagonalMatrix(new Parameter.Default(new double[]{1.5, 20.0}));

        DiffusionProcessDelegate diffusionProcessDelegate
                = new OUDiffusionModelDelegate(treeModel, diffusionModelFactor,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateFactors = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelFactor, rootPriorFactor,
                rateTransformation, rateModel, false);

        dataModelFactor.setLikelihoodDelegate(likelihoodDelegateFactors);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodFactors
                = new TreeDataLikelihood(likelihoodDelegateFactors, treeModel, rateModel);

        String sf = dataModelFactor.getReport();
        int indLikBegF = sf.indexOf("logMultiVariateNormalDensity = ") + 31;
        int indLikEndF = sf.indexOf("\n", indLikBegF);
        char[] logDatumLikelihoodCharF = new char[indLikEndF - indLikBegF + 1];
        sf.getChars(indLikBegF, indLikEndF, logDatumLikelihoodCharF, 0);
        double logDatumLikelihoodFactor = Double.parseDouble(String.valueOf(logDatumLikelihoodCharF));

        double likelihoodFactorData = dataLikelihoodFactors.getLogLikelihood();
        double likelihoodFactorDiffusion = dataModelFactor.getLogLikelihood();


        assertEquals("likelihoodDiagonalOURelaxedFactor",
                format.format(logDatumLikelihoodFactor),
                format.format(likelihoodFactorData + likelihoodFactorDiffusion));

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModelFactor, dataModelFactor, rootPrior, rateTransformation, likelihoodDelegateFactors);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihoodFactors, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{1.195190236152443,-1.1789243642241476,1.299455400150337,-1.0647393453748792,1.4512482170336167,-1.1491156455961575,1.7445630013416964,-0.9846801795321203,1.4453824289220822,-1.0472280568477315,1.5386906072447606,-0.9850032177835296};
        double[] traits = parseVectorLine(treeTrait[0].getTraitString(treeModel, null), ",");
        for (int i = 0; i < traits.length; i++) {
            assertEquals(format.format(expectedTraits[i]), format.format(traits[i]));
        }
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

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateFactors = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelFactor, rootPriorFactor,
                rateTransformation, rateModel, false);

        dataModelFactor.setLikelihoodDelegate(likelihoodDelegateFactors);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodFactors
                = new TreeDataLikelihood(likelihoodDelegateFactors, treeModel, rateModel);

        String sf = dataModelFactor.getReport();
        int indLikBegF = sf.indexOf("logMultiVariateNormalDensity = ") + 31;
        int indLikEndF = sf.indexOf("\n", indLikBegF);
        char[] logDatumLikelihoodCharF = new char[indLikEndF - indLikBegF + 1];
        sf.getChars(indLikBegF, indLikEndF, logDatumLikelihoodCharF, 0);
        double logDatumLikelihoodFactor = Double.parseDouble(String.valueOf(logDatumLikelihoodCharF));

        double likelihoodFactorData = dataLikelihoodFactors.getLogLikelihood();
        double likelihoodFactorDiffusion = dataModelFactor.getLogLikelihood();


        assertEquals("likelihoodFullOUFactor",
                format.format(logDatumLikelihoodFactor),
                format.format(likelihoodFactorData + likelihoodFactorDiffusion));

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModelFactor, dataModelFactor, rootPrior, rateTransformation, likelihoodDelegateFactors);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihoodFactors, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{0.4653207352928792,1.8651780543656606,0.4504975010414803,1.9791256873553373,0.5928597377048458,1.8943915189612304,1.0019098658641454,2.0470937144421897,0.7411887507380247,1.9911336598314213,0.849042656244048,2.0371227193171015};
        double[] traits = parseVectorLine(treeTrait[0].getTraitString(treeModel, null), ",");
        for (int i = 0; i < traits.length; i++) {
            assertEquals(format.format(expectedTraits[i]), format.format(traits[i]));
        }
    }

    public void testLikelihoodFullOURelaxedFactor() {
        System.out.println("\nTest Likelihood using full Relaxed OU and factor:");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        ArbitraryBranchRates.BranchRateTransform transform = make(false, false);
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

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateFactors = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelFactor, rootPriorFactor,
                rateTransformation, rateModel, false);

        dataModelFactor.setLikelihoodDelegate(likelihoodDelegateFactors);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodFactors
                = new TreeDataLikelihood(likelihoodDelegateFactors, treeModel, rateModel);

        String sf = dataModelFactor.getReport();
        int indLikBegF = sf.indexOf("logMultiVariateNormalDensity = ") + 31;
        int indLikEndF = sf.indexOf("\n", indLikBegF);
        char[] logDatumLikelihoodCharF = new char[indLikEndF - indLikBegF + 1];
        sf.getChars(indLikBegF, indLikEndF, logDatumLikelihoodCharF, 0);
        double logDatumLikelihoodFactor = Double.parseDouble(String.valueOf(logDatumLikelihoodCharF));

        double likelihoodFactorData = dataLikelihoodFactors.getLogLikelihood();
        double likelihoodFactorDiffusion = dataModelFactor.getLogLikelihood();


        assertEquals("likelihoodFullRelaxedOUFactor",
                format.format(logDatumLikelihoodFactor),
                format.format(likelihoodFactorData + likelihoodFactorDiffusion));

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModelFactor, dataModelFactor, rootPrior, rateTransformation, likelihoodDelegateFactors);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihoodFactors, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{0.5771219814883327,1.4228543643506777,0.6061440030287131,1.54565271965551,0.7505232714938799,1.4621507423823044,1.1211158357858184,1.6081573398162914,0.849513804864352,1.5512812784344256,0.956280303218616,1.5995936766474121};
        double[] traits = parseVectorLine(treeTrait[0].getTraitString(treeModel, null), ",");
        for (int i = 0; i < traits.length; i++) {
            assertEquals(format.format(expectedTraits[i]), format.format(traits[i]));
        }
    }

    public void testLikelihoodFullDiagonalOUFactor() {
        System.out.println("\nTest Likelihood comparing full and diagonal OU and factor:");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        ArbitraryBranchRates.BranchRateTransform transform = make(false, false);
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

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

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

    private double[] parseVectorLine(String s, String sep) {
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

    private void testCMeans(String s, String name, double[] partials) {
        int offset = 0;
        int indBeg = 0;
        for (int tip = 0; tip < nTips; tip++) {
            indBeg = s.indexOf(name, indBeg + 1) + name.length() + 4;
            int indEnd = s.indexOf("]", indBeg);
            double[] vector = parseVector(s.substring(indBeg, indEnd), ",");
            for (int i = 0; i < vector.length; i++) {
//                System.out.println("cMean Mat: " + vector[i]);
//                System.out.println("cMean preorder: " + partials[offset + i]);
                assertEquals("cMean " + tip + "; " + i,
                        format.format(partials[offset + i]),
                        format.format(vector[i]));
            }
            offset += dimTrait + 2 * dimTrait * dimTrait + 1;
        }
    }

    private void testCVariances(String s, String name, double[] partials) {
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
            offset += dimTrait + 2 * dimTrait * dimTrait + 1;
        }
    }
}