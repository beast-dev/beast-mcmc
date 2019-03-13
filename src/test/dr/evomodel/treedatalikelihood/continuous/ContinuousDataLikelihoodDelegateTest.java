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
import dr.math.matrixAlgebra.Vector;
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
    private ConjugateRootTraitPrior rootPriorInf;
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

        // Root prior Inf
        String sInf = "<beast>\n" +
                "    <conjugateRootPrior>\n" +
                "        <meanParameter>\n" +
                "            <parameter id=\"meanRoot\"  value=\"-1.0 -3.0 -2.0\"/>\n" +
                "        </meanParameter>\n" +
                "        <priorSampleSize>\n" +
                "            <parameter id=\"sampleSizeRoot\" value=\"Infinity\"/>\n" +
                "        </priorSampleSize>\n" +
                "    </conjugateRootPrior>\n" +
                "</beast>";
        XMLParser parserInf = new XMLParser(true, true, true, null);
        parserInf.addXMLObjectParser(new AttributeParser());
        parserInf.addXMLObjectParser(new ParameterParser());
        parserInf.parse(new StringReader(sInf), true);
        rootPriorInf = ConjugateRootTraitPrior.parseConjugateRootTraitPrior(parserInf.getRoot(), dimTrait);

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
        Parameter factorPrecisionParameters = new Parameter.Default("factorPrecision", new double[]{1.0, 5.0, 0.5});

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

        testLikelihood("likelihoodBM", dataLikelihood);

        // Conditional moments (preorder)
        testConditionalMoments(dataLikelihood, likelihoodDelegate);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihood, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{-1.0, 2.0, 0.0, 0.45807521679597646, 2.6505355982097605, 3.4693334367360538, 0.5, 2.64206285585883, 5.5, 2.0, 5.0, -8.0, 11.0, 1.0, -1.5, 1.0, 2.5, 4.0};
        double[] traits = parseVectorLine(treeTrait[0].getTraitString(treeModel, null), ",");
        System.err.println(new Vector(traits));
        for (int i = 0; i < traits.length; i++) {
            assertEquals(format.format(expectedTraits[i]), format.format(traits[i]));
        }

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

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

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
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihood, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{-1.0, 2.0, 0.0, 0.5457621072639138, 3.28662834718796, 3.2939596558001845, 0.5, 1.0742799493604265, 5.5, 2.0, 5.0, -8.0, 11.0, 1.0, -1.5, 1.0, 2.5, 4.0};
        double[] traits = parseVectorLine(treeTrait[0].getTraitString(treeModel, null), ",");
        for (int i = 0; i < traits.length; i++) {
            assertEquals(format.format(expectedTraits[i]), format.format(traits[i]));
        }

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

        testLikelihood("likelihoodDriftRelaxed", dataLikelihood);

        // Conditional moments (preorder)
        testConditionalMoments(dataLikelihood, likelihoodDelegate);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihood, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{-1.0, 2.0, 0.0, 2.843948876154644, 10.866053719140933, 3.467579698926694, 0.5, 12.000214659757933, 5.5, 2.0, 5.0, -8.0, 11.0, 1.0, -1.5, 1.0, 2.5, 4.0};
        double[] traits = parseVectorLine(treeTrait[0].getTraitString(treeModel, null), ",");
        for (int i = 0; i < traits.length; i++) {
            assertEquals(format.format(expectedTraits[i]), format.format(traits[i]));
        }

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

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

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
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihood, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{-1.0, 2.0, 0.0, 1.0369622398437415, 2.065450266793184, 0.6174755164694558, 0.5, 2.0829935706195615, 5.5, 2.0, 5.0, -8.0, 11.0, 1.0, -1.5, 1.0, 2.5, 4.0};
        double[] traits = parseVectorLine(treeTrait[0].getTraitString(treeModel, null), ",");
        for (int i = 0; i < traits.length; i++) {
            assertEquals(format.format(expectedTraits[i]), format.format(traits[i]));
        }

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

        testLikelihood("likelihoodDiagonalOURelaxed", dataLikelihood);

        // Conditional moments (preorder)
        testConditionalMoments(dataLikelihood, likelihoodDelegate);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihood, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{-1.0, 2.0, 0.0, 1.811803424441062, 0.6837595819961084, -1.0607909328094163, 0.5, 3.8623525502275142, 5.5, 2.0, 5.0, -8.0, 11.0, 1.0, -1.5, 1.0, 2.5, 4.0};
        double[] traits = parseVectorLine(treeTrait[0].getTraitString(treeModel, null), ",");
        for (int i = 0; i < traits.length; i++) {
            assertEquals(format.format(expectedTraits[i]), format.format(traits[i]));
        }

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

        testLikelihood("likelihoodFullOU", dataLikelihood);

        // Conditional moments (preorder)
        testConditionalMoments(dataLikelihood, likelihoodDelegate);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihood, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{-1.0, 2.0, 0.0, 1.0427958776637916, 2.060317467842193, 0.5916377446549433, 0.5, 2.07249828895442, 5.5, 2.0, 5.0, -8.0, 11.0, 1.0, -1.5, 1.0, 2.5, 4.0};
        double[] traits = parseVectorLine(treeTrait[0].getTraitString(treeModel, null), ",");
        for (int i = 0; i < traits.length; i++) {
            assertEquals(format.format(expectedTraits[i]), format.format(traits[i]));
        }

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

        testLikelihood("likelihoodFullOURelaxed", dataLikelihood);

        // Conditional moments (preorder)
        testConditionalMoments(dataLikelihood, likelihoodDelegate);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihood, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{-1.0, 2.0, 0.0, 1.6349449153945943, 2.8676718538313635, -1.0653412418514505, 0.5, 3.3661883786009166, 5.5, 2.0, 5.0, -8.0, 11.0, 1.0, -1.5, 1.0, 2.5, 4.0};
        double[] traits = parseVectorLine(treeTrait[0].getTraitString(treeModel, null), ",");
        for (int i = 0; i < traits.length; i++) {
            assertEquals(format.format(expectedTraits[i]), format.format(traits[i]));
        }

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

        testLikelihood("likelihoodFullIOU", dataLikelihood);

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

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

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

        testLikelihood("likelihoodBMFactor", dataModelFactor, dataLikelihoodFactors);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModelFactor, dataModelFactor, rootPrior, rateTransformation, likelihoodDelegateFactors);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihoodFactors, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{0.6002879987080073, 1.3630884580519484, 0.5250449300511655, 1.4853676908300644, 0.6673202215955497, 1.399820047380221, 1.0853554355129353, 1.6054879123935393, 0.4495494080256063, 1.4427296475118248, 0.8750789069500045, 1.8099596179292183};
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

        testLikelihood("likelihoodDriftFactor", dataModelFactor, dataLikelihoodFactors);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModelFactor, dataModelFactor, rootPrior, rateTransformation, likelihoodDelegateFactors);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihoodFactors, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{1.5058510863259034, -2.344107747791032, 1.415239714927795, -2.225937980916329, 1.5639840062954773, -2.3082612693286513, 1.9875205911751028, -2.1049011248405525, 1.3355460225282372, -2.2848471441564056, 1.742347318026791, -1.940903337116235};
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

        testLikelihood("likelihoodDriftRelaxedFactor", dataModelFactor, dataLikelihoodFactors);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModelFactor, dataModelFactor, rootPrior, rateTransformation, likelihoodDelegateFactors);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihoodFactors, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{0.21992781609528125, 1.271388273711557, 0.40761548539751596, 1.3682648770877144, 0.6599021787120436, 1.2830636141108613, 1.1488658943588324, 1.472103688153391, 0.8971632986744889, 1.20748933414854, 1.603739823726808, 1.4761482401796842};
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

        testLikelihood("likelihoodDiagonalOUFactor", dataModelFactor, dataLikelihoodFactors);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModelFactor, dataModelFactor, rootPrior, rateTransformation, likelihoodDelegateFactors);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihoodFactors, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{1.3270345780274333, -1.5589839744569975, 1.241407854756886, -1.4525648723106128, 1.388017192005544, -1.533399261149814, 1.8040948421311085, -1.4189758121385794, 1.1408165195832969, -1.4607180451268982, 1.6048925583434688, -1.4333922414628846};
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

        testLikelihood("likelihoodDiagonalOURelaxedFactor", dataModelFactor, dataLikelihoodFactors);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModelFactor, dataModelFactor, rootPrior, rateTransformation, likelihoodDelegateFactors);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihoodFactors, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{1.2546097113922914, -1.1761389606670978, 1.305611773283861, -1.0644815941127401, 1.4571577864569687, -1.1477885449972944, 1.749551506462585, -0.9890375857170963, 1.0763987351136657, -1.0671848958534547, 1.5276137550128892, -0.9822950795368887};
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

        testLikelihood("likelihoodFullOUFactor", dataModelFactor, dataLikelihoodFactors);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModelFactor, dataModelFactor, rootPrior, rateTransformation, likelihoodDelegateFactors);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihoodFactors, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{0.4889263054598222, 1.866143125522109, 0.41845209107775877, 1.978457443711536, 0.5589398189015322, 1.8942177991552116, 0.9699471556784252, 2.0423474270630155, 0.3288819110219145, 1.9759942582707206, 0.8081782260054755, 2.038299849681893};
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

        testLikelihood("likelihoodFullRelaxedOUFactor", dataModelFactor, dataLikelihoodFactors);

        // Conditional simulations
        MathUtils.setSeed(17890826);
        ProcessSimulationDelegate simulationDelegate =
                new MultivariateConditionalOnTipsRealizedDelegate("dataModel", treeModel,
                        diffusionModelFactor, dataModelFactor, rootPrior, rateTransformation, likelihoodDelegateFactors);
        ProcessSimulation simulationProcess = new ProcessSimulation(dataLikelihoodFactors, simulationDelegate);
        simulationProcess.cacheSimulatedTraits(null);
        TreeTrait[] treeTrait = simulationProcess.getTreeTraits();

        double[] expectedTraits = new double[]{0.6074917696668031, 1.4240248941610945, 0.5818653246406664, 1.545237778993696, 0.7248840308905077, 1.4623057820376757, 1.0961030597302799, 1.603694717986661, 0.44280937767720896, 1.5374906898020686, 0.920698984735896, 1.6011019734876784};
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
            offset += dimTrait + 2 * dimTrait * dimTrait + 1;
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
            offset += dimTrait + 2 * dimTrait * dimTrait + 1;
        }
    }

    private double getLogDatumLikelihood(TreeDataLikelihood dataLikelihood) {
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

    private void testLikelihood(String message, IntegratedFactorAnalysisLikelihood dataModelFactor, TreeDataLikelihood dataLikelihoodFactors) {
        String sf = dataModelFactor.getReport();
        int indLikBegF = sf.indexOf("logMultiVariateNormalDensity = ") + 31;
        int indLikEndF = sf.indexOf("\n", indLikBegF);
        char[] logDatumLikelihoodCharF = new char[indLikEndF - indLikBegF + 1];
        sf.getChars(indLikBegF, indLikEndF, logDatumLikelihoodCharF, 0);
        double logDatumLikelihoodFactor = Double.parseDouble(String.valueOf(logDatumLikelihoodCharF));

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
}