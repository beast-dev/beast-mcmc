package test.dr.calibration;


import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.coalescent.DemographicModel;
import dr.evomodel.operators.SubtreeSlideOperator;
import dr.evomodel.speciation.BirthDeathGernhard08Model;
import dr.evomodel.speciation.SpeciationLikelihood;
import dr.evomodel.speciation.SpeciationModel;
import dr.evomodel.tree.*;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.evomodelxml.tree.TreeModelParser;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.loggers.ArrayLogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.BooleanLikelihood;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.inference.trace.ArrayTraceList;
import dr.inference.trace.Trace;
import dr.inference.trace.TraceCorrelation;
import dr.inferencexml.model.CompoundLikelihoodParser;
import dr.math.distributions.LogNormalDistribution;
import dr.util.NumberFormatter;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * @auther Walter Xie
 */

public class TestCalibratedYuleModel {
    protected static final String TL = "TL";
    protected static final String TREE_HEIGHT = TreeModel.TREE_MODEL + "." + TreeModelParser.ROOT_HEIGHT;
    //    private final int treeSize;
    private final BufferedWriter out;
    Taxa taxa;
    private static final double M = 1 + Math.log(10);

    public TestCalibratedYuleModel(int treeSize, double S, int chainLength, BufferedWriter out) throws Exception {
//        this.treeSize = treeSize;
        this.out = out;
//        out.write(Integer.toString(treeSize) + "\t");

        TreeModel treeModel = createTreeModel(treeSize);

        Parameter brParameter = new Parameter.Default("birthRate", 2.0, 0.0, 100.0);

        OperatorSchedule schedule = new SimpleOperatorSchedule();
        MCMCOperator operator = new SubtreeSlideOperator(treeModel, 10, 1, true, false, false, false, CoercionMode.COERCION_ON);
        schedule.addOperator(operator);

        operator = new ScaleOperator(brParameter, 0.5);
        operator.setWeight(1.0);
        schedule.addOperator(operator);

        System.out.println("treeModel = " + Tree.Utils.newickNoLengths(treeModel));
//        out.write("\t");
//        out.write("treeModel = \t");
        out.write(Tree.Utils.newickNoLengths(treeModel));
        out.write("\t");
        yuleTester(treeModel, schedule, brParameter, S, chainLength);

    }


    protected TreeModel createTreeModel(int treeSize) throws Exception {
        taxa = new Taxa();
        for (int i = 0; i < treeSize; i++) {
            taxa.addTaxon(new Taxon("T" + Integer.toString(i)));
        }

        //System.out.println("taxaSubSet_size = " + taxaSubSet.getTaxonCount());

        Parameter popSize = new Parameter.Default(treeSize);
        popSize.setId(ConstantPopulationModelParser.POPULATION_SIZE);
        ConstantPopulationModel startingTree = new ConstantPopulationModel(popSize, Units.Type.YEARS);

        Tree tree = calibration(taxa, startingTree);

        return new TreeModel(tree);//treeModel
    }

    private Tree calibration(final TaxonList taxa, DemographicModel demoModel) throws Exception {

        dr.evolution.coalescent.CoalescentSimulator simulator = new dr.evolution.coalescent.CoalescentSimulator();

        DemographicFunction demoFunction = demoModel.getDemographicFunction();

        SimpleNode[] firstHalfNodes = new SimpleNode[taxa.getTaxonCount() / 2];
        SimpleNode[] secondHalfNodes = new SimpleNode[taxa.getTaxonCount() - taxa.getTaxonCount() / 2];
        for (int i = 0; i < firstHalfNodes.length; i++) {
            firstHalfNodes[i] = new SimpleNode();
            firstHalfNodes[i].setTaxon(taxa.getTaxon(i));
        }
        for (int i = 0; i < secondHalfNodes.length; i++) {
            secondHalfNodes[i] = new SimpleNode();
            secondHalfNodes[i].setTaxon(taxa.getTaxon(i + taxa.getTaxonCount() / 2));
        }

        SimpleNode firstHalfRootNode = simulator.simulateCoalescent(firstHalfNodes, demoFunction);

        SimpleNode[] restNodes = simulator.simulateCoalescent(secondHalfNodes, demoFunction, 0, firstHalfRootNode.getHeight());

        SimpleNode[] together = new SimpleNode[restNodes.length + 1];
        for (int i = 0; i < restNodes.length; i++) {
            together[i + 1] = restNodes[i];
        }
        together[0] = firstHalfRootNode;

        SimpleNode root = simulator.simulateCoalescent(together, demoFunction);

        Tree tree = new SimpleTree(root);

        return tree;
    }

    private void yuleTester(TreeModel treeModel, OperatorSchedule schedule, Parameter brParameter, double S, int chainLength)
            throws IOException, Tree.MissingTaxonException {

        MCMC mcmc = new MCMC("mcmc1");
        MCMCOptions options = new MCMCOptions();
        options.setChainLength(chainLength);
        options.setUseCoercion(true);
        options.setCoercionDelay((int)(options.getChainLength() / 100));
        options.setTemperature(1.0);
        options.setFullEvaluationCount(2000);

        TreeLengthStatistic tls = new TreeLengthStatistic(TL, treeModel);
        TreeHeightStatistic rootHeight = new TreeHeightStatistic(TREE_HEIGHT, treeModel);

        SpeciationModel speciationModel = new BirthDeathGernhard08Model("yule", brParameter, null, null,
                BirthDeathGernhard08Model.TreeType.UNSCALED, Units.Type.SUBSTITUTIONS, false);
        Likelihood speciationLikelihood = new SpeciationLikelihood(treeModel, speciationModel, "yule.like");

        Taxa halfTaxa = new Taxa();
        for (int i = 0; i < taxa.getTaxonCount() / 2; i++) {
            halfTaxa.addTaxon(new Taxon("T" + Integer.toString(i)));
        }

        TMRCAStatistic tmrca = new TMRCAStatistic("tmrca(halfTaxa)", treeModel, halfTaxa, false, false);
        DistributionLikelihood logNormalLikelihood = new DistributionLikelihood(
                new LogNormalDistribution(M, S), 0); // meanInRealSpace="false"
        logNormalLikelihood.addData(tmrca);

        MonophylyStatistic monophylyStatistic = new MonophylyStatistic("monophyly(halfTaxa)", treeModel, halfTaxa, null);
        BooleanLikelihood booleanLikelihood = new BooleanLikelihood();
        booleanLikelihood.addData(monophylyStatistic);

        //CompoundLikelihood
        List<Likelihood> likelihoods = new ArrayList<Likelihood>();
        likelihoods.add(speciationLikelihood);
        likelihoods.add(logNormalLikelihood);
        likelihoods.add(booleanLikelihood);
        Likelihood prior = new CompoundLikelihood(0, likelihoods);
        prior.setId(CompoundLikelihoodParser.PRIOR);

        ArrayLogFormatter logformatter = new ArrayLogFormatter(false);

        MCLogger[] loggers = new MCLogger[1];
        loggers[0] = new MCLogger(logformatter, (int)(options.getChainLength() / 10000), false);
        loggers[0].add(speciationLikelihood);
        loggers[0].add(rootHeight);
        loggers[0].add(tmrca);
        loggers[0].add(tls);
        loggers[0].add(brParameter);

        mcmc.setShowOperatorAnalysis(false);

        mcmc.init(options, prior, schedule, loggers);
        mcmc.run();

        List<Trace> traces = logformatter.getTraces();
        ArrayTraceList traceList = new ArrayTraceList("yuleModelTest", traces, 1000);

        for (int i = 1; i < traces.size(); i++) {
            traceList.analyseTrace(i);
        }

        NumberFormatter formatter = new NumberFormatter(8);

        TraceCorrelation tlStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(TL));
        TraceCorrelation treeHeightStats = traceList.getCorrelationStatistics(traceList.getTraceIndex("tmrca(halfTaxa)"));
//        out.write("tmrcaHeight = \t");
        out.write(formatter.format(treeHeightStats.getMean()));
        out.write("\t");

        double expectedNodeHeight = Math.pow(Math.E, (M + (Math.pow(S, 2) / 2)));
//        out.write("expectation = \t");
        out.write(formatter.format(expectedNodeHeight));
        out.write("\t");

        double error = Math.abs((treeHeightStats.getMean() - expectedNodeHeight) / expectedNodeHeight);
        NumberFormat percentFormatter = NumberFormat.getNumberInstance();
        percentFormatter.setMinimumFractionDigits(5);
        percentFormatter.setMinimumFractionDigits(5);
//        out.write("error = \t");
        out.write(percentFormatter.format(error));
        out.write("\t");
//        out.write("tl.ess = \t");
        out.write(Double.toString(tlStats.getESS()));

        System.out.println("tmrcaHeight = " + formatter.format(treeHeightStats.getMean())
                + ";  expectation = " + formatter.format(expectedNodeHeight)
                + ";  error = " + percentFormatter.format(error)
                + ";  tl.ess = " + tlStats.getESS());

    }


    public static void main(String[] args) {
        try {
            System.out.println("M = " + M);
            BufferedWriter out = new BufferedWriter(new FileWriter("TestCalibratedYuleModel.txt"));
            out.write("M = \t" + M);
            out.newLine();
            out.write("chainLeng\ttreeSize\tS\ttreeModel\ttmrcaHeight\texpectation\terror\ttl.ess");
            out.newLine();

            int[] taxaSchedule = new int[]{4, 8, 16, 32, 48, 64};
            int[] chainLengthSchedule = new int[]{2000000, 2000000, 4000000, 8000000, 10000000, 10000000};

            if (taxaSchedule.length != chainLengthSchedule.length) throw new Exception ();

            double[] S_Schedule = new double[]{0.05, 0.1, 0.2, 0.4};

            for (double S : S_Schedule) {

                for (int i = 0; i < taxaSchedule.length; i++) {

                    int n = taxaSchedule[i];
                    int chainLength = chainLengthSchedule[i];

                    System.out.print("chainLeng = " + chainLengthSchedule[i] + "\tS = " + S + "\ttreeSize = " + n + "\t");
                    out.write(chainLengthSchedule[i] + "\t" + S + "\t" + n + "\t");
                    TestCalibratedYuleModel testCalibratedYuleModel = new TestCalibratedYuleModel(n, S, chainLength, out);
                    out.newLine();
                }
            }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
