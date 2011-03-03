package test.dr.calibration;

import dr.evolution.coalescent.CoalescentSimulator;
import dr.evolution.coalescent.ConstantPopulation;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.operators.SubtreeSlideOperator;
import dr.evomodel.speciation.BirthDeathGernhard08Model;
import dr.evomodel.speciation.SpeciationLikelihood;
import dr.evomodel.speciation.SpeciationModel;
import dr.evomodel.tree.TreeHeightStatistic;
import dr.evomodel.tree.TreeLengthStatistic;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.evomodelxml.tree.TreeModelParser;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.loggers.ArrayLogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
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

    public TestCalibratedYuleModel(int treeSize, double S, BufferedWriter out) throws IOException {
//        this.treeSize = treeSize;
        this.out = out;
        out.write(Integer.toString(treeSize) + "\t");

        TreeModel treeModel = createTreeModel(treeSize);

        Parameter brParameter = new Parameter.Default("birthRate", 2.0, 0.0, 100.0);

        OperatorSchedule schedule = new SimpleOperatorSchedule();
        MCMCOperator operator = new SubtreeSlideOperator(treeModel, 10, 1, true, false, false, false, CoercionMode.COERCION_ON);
        schedule.addOperator(operator);

        operator = new ScaleOperator(brParameter, 0.5);
        operator.setWeight(1.0);
        schedule.addOperator(operator);

        yuleTester(treeModel, schedule, brParameter, S);

    }


    protected TreeModel createTreeModel(int treeSize) {
        Taxa taxonList = new Taxa();
        for (int i = 0; i < treeSize; i++) {
            taxonList.addTaxon(new Taxon("T" + Integer.toString(i)));
        }

        Parameter popSize = new Parameter.Default(treeSize);
        popSize.setId(ConstantPopulationModelParser.POPULATION_SIZE);
        ConstantPopulationModel startingTree = new ConstantPopulationModel(popSize, Units.Type.YEARS);
        ConstantPopulation constant = (ConstantPopulation) startingTree.getDemographicFunction();


        CoalescentSimulator simulator = new CoalescentSimulator();
        Tree tree = simulator.simulateTree(taxonList, constant);

        return new TreeModel(tree);//treeModel
    }

    private void yuleTester(TreeModel treeModel, OperatorSchedule schedule, Parameter brParameter, double S) throws IOException {

        MCMC mcmc = new MCMC("mcmc1");
        MCMCOptions options = new MCMCOptions();
        options.setChainLength(20000000);
        options.setUseCoercion(true);
        options.setCoercionDelay(options.getChainLength() / 100);
        options.setTemperature(1.0);
        options.setFullEvaluationCount(2000);

        TreeLengthStatistic tls = new TreeLengthStatistic(TL, treeModel);
        TreeHeightStatistic rootHeight = new TreeHeightStatistic(TREE_HEIGHT, treeModel);

        SpeciationModel speciationModel = new BirthDeathGernhard08Model("yule", brParameter, null, null,
                BirthDeathGernhard08Model.TreeType.UNSCALED, Units.Type.SUBSTITUTIONS, false);
        Likelihood speciationLikelihood = new SpeciationLikelihood(treeModel, speciationModel, "yule.like");

        DistributionLikelihood logNormalLikelihood = new DistributionLikelihood(
                new LogNormalDistribution(1.0, S), 0); // meanInRealSpace="false"
        logNormalLikelihood.addData(treeModel.getRootHeightParameter());

        //CompoundLikelihood
        List<Likelihood> likelihoods = new ArrayList<Likelihood>();
        likelihoods.add(speciationLikelihood);
        likelihoods.add(logNormalLikelihood);
        Likelihood prior = new CompoundLikelihood(0, likelihoods);
        prior.setId(CompoundLikelihoodParser.PRIOR);

        ArrayLogFormatter logformatter = new ArrayLogFormatter(false);

        MCLogger[] loggers = new MCLogger[1];
        loggers[0] = new MCLogger(logformatter, options.getChainLength() / 10000, false);
        loggers[0].add(speciationLikelihood);
        loggers[0].add(rootHeight);
        loggers[0].add(tls);
        loggers[0].add(brParameter);

//        loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), options.getChainLength() / 100000, false);
//        loggers[1].add(speciationLikelihood);
//        loggers[1].add(rootHeight);
//        loggers[1].add(tls);

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
        TraceCorrelation treeHeightStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(TREE_HEIGHT));
        out.write(formatter.format(treeHeightStats.getMean()));
        out.write("\t");

        double expectedRootHeight = Math.pow(Math.E, (1 + (Math.pow(0.2, 2) / 2)));
        out.write(formatter.format(expectedRootHeight));
        out.write("\t");

        double error = Math.abs( (treeHeightStats.getMean() - expectedRootHeight) / expectedRootHeight );
        NumberFormat percentFormatter = NumberFormat.getNumberInstance();
        percentFormatter.setMinimumFractionDigits(5);
        percentFormatter.setMinimumFractionDigits(5);
        out.write(percentFormatter.format(error));
        out.write("\t");
        out.write(Double.toString(tlStats.getESS()));

        System.out.println("rootHeight = " + formatter.format(treeHeightStats.getMean())
                + ";  expectation = " + formatter.format(expectedRootHeight)
                + ";  error = " + percentFormatter.format(error)
                + ";  tl.ess = " + tlStats.getESS());

    }


    public static void main(String[] args) {
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter("TestCalibratedYuleModel.txt"));
            out.write("treeSize\trootHeight\texpectation\terror");
            out.newLine();

            int[] taxaSchedule = new int[] {4,6,8,12,16,24,32,48,64,96,128};
            double[] S_Schedule = new double[] {0.05, 0.1, 0.2, 0.4};

            for (double S : S_Schedule) {
            for (int i : taxaSchedule) {
                System.out.print("tree size = " + i + "\t" + "S=" + S + "\t");
                out.write("tree size = " + i + "\t" + "S=" + S + "\t");
                TestCalibratedYuleModel testCalibratedYuleModel = new TestCalibratedYuleModel(i, S, out);
                out.newLine();
            }
            }
            out.close();
        } catch (IOException e) {
        }
    }
}
