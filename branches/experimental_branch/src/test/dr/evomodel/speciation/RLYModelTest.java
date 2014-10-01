package test.dr.evomodel.speciation;

import dr.evolution.io.NewickImporter;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.util.Units;
import dr.evomodel.operators.TreeBitRandomWalkOperator;
import dr.evomodel.speciation.RandomLocalYuleModel;
import dr.evomodel.speciation.SpeciationLikelihood;
import dr.evomodel.speciation.SpeciationModel;
import dr.evomodel.tree.TreeHeightStatistic;
import dr.evomodel.tree.TreeLengthStatistic;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.OldCoalescentSimulatorParser;
import dr.inference.loggers.ArrayLogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.BitFlipOperator;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.inference.trace.ArrayTraceList;
import dr.inference.trace.Trace;
import dr.inference.trace.TraceCorrelation;
import dr.math.MathUtils;
import junit.framework.Test;
import junit.framework.TestSuite;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.util.List;

/**
 * YuleModel Tester.
 *
 * @author Alexei Drummond
 * @version 1.0
 * @since <pre>08/26/2007</pre>
 */
public class RLYModelTest extends TraceCorrelationAssert {

    static final String TL = "TL";
    static final String TREE_HEIGHT = OldCoalescentSimulatorParser.ROOT_HEIGHT;
    static final String birthRateIndicator = "birthRateIndicator";
    static final String birthRate = "birthRate";

    private FlexibleTree tree;

    public RLYModelTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        MathUtils.setSeed(666);

        NewickImporter importer = new NewickImporter(
                "(((((A:1.0,B:1.0):1.0,C:2.0),D:3.0):1.0, E:4.0),F:5.0);");
        tree = (FlexibleTree) importer.importTree(null);
    }

    public void testTreeBitRandomWalk() {

        TreeModel treeModel = new TreeModel("treeModel", tree);

        Parameter I = treeModel.createNodeTraitsParameter(
                birthRateIndicator, new double[]{1});

        Parameter b = treeModel.createNodeTraitsParameter(
                birthRate, new double[]{1});

        OperatorSchedule schedule = new SimpleOperatorSchedule();

        TreeBitRandomWalkOperator tbrw =
                new TreeBitRandomWalkOperator(treeModel, birthRateIndicator, birthRate, 1.0, 4, true);
        BitFlipOperator bfo = new BitFlipOperator(I, 1.0, true);

        schedule.addOperator(tbrw);
        schedule.addOperator(bfo);

        randomLocalYuleTester(treeModel, I, b, schedule);

    }

    private void randomLocalYuleTester(TreeModel treeModel, Parameter I, Parameter b, OperatorSchedule schedule) {

        MCMC mcmc = new MCMC("mcmc1");
        MCMCOptions options = new MCMCOptions(1000000);

        TreeLengthStatistic tls = new TreeLengthStatistic(TL, treeModel);
        TreeHeightStatistic rootHeight = new TreeHeightStatistic(TREE_HEIGHT, treeModel);

        Parameter m = new Parameter.Default("m", 1.0, 0.0, Double.MAX_VALUE);

        SpeciationModel speciationModel = new RandomLocalYuleModel(b, I, m, false, Units.Type.YEARS, 4);

        Likelihood likelihood = new SpeciationLikelihood(treeModel, speciationModel, "randomYule.like");

        ArrayLogFormatter formatter = new ArrayLogFormatter(false);

        MCLogger[] loggers = new MCLogger[2];
        loggers[0] = new MCLogger(formatter, 100, false);
        loggers[0].add(likelihood);
        loggers[0].add(rootHeight);
        loggers[0].add(tls);
        loggers[0].add(I);

        loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), 100000, false);
        loggers[1].add(likelihood);
        loggers[1].add(rootHeight);
        loggers[1].add(tls);
        loggers[1].add(I);

        mcmc.setShowOperatorAnalysis(true);

        mcmc.init(options, likelihood, schedule, loggers);

        mcmc.run();

        List<Trace> traces = formatter.getTraces();
        ArrayTraceList traceList = new ArrayTraceList("yuleModelTest", traces, 0);

        for (int i = 1; i < traces.size(); i++) {
            traceList.analyseTrace(i);
        }

        TraceCorrelation tlStats =
                traceList.getCorrelationStatistics(traceList.getTraceIndex("root." + birthRateIndicator));

        System.out.println("mean = " + tlStats.getMean());
        System.out.println("expected mean = 0.5");

        assertExpectation("root." + birthRateIndicator, tlStats, 0.5);
    }

    public static Test suite() {
        return new TestSuite(RLYModelTest.class);
    }
}
