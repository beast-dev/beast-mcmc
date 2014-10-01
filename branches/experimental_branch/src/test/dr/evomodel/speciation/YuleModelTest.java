package test.dr.evomodel.speciation;

import dr.evolution.io.NewickImporter;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.util.Units;
import dr.evomodel.operators.SubtreeSlideOperator;
import dr.evomodel.speciation.BirthDeathGernhard08Model;
import dr.evomodel.speciation.SpeciationLikelihood;
import dr.evomodel.speciation.SpeciationModel;
import dr.evomodel.tree.TreeHeightStatistic;
import dr.evomodel.tree.TreeLengthStatistic;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.tree.TreeModelParser;
import dr.inference.loggers.ArrayLogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
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
public class YuleModelTest extends TraceCorrelationAssert {

    static final String TL = "TL";
    static final String TREE_HEIGHT = TreeModelParser.ROOT_HEIGHT;

    private FlexibleTree tree;

    public YuleModelTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        MathUtils.setSeed(666);
        
        NewickImporter importer = new NewickImporter("((1:1.0,2:1.0):1.0,(3:1.0,4:1.0):1.0);");
        tree = (FlexibleTree) importer.importTree(null);
    }

    public void testYuleWithSubtreeSlide() {

        TreeModel treeModel = new TreeModel("treeModel", tree);

        OperatorSchedule schedule = new SimpleOperatorSchedule();
        MCMCOperator operator =
                new SubtreeSlideOperator(treeModel, 1, 1, true, false, false, false, CoercionMode.COERCION_ON);
        schedule.addOperator(operator);

        yuleTester(treeModel, schedule);

    }

//    public void testYuleWithWideExchange() {
//
//        TreeModel treeModel = new TreeModel("treeModel", tree);

        // Doesn't compile...
  //      yuleTester(treeModel, ExchangeOperatorTest.getWideExchangeSchedule(treeModel));
//    }

    private void yuleTester(TreeModel treeModel, OperatorSchedule schedule) {

        MCMC mcmc = new MCMC("mcmc1");
        MCMCOptions options = new MCMCOptions(1000000);

        TreeLengthStatistic tls = new TreeLengthStatistic(TL, treeModel);
        TreeHeightStatistic rootHeight = new TreeHeightStatistic(TREE_HEIGHT, treeModel);

        Parameter b = new Parameter.Default("b", 2.0, 0.0, Double.MAX_VALUE);
        Parameter d = new Parameter.Default("d", 0.0, 0.0, Double.MAX_VALUE);

        SpeciationModel speciationModel = new BirthDeathGernhard08Model(b, d, null, BirthDeathGernhard08Model.TreeType.TIMESONLY,
                Units.Type.YEARS);
        Likelihood likelihood = new SpeciationLikelihood(treeModel, speciationModel, "yule.like");

        ArrayLogFormatter formatter = new ArrayLogFormatter(false);

        MCLogger[] loggers = new MCLogger[2];
        loggers[0] = new MCLogger(formatter, 100, false);
        loggers[0].add(likelihood);
        loggers[0].add(rootHeight);
        loggers[0].add(tls);

        loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), 100000, false);
        loggers[1].add(likelihood);
        loggers[1].add(rootHeight);
        loggers[1].add(tls);

        mcmc.setShowOperatorAnalysis(true);

        mcmc.init(options, likelihood, schedule, loggers);

        mcmc.run();

        List<Trace> traces = formatter.getTraces();
        ArrayTraceList traceList = new ArrayTraceList("yuleModelTest", traces, 0);

        for (int i = 1; i < traces.size(); i++) {
            traceList.analyseTrace(i);
        }

        // expectation of root height for 4 tips and lambda = 2
        // rootHeight = 0.541666
        // TL = 1.5

        TraceCorrelation tlStats =
                traceList.getCorrelationStatistics(traceList.getTraceIndex(TL));

        assertExpectation(TL, tlStats, 1.5);

        TraceCorrelation treeHeightStats =
                traceList.getCorrelationStatistics(traceList.getTraceIndex(TREE_HEIGHT));

        assertExpectation(TREE_HEIGHT, treeHeightStats, 0.5416666);


    }

    public static Test suite() {
        return new TestSuite(YuleModelTest.class);
    }
}
