package test.dr.evomodel.arg.operators;

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evomodel.arg.ARGLogger;
import dr.evomodel.arg.ARGModel;
import dr.evomodel.arg.ARGReassortmentNodeCountStatistic;
import dr.evomodel.arg.coalescent.ARGUniformPrior;
import dr.evomodelxml.tree.TreeModelParser;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.loggers.ArrayLogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.ScaleOperator;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.inference.trace.ArrayTraceList;
import dr.inference.trace.Trace;
import dr.inference.trace.TraceCorrelation;
import dr.math.distributions.GammaDistribution;
import dr.math.distributions.PoissonDistribution;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc Suchard
 */

public class ARGAddRemoveOperatorTest extends TraceCorrelationAssert {

    public ARGAddRemoveOperatorTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        NewickImporter importer = new NewickImporter(
                "(((A:1.0,B:1.0):1.0,C:2.0):1.0,D:3.0);");
        arg4 = new ARGModel(importer.importTree(null));
        arg4.setupHeightBounds();
        arg4.addLikelihoodCalculator(null);
        arg4.addLikelihoodCalculator(null);

        importer = new NewickImporter(
                "((((A:1.0,B:1.0):1.0,C:2.0):1.0,D:3.0):1.0, E:4.0);");
        ARGModel arg5 = new ARGModel(importer.importTree(null));

        importer = new NewickImporter(
                "(((((A:1.0,B:1.0):1.0,C:2.0):1.0,D:3.0):1.0, E:4.0),F:5.0);");
        ARGModel arg6 = new ARGModel(importer.importTree(null));
    }

    // 4 taxa args
    public void testFlatPrior4() throws IOException, Importer.ImportException {
        flatPriorTester(arg4, 2000000, 1000, 2.0, 100.0, 0.5, 3);
    }

    private void flatPriorTester(ARGModel arg, int chainLength, int sampleTreeEvery,
                                 double nodeCountSetting, double rootHeightAlpha, double rootHeightBeta, int maxCount)
            throws IOException, Importer.ImportException {

        MCMC mcmc = new MCMC("mcmc1");
        MCMCOptions options = new MCMCOptions(chainLength);

//        double nodeCountSetting = 2.0;
//        double rootHeightAlpha = 100;
//        double rootHeightBeta = 0.5;

        OperatorSchedule schedule = getSchedule(arg);

        ARGUniformPrior uniformPrior = new ARGUniformPrior(arg, maxCount, arg.getExternalNodeCount());

        PoissonDistribution poisson = new PoissonDistribution(nodeCountSetting);

        DistributionLikelihood nodeCountPrior = new DistributionLikelihood(poisson, 0.0);
        ARGReassortmentNodeCountStatistic nodeCountStatistic = new ARGReassortmentNodeCountStatistic("nodeCount", arg);
        nodeCountPrior.addData(nodeCountStatistic);

        DistributionLikelihood rootPrior = new DistributionLikelihood(new GammaDistribution(rootHeightAlpha, rootHeightBeta), 0.0);
        CompoundParameter rootHeight = (CompoundParameter) arg.createNodeHeightsParameter(true, false, false);
        rootPrior.addData(rootHeight);

        List<Likelihood> likelihoods = new ArrayList<Likelihood>();
        likelihoods.add(uniformPrior);
        likelihoods.add(rootPrior);
        likelihoods.add(nodeCountPrior);
        
        CompoundLikelihood compoundLikelihood = new CompoundLikelihood(1, likelihoods);
        compoundLikelihood.setId("likelihood1");

        MCLogger[] loggers = new MCLogger[3];

        loggers[0] = new MCLogger(new TabDelimitedFormatter(System.out), 10000, false);
        loggers[0].add(compoundLikelihood);
        loggers[0].add(arg);

        File file = new File("test.args");
        file.deleteOnExit();
        FileOutputStream out = new FileOutputStream(file);

        loggers[1] = new ARGLogger(arg, new TabDelimitedFormatter(out), sampleTreeEvery, "test");

        ArrayLogFormatter formatter = new ArrayLogFormatter(false);

        loggers[2] = new MCLogger(formatter, sampleTreeEvery, false);
        loggers[2].add(arg);
        arg.getRootHeightParameter().setId("root");
        loggers[2].add(arg.getRootHeightParameter());

        mcmc.setShowOperatorAnalysis(true);

        mcmc.init(options, compoundLikelihood, schedule, loggers);

        mcmc.run();
        out.flush();
        out.close();

        List<Trace> traces = formatter.getTraces();

//        Set<String> uniqueTrees = new HashSet<String>();
//
//        NexusImporter importer = new NexusImporter(new FileReader(file));
//        while (importer.hasTree()) {
//            Tree t = importer.importNextTree();
//            uniqueTrees.add(Tree.Utils.uniqueNewick(t, t.getRoot()));
//        }
//
//        TestCase.assertEquals(numTopologies, uniqueTrees.size());            List<Trace> traces = formatter.getTraces();

        ArrayTraceList traceList = new ArrayTraceList("ARGTest", traces, 0);

        for (int i = 1; i < traces.size(); i++) {
            traceList.analyseTrace(i);
        }

        TraceCorrelation nodeCountStats = traceList.getCorrelationStatistics(1);
        TraceCorrelation rootHeightStats = traceList.getCorrelationStatistics(4);

        assertExpectation("nodeCount", nodeCountStats, poisson.truncatedMean(maxCount));
        assertExpectation(TreeModelParser.ROOT_HEIGHT, rootHeightStats, rootHeightAlpha * rootHeightBeta);

    }

    public static OperatorSchedule getSchedule(ARGModel arg) {

        CompoundParameter rootHeight = (CompoundParameter) arg.createNodeHeightsParameter(true, false, false);
        CompoundParameter internalHeights = (CompoundParameter) arg.createNodeHeightsParameter(false, true, false);
        //CompoundParameter allInternalNodeHeights = (CompoundParameter) arg.createNodeHeightsParameter(true, true, false);
//        CompoundParameter rates = (CompoundParameter) arg.createNodeRatesParameter(false, true, true);

//        ARGAddRemoveEventOperator operator1 = new ARGAddRemoveEventOperator(arg, 5, 0.5,
//                CoercionMode.COERCION_ON, internalHeights, allInternalNodeHeights, rates, 0.9, null,-1);

        ScaleOperator operator2 = new ScaleOperator(rootHeight, 0.75, CoercionMode.COERCION_ON, 5);
        ScaleOperator operator3 = new ScaleOperator(internalHeights, 0.75, CoercionMode.COERCION_ON, 10);

        OperatorSchedule schedule = new SimpleOperatorSchedule();
//        schedule.addOperator(operator1);
        schedule.addOperator(operator2);
        schedule.addOperator(operator3);

        return schedule;

//		<scaleOperator id="rootOperator" scaleFactor="0.5"
//			weight="10">
//			<parameter idref="argModel.rootHeight" />
//		</scaleOperator>
//
//		<scaleOperator scaleFactor="0.95" weight="10">
//			<parameter idref="argModel.internalNodeHeights" />
//		</scaleOperator>

//            <ARGEventOperator weight="5" addProbability="0.5"
//			autoOptimize="false">
//			<argTreeModel idref="argModel" />
//			<internalNodes>
//				<parameter idref="argModel.internalNodeHeights" />
//			</internalNodes>
//			<internalNodesPlusRoot>
//				<parameter idref="argModel.allInternalNodeHeights" />
//			</internalNodesPlusRoot>
//			<nodeRates>
//				<parameter idref="argModel.rates" />
//			</nodeRates>
//		</ARGEventOperator>
    }


    private ARGModel arg4;
}
