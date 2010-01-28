package test.dr.integration;

import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.operators.SubtreeSlideOperator;
import dr.evomodel.operators.WilsonBalding;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HKY;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.evomodelxml.HKYParser;
import dr.inference.loggers.ArrayLogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.inference.prior.Prior;
import dr.inference.trace.ArrayTraceList;
import dr.inference.trace.Trace;
import dr.inference.trace.TraceCorrelation;
import dr.math.MathUtils;
import junit.framework.Test;
import junit.framework.TestSuite;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.util.List;

/**
 * @author Walter Xie
 * convert testLikelihood.xml in the folder /example
 */

public class LikelihoodTest extends TraceCorrelationAssert {

    private TreeModel treeModel;

    public LikelihoodTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        MathUtils.setSeed(666);

        createAlignment(HOMINID_TAXON_SEQUENCE, Nucleotides.INSTANCE);

    }

    public void testNewickTree() {

        SimpleNode[] nodes = new SimpleNode[10];
        for (int n=0; n < 10; n++) {
            nodes[n] = new SimpleNode();
        }

//        nodes[0].setHeight(0);
        nodes[0].setTaxon(taxa[0]); // human

        nodes[1].setTaxon(taxa[1]); // chimp

        nodes[2].setTaxon(taxa[2]); // bonobo

        nodes[3].setHeight(0.010772);
        nodes[3].addChild(nodes[1]);
        nodes[3].addChild(nodes[2]);

        nodes[4].setHeight(0.024003);
        nodes[4].addChild(nodes[0]);
        nodes[4].addChild(nodes[3]);

        nodes[5].setTaxon(taxa[3]); // gorilla

        nodes[6].setHeight(0.036038);
        nodes[6].addChild(nodes[4]);
        nodes[6].addChild(nodes[5]);

        nodes[7].setTaxon(taxa[4]); // orangutan

        nodes[8].setHeight(0.069125);
        nodes[8].addChild(nodes[6]);
        nodes[8].addChild(nodes[7]);

        nodes[9].setTaxon(taxa[5]); // siamang

        SimpleNode root = new SimpleNode();
        root.setHeight(0.099582);
        root.addChild(nodes[8]);
        root.addChild(nodes[9]);

        Tree tree = new SimpleTree(root);
        tree.setUnits(Units.Type.YEARS);

        treeModel = new TreeModel(tree); //treeModel

        String expectedNewickTree = "((((human:0.024003,(chimp:0.010772,bonobo:0.010772):0.013231):0.012035," +
                "gorilla:0.036038):0.033087,orangutan:0.069125):0.030457,siamang:0.099582);";
        assertEquals("Fail to covert the correct tree !!!", expectedNewickTree, Tree.Utils.newick(treeModel, 6));
    }



    public void testLikelihood() {
/*        // Sub model
        Parameter freqs = new Parameter.Default(alignment.getStateFrequencies());//new double[]{0.25, 0.25, 0.25, 0.25});
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, f);

        //siteModel
        GammaSiteModel siteModel = new GammaSiteModel(hky);
        Parameter mu = new Parameter.Default(GammaSiteModel.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        siteModel.setMutationRateParameter(mu);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, null,
                false, false, true, false, false);
        treeLikelihood.setId("treeLikelihood");

        // Operators
        OperatorSchedule schedule = new SimpleOperatorSchedule();

        MCMCOperator operator = new ScaleOperator(kappa, 0.5);
        operator.setWeight(1.0);
        schedule.addOperator(operator);

//        Parameter rootParameter = treeModel.createNodeHeightsParameter(true, false, false);
//        ScaleOperator scaleOperator = new ScaleOperator(rootParameter, 0.75, CoercionMode.COERCION_ON, 1.0);

        Parameter rootHeight = treeModel.getRootHeightParameter();
        rootHeight.setId(TREE_HEIGHT);
        operator = new ScaleOperator(rootHeight, 0.5);
        operator.setWeight(1.0);
        schedule.addOperator(operator);

        Parameter internalHeights = treeModel.createNodeHeightsParameter(false, true, false);
        operator = new UniformOperator(internalHeights, 10.0);
        schedule.addOperator(operator);

        operator = new SubtreeSlideOperator(treeModel, 1, 1, true, false, false, false, CoercionMode.COERCION_ON);
        schedule.addOperator(operator);

        operator = new ExchangeOperator(ExchangeOperator.NARROW, treeModel, 1.0);
//        operator.doOperation();
        schedule.addOperator(operator);

        operator = new ExchangeOperator(ExchangeOperator.WIDE, treeModel, 1.0);
//        operator.doOperation();
        schedule.addOperator(operator);

        operator = new WilsonBalding(treeModel, 1.0);
//        operator.doOperation();
        schedule.addOperator(operator);

        // Log
        ArrayLogFormatter formatter = new ArrayLogFormatter(false);

        MCLogger[] loggers = new MCLogger[2];
        loggers[0] = new MCLogger(formatter, 1000, false);
        loggers[0].add(treeLikelihood);
        loggers[0].add(rootHeight);
        loggers[0].add(kappa);

        loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), 100000, false);
        loggers[1].add(treeLikelihood);
        loggers[1].add(rootHeight);
        loggers[1].add(kappa);

        // MCMC
        MCMC mcmc = new MCMC("mcmc1");
        MCMCOptions options = new MCMCOptions();
        options.setChainLength(10000000);
        options.setUseCoercion(true); // autoOptimize = true
        options.setCoercionDelay(100);
        options.setTemperature(1.0);
        options.setFullEvaluationCount(2000);

        mcmc.setShowOperatorAnalysis(true);
        mcmc.init(options, treeLikelihood, Prior.UNIFORM_PRIOR, schedule, loggers);
        mcmc.run();

        // time
        System.out.println(mcmc.getTimer().toString());

        // Tracer
        List<Trace> traces = formatter.getTraces();
        ArrayTraceList traceList = new ArrayTraceList("LikelihoodTest", traces, 0);

        for (int i = 1; i < traces.size(); i++) {
            traceList.analyseTrace(i);
        }

//      <expectation name="likelihood" value="-1815.75"/>
//		<expectation name="treeModel.rootHeight" value="6.42048E-2"/>
//		<expectation name="hky.kappa" value="32.8941"/>

        TraceCorrelation likelihoodStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(TreeLikelihood.TREE_LIKELIHOOD));
        assertExpectation(TreeLikelihood.TREE_LIKELIHOOD, likelihoodStats, -1815.75);

        TraceCorrelation treeHeightStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(TREE_HEIGHT));
        assertExpectation(TREE_HEIGHT, treeHeightStats, 6.42048E-2);

        TraceCorrelation kappaStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(HKYParser.KAPPA));
        assertExpectation(HKYParser.KAPPA, kappaStats, 32.8941);  */
    }

    public static Test suite() {
        return new TestSuite(LikelihoodTest.class);
    }
}
