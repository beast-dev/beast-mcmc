package test.dr.integration;

import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.GeneralDataType;
import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.operators.SubtreeSlideOperator;
import dr.evomodel.operators.WilsonBalding;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.GeneralSubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.evomodelxml.HKYParser;
import dr.evomodelxml.TreeModelParser;
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
import junit.framework.Test;
import junit.framework.TestSuite;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Walter Xie
 * convert testGeneralSubstitutionModel.xml in the folder /example
 */
public class GeneralSubsitutionModelTest extends TraceCorrelationAssert {

    private static final String TREE_HEIGHT = TreeModel.TREE_MODEL + "." + TreeModelParser.ROOT_HEIGHT;

    private TreeModel treeModel;
    private SimpleAlignment alignment;
    private GeneralDataType dataType;

    public GeneralSubsitutionModelTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        alignment = createAlignment();

        List<String> states = new ArrayList<String>();
        states.addAll(Arrays.asList("A", "C", "G", "T"));
        dataType = new GeneralDataType(states);
        dataType.addAmbiguity("-", new String[] {"ACGT"});

        alignment.setDataType(dataType);

        treeModel = createTree("(((((chimp:0.010464222027296717,bonobo:0.010464222027296717):0.010716369046616688," +
                "human:0.021180591073913405):0.010988083344422011,gorilla:0.032168674418335416):0.022421978632286572," +
                "orangutan:0.05459065305062199):0.009576302472349953,siamang:0.06416695552297194);");
        
    }


    public void testMCMC() {

        // Sub model
        FrequencyModel freqModel = new FrequencyModel(dataType, alignment.getStateFrequencies());
        Parameter ratesPara = new Parameter.Default(GeneralSubstitutionModel.RATES, 5, 1.0); // dimension="5" value="1.0"
        GeneralSubstitutionModel generalSubstitutionModel = new GeneralSubstitutionModel(dataType, freqModel, ratesPara, 5); // relativeTo="5"

        //siteModel
        GammaSiteModel siteModel = new GammaSiteModel(generalSubstitutionModel);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);
        
        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, null,
                false, false, true, false, false);
        treeLikelihood.setId("treeLikelihood");

        // Operators
        OperatorSchedule schedule = new SimpleOperatorSchedule();

        MCMCOperator operator = new ScaleOperator(ratesPara, 0.5);
        operator.setWeight(1.0);
        schedule.addOperator(operator);

//        operator = new ScaleOperator(rootHeight, 0.5);
//        operator.setWeight(1.0);
//        schedule.addOperator(operator);

//        Parameter rootParameter = treeModel.createNodeHeightsParameter(true, false, false);
//        ScaleOperator scaleOperator = new ScaleOperator(rootParameter, 0.75, CoercionMode.COERCION_ON, 1.0);

        Parameter rootHeight = treeModel.getRootHeightParameter();
        operator = new ScaleOperator(rootHeight, 0.5);
        operator.setWeight(1.0);
        schedule.addOperator(operator);

        Parameter internalHeights = treeModel.createNodeHeightsParameter(false, true, false);
        operator = new UniformOperator(internalHeights, 1.0);
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
        loggers[0].add(ratesPara);

        loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), 100000, false);
        loggers[1].add(treeLikelihood);
        loggers[1].add(rootHeight);
        loggers[1].add(ratesPara);

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
        mcmc.getTimer();

        // Tracer
        List<Trace> traces = formatter.getTraces();
        ArrayTraceList traceList = new ArrayTraceList("MCMCTest", traces, 0);

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
        assertExpectation(HKYParser.KAPPA, kappaStats, 32.8941);
    }

    public static Test suite() {
        return new TestSuite(GeneralSubsitutionModelTest.class);
    }
}