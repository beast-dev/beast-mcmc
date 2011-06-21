package test.dr.evomodel.substmodel;

import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.GeneralDataType;
import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.operators.SubtreeSlideOperator;
import dr.evomodel.operators.WilsonBalding;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.GeneralSubstitutionModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.evomodelxml.sitemodel.GammaSiteModelParser;
import dr.evomodelxml.substmodel.GeneralSubstitutionModelParser;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.inference.loggers.ArrayLogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.inference.trace.ArrayTraceList;
import dr.inference.trace.Trace;
import dr.inference.trace.TraceCorrelation;
import dr.math.MathUtils;
import junit.framework.Test;
import junit.framework.TestSuite;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Walter Xie
 * convert testGeneralSubstitutionModel.xml in the folder /example
 */
public class GeneralSubstitutionModelTest extends TraceCorrelationAssert {

    private GeneralDataType dataType;

    public GeneralSubstitutionModelTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        MathUtils.setSeed(666);

        List<String> states = new ArrayList<String>();
        states.add("A");
        states.add("C");
        states.add("G");
        states.add("T");
        //states.addAll(Arrays.asList("A", "C", "G", "T"));
        dataType = new GeneralDataType(states);
        dataType.addAmbiguity("-", new String[] {"A", "C", "G", "T"});
        dataType.addAmbiguity("?", new String[] {"A", "C", "G", "T"});

        createAlignment(PRIMATES_TAXON_SEQUENCE, dataType);


        createRandomInitialTree(0.0001); // popSize

//        createSpecifiedTree("(((((chimp:0.010464222027296717,bonobo:0.010464222027296717):0.010716369046616688," +
//                "human:0.021180591073913405):0.010988083344422011,gorilla:0.032168674418335416):0.022421978632286572," +
//                "orangutan:0.05459065305062199):0.009576302472349953,siamang:0.06416695552297194);");
        
    }


    public void testGeneralSubstitutionModel() {

        // Sub model
        FrequencyModel freqModel = new FrequencyModel(dataType, alignment.getStateFrequencies());
        Parameter ratesPara = new Parameter.Default(GeneralSubstitutionModelParser.RATES, 5, 1.0); // dimension="5" value="1.0"
        GeneralSubstitutionModel generalSubstitutionModel = new GeneralSubstitutionModel(dataType, freqModel, ratesPara, 4); // relativeTo="5"

        //siteModel
        GammaSiteModel siteModel = new GammaSiteModel(generalSubstitutionModel);
        Parameter mu = new Parameter.Default(GammaSiteModelParser.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        siteModel.setMutationRateParameter(mu);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);
        
        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, null,
                false, false, true, false, false);
        treeLikelihood.setId(TreeLikelihoodParser.TREE_LIKELIHOOD);

        // Operators
        OperatorSchedule schedule = new SimpleOperatorSchedule();

        MCMCOperator operator = new ScaleOperator(ratesPara, 0.5);
        operator.setWeight(1.0);
        schedule.addOperator(operator);

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
        mcmc.init(options, treeLikelihood, schedule, loggers);
        mcmc.run();
        
        // time
        System.out.println(mcmc.getTimer().toString());

        // Tracer
        List<Trace> traces = formatter.getTraces();
        ArrayTraceList traceList = new ArrayTraceList("GeneralSubstitutionModelTest", traces, 0);

        for (int i = 1; i < traces.size(); i++) {
            traceList.analyseTrace(i);
        }

//      <expectation name="likelihood" value="-1815.75"/>
//		<expectation name="treeModel.rootHeight" value="6.42048E-2"/>
//		<expectation name="rateAC" value="6.08986E-2"/>

        TraceCorrelation likelihoodStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(TreeLikelihoodParser.TREE_LIKELIHOOD));
        assertExpectation(TreeLikelihoodParser.TREE_LIKELIHOOD, likelihoodStats, -1815.75);

        TraceCorrelation treeHeightStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(TREE_HEIGHT));
        assertExpectation(TREE_HEIGHT, treeHeightStats, 0.0640787258170083);

        TraceCorrelation rateACStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(GeneralSubstitutionModelParser.RATES + "1"));
        assertExpectation(GeneralSubstitutionModelParser.RATES + "1", rateACStats, 0.061071756742081366);
    }

    public static Test suite() {
        return new TestSuite(GeneralSubstitutionModelTest.class);
    }
}