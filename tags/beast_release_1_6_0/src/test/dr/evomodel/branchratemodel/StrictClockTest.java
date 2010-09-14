package test.dr.evomodel.branchratemodel;

import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.coalescent.CoalescentLikelihood;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.operators.SubtreeSlideOperator;
import dr.evomodel.operators.WilsonBalding;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HKY;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.evomodelxml.sitemodel.GammaSiteModelParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.inference.loggers.ArrayLogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.TabDelimitedFormatter;
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
import dr.math.MathUtils;
import junit.framework.Test;
import junit.framework.TestSuite;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Walter Xie
 * convert testStrictClock.xml in the folder /example
 */
public class StrictClockTest extends TraceCorrelationAssert {

    public StrictClockTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        MathUtils.setSeed(666);

        createAlignment(DENGUE4_TAXON_SEQUENCE, Nucleotides.INSTANCE);
    }


    public void testStrictClock() throws Exception {
        Parameter popSize = new Parameter.Default(ConstantPopulationModelParser.POPULATION_SIZE, 380.0, 0, 38000.0);
        ConstantPopulationModel constantModel = createRandomInitialTree(popSize);

        CoalescentLikelihood coalescent = new CoalescentLikelihood(treeModel, null, new ArrayList<TaxonList>(), constantModel);
        coalescent.setId("coalescent");

        // clock model
        Parameter rateParameter =  new Parameter.Default(StrictClockBranchRates.RATE, 2.3E-5, 0, 100.0);

        StrictClockBranchRates branchRateModel = new StrictClockBranchRates(rateParameter);

        // Sub model
        Parameter freqs = new Parameter.Default(alignment.getStateFrequencies());
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 0, 100.0);

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, f);

        //siteModel
        GammaSiteModel siteModel = new GammaSiteModel(hky);
        Parameter mu = new Parameter.Default(GammaSiteModelParser.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        siteModel.setMutationRateParameter(mu);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, branchRateModel, null,
                false, false, true, false, false);
        treeLikelihood.setId(TreeLikelihoodParser.TREE_LIKELIHOOD);

        // Operators
        OperatorSchedule schedule = new SimpleOperatorSchedule();

        MCMCOperator operator = new ScaleOperator(kappa, 0.75);
        operator.setWeight(1.0);
        schedule.addOperator(operator);

        operator = new ScaleOperator(rateParameter, 0.75);
        operator.setWeight(3.0);
        schedule.addOperator(operator);

        Parameter allInternalHeights = treeModel.createNodeHeightsParameter(true, true, false);
        operator = new UpDownOperator(new Scalable[]{new Scalable.Default(rateParameter)},
                new Scalable[] {new Scalable.Default(allInternalHeights)}, 0.75, 3.0, CoercionMode.COERCION_ON);
        schedule.addOperator(operator);

        operator = new ScaleOperator(popSize, 0.75);
        operator.setWeight(3.0);
        schedule.addOperator(operator);

        Parameter rootHeight = treeModel.getRootHeightParameter();
        rootHeight.setId(TREE_HEIGHT);
        operator = new ScaleOperator(rootHeight, 0.75);
        operator.setWeight(3.0);
        schedule.addOperator(operator);

        Parameter internalHeights = treeModel.createNodeHeightsParameter(false, true, false);
        operator = new UniformOperator(internalHeights, 30.0);
        schedule.addOperator(operator);

        operator = new SubtreeSlideOperator(treeModel, 15.0, 1.0, true, false, false, false, CoercionMode.COERCION_ON);
        schedule.addOperator(operator);

        operator = new ExchangeOperator(ExchangeOperator.NARROW, treeModel, 15.0);
//        operator.doOperation();
        schedule.addOperator(operator);

        operator = new ExchangeOperator(ExchangeOperator.WIDE, treeModel, 3.0);
//        operator.doOperation();
        schedule.addOperator(operator);

        operator = new WilsonBalding(treeModel, 3.0);
//        operator.doOperation();
        schedule.addOperator(operator);

        //CompoundLikelihood
        List<Likelihood> likelihoods = new ArrayList<Likelihood>();        
        likelihoods.add(coalescent);
        Likelihood prior = new CompoundLikelihood(0, likelihoods);
        prior.setId(CompoundLikelihoodParser.PRIOR);

        likelihoods.clear();
        likelihoods.add(treeLikelihood);
        Likelihood likelihood = new CompoundLikelihood(-1, likelihoods);

        likelihoods.clear();
        likelihoods.add(prior);
        likelihoods.add(likelihood);
        Likelihood posterior = new CompoundLikelihood(0, likelihoods);
        posterior.setId(CompoundLikelihoodParser.POSTERIOR);

        // Log
        ArrayLogFormatter formatter = new ArrayLogFormatter(false);

        MCLogger[] loggers = new MCLogger[2];
        loggers[0] = new MCLogger(formatter, 500, false);
        loggers[0].add(posterior);
        loggers[0].add(treeLikelihood);
        loggers[0].add(rootHeight);
        loggers[0].add(rateParameter);
        loggers[0].add(popSize);
        loggers[0].add(kappa);
        loggers[0].add(coalescent);

        loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), 10000, false);
        loggers[1].add(posterior);
        loggers[1].add(treeLikelihood);
        loggers[1].add(rootHeight);
        loggers[1].add(rateParameter);
        loggers[1].add(coalescent);

        // MCMC
        MCMC mcmc = new MCMC("mcmc1");
        MCMCOptions options = new MCMCOptions();
        options.setChainLength(1000000);
        options.setUseCoercion(true); // autoOptimize = true
        options.setCoercionDelay(100);
        options.setTemperature(1.0);
        options.setFullEvaluationCount(2000);

        mcmc.setShowOperatorAnalysis(true);
        mcmc.init(options, posterior, schedule, loggers);
        mcmc.run();

        // time
        System.out.println(mcmc.getTimer().toString());

        // Tracer
        List<Trace> traces = formatter.getTraces();
        ArrayTraceList traceList = new ArrayTraceList("RandomLocalClockTest", traces, 0);

        for (int i = 1; i < traces.size(); i++) {
            traceList.analyseTrace(i);
        }

//        <expectation name="posterior" value="-3928.71"/>
//        <expectation name="clock.rate" value="8.04835E-4"/>
//        <expectation name="constant.popSize" value="37.3762"/>
//        <expectation name="hky.kappa" value="18.2782"/>
//        <expectation name="treeModel.rootHeight" value="69.0580"/>
//        <expectation name="treeLikelihood" value="-3856.59"/>
//        <expectation name="coalescent" value="-72.1285"/>

        TraceCorrelation likelihoodStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(CompoundLikelihoodParser.POSTERIOR));
        assertExpectation(CompoundLikelihoodParser.POSTERIOR, likelihoodStats, -3928.71);

        likelihoodStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(TreeLikelihoodParser.TREE_LIKELIHOOD));
        assertExpectation(TreeLikelihoodParser.TREE_LIKELIHOOD, likelihoodStats, -3856.59);

        TraceCorrelation treeHeightStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(TREE_HEIGHT));
        assertExpectation(TREE_HEIGHT, treeHeightStats, 69.0580);

        TraceCorrelation kappaStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(HKYParser.KAPPA));
        assertExpectation(HKYParser.KAPPA, kappaStats, 18.2782);

        TraceCorrelation rateStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(StrictClockBranchRates.RATE));
        assertExpectation(StrictClockBranchRates.RATE, rateStats, 8.04835E-4);        

        TraceCorrelation popStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(ConstantPopulationModelParser.POPULATION_SIZE));
        assertExpectation(ConstantPopulationModelParser.POPULATION_SIZE, popStats, 37.3762);

        TraceCorrelation coalescentStats = traceList.getCorrelationStatistics(traceList.getTraceIndex("coalescent"));
        assertExpectation("coalescent", coalescentStats, -72.1285);
    }

    public static Test suite() {
        return new TestSuite(StrictClockTest.class);
    }
}

