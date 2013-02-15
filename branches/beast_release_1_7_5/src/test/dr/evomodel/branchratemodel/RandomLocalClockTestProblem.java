package test.dr.evomodel.branchratemodel;

import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.RandomLocalClockModel;
import dr.evomodel.coalescent.CoalescentLikelihood;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.operators.SubtreeSlideOperator;
import dr.evomodel.operators.WilsonBalding;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HKY;
import dr.evomodel.tree.RateCovarianceStatistic;
import dr.evomodel.tree.RateStatistic;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.evomodelxml.branchratemodel.RandomLocalClockModelParser;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.evomodelxml.sitemodel.GammaSiteModelParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.evomodelxml.tree.RateStatisticParser;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.loggers.ArrayLogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.*;
import dr.inference.operators.*;
import dr.inference.trace.ArrayTraceList;
import dr.inference.trace.Trace;
import dr.inference.trace.TraceCorrelation;
import dr.inferencexml.model.CompoundLikelihoodParser;
import dr.math.MathUtils;
import dr.math.distributions.GammaDistribution;
import dr.math.distributions.PoissonDistribution;
import junit.framework.Test;
import junit.framework.TestSuite;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Walter Xie
 * convert testRandomLocalClock.xml in the folder /example
 * the Random Local Clock model has problem
 */

public class RandomLocalClockTestProblem extends TraceCorrelationAssert {

    public RandomLocalClockTestProblem(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        MathUtils.setSeed(666);

        createAlignment(PRIMATES_TAXON_SEQUENCE, Nucleotides.INSTANCE);
    }


    public void testRandomLocalClock() throws Exception {
        Parameter popSize = new Parameter.Default(ConstantPopulationModelParser.POPULATION_SIZE, 0.077, 0, Double.POSITIVE_INFINITY);
        ConstantPopulationModel constantModel = createRandomInitialTree(popSize);

        CoalescentLikelihood coalescent = new CoalescentLikelihood(treeModel, null, new ArrayList<TaxonList>(), constantModel);
        coalescent.setId("coalescent");

        // clock model
        Parameter ratesParameter = new Parameter.Default(RandomLocalClockModelParser.RATES, 10, 1);
        Parameter rateIndicatorParameter = new Parameter.Default(RandomLocalClockModelParser.RATE_INDICATORS, 10, 1);
        Parameter meanRateParameter = new Parameter.Default(RandomLocalClockModelParser.CLOCK_RATE, 1, 1.0);

        RandomLocalClockModel branchRateModel = new RandomLocalClockModel(treeModel, meanRateParameter,
                rateIndicatorParameter, ratesParameter, false);

        SumStatistic rateChanges = new SumStatistic("rateChangeCount", true);
        rateChanges.addStatistic(rateIndicatorParameter);

        RateStatistic meanRate = new RateStatistic("meanRate", treeModel, branchRateModel, true, true, RateStatisticParser.MEAN);
        RateStatistic coefficientOfVariation = new RateStatistic(RateStatisticParser.COEFFICIENT_OF_VARIATION, treeModel, branchRateModel,
                true, true, RateStatisticParser.COEFFICIENT_OF_VARIATION);
        RateCovarianceStatistic covariance = new RateCovarianceStatistic("covariance", treeModel, branchRateModel);

        // Sub model
        Parameter freqs = new Parameter.Default(alignment.getStateFrequencies());
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 0, Double.POSITIVE_INFINITY);

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

        operator = new ScaleOperator(ratesParameter, 0.75);
        operator.setWeight(10.0);
        schedule.addOperator(operator);

        operator = new BitFlipOperator(rateIndicatorParameter, 15.0, true);
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

        operator = new SubtreeSlideOperator(treeModel, 15.0, 0.0077, true, false, false, false, CoercionMode.COERCION_ON);
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
        OneOnXPrior likelihood1 = new OneOnXPrior();
        likelihood1.addData(popSize);
        OneOnXPrior likelihood2 = new OneOnXPrior();
        likelihood2.addData(kappa);
        DistributionLikelihood likelihood3 = new DistributionLikelihood(new GammaDistribution(0.5, 2.0), 0.0);
        likelihood3.addData(ratesParameter);
        DistributionLikelihood likelihood4 = new DistributionLikelihood(new PoissonDistribution(1.0), 0.0);
        likelihood4.addData(rateChanges);

        List<Likelihood> likelihoods = new ArrayList<Likelihood>();
        likelihoods.add(likelihood1);
        likelihoods.add(likelihood2);
        likelihoods.add(likelihood3);
        likelihoods.add(likelihood4);
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
        loggers[0] = new MCLogger(formatter, 1000, false);
        loggers[0].add(posterior);
        loggers[0].add(prior);
        loggers[0].add(treeLikelihood);
        loggers[0].add(rootHeight);
        loggers[0].add(kappa);
//        loggers[0].add(meanRate);
        loggers[0].add(rateChanges);
        loggers[0].add(coefficientOfVariation);
        loggers[0].add(covariance);
        loggers[0].add(popSize);
        loggers[0].add(coalescent);         

        loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), 10000, false);
        loggers[1].add(posterior);
        loggers[1].add(treeLikelihood);
        loggers[1].add(rootHeight);
        loggers[1].add(meanRate);
        loggers[1].add(rateChanges);

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

//        <expectation name="posterior" value="-1818.26"/>
//        <expectation name="prior" value="-2.70143"/>
//        <expectation name="likelihood" value="-1815.56"/>
//        <expectation name="treeModel.rootHeight" value="6.363E-2"/>
//        <expectation name="constant.popSize" value="9.67405E-2"/>
//        <expectation name="hky.kappa" value="30.0394"/>
//        <expectation name="coefficientOfVariation" value="7.02408E-2"/>
//        covariance 0.47952
//        <expectation name="rateChangeCount" value="0.40786"/>
//        <expectation name="coalescent" value="7.29521"/>

        TraceCorrelation likelihoodStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(CompoundLikelihoodParser.POSTERIOR));
        assertExpectation(CompoundLikelihoodParser.POSTERIOR, likelihoodStats, -1818.26);

        likelihoodStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(CompoundLikelihoodParser.PRIOR));
        assertExpectation(CompoundLikelihoodParser.PRIOR, likelihoodStats, -2.70143);

        likelihoodStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(TreeLikelihoodParser.TREE_LIKELIHOOD));
        assertExpectation(TreeLikelihoodParser.TREE_LIKELIHOOD, likelihoodStats, -1815.56);

        TraceCorrelation treeHeightStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(TREE_HEIGHT));
        assertExpectation(TREE_HEIGHT, treeHeightStats, 6.363E-2);

        TraceCorrelation kappaStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(HKYParser.KAPPA));
        assertExpectation(HKYParser.KAPPA, kappaStats, 30.0394);

        TraceCorrelation rateChangeStats = traceList.getCorrelationStatistics(traceList.getTraceIndex("rateChangeCount"));
        assertExpectation("rateChangeCount", rateChangeStats, 0.40786);

        TraceCorrelation coefficientOfVariationStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(RateStatisticParser.COEFFICIENT_OF_VARIATION));
        assertExpectation(RateStatisticParser.COEFFICIENT_OF_VARIATION, coefficientOfVariationStats, 7.02408E-2);

        TraceCorrelation covarianceStats = traceList.getCorrelationStatistics(traceList.getTraceIndex("covariance"));
        assertExpectation("covariance", covarianceStats, 0.47952);

        TraceCorrelation popStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(ConstantPopulationModelParser.POPULATION_SIZE));
        assertExpectation(ConstantPopulationModelParser.POPULATION_SIZE, popStats, 9.67405E-2);

        TraceCorrelation coalescentStats = traceList.getCorrelationStatistics(traceList.getTraceIndex("coalescent"));
        assertExpectation("coalescent", coalescentStats, 7.29521);
    }

    public static Test suite() {
        return new TestSuite(RandomLocalClockTestProblem.class);
    }
}
