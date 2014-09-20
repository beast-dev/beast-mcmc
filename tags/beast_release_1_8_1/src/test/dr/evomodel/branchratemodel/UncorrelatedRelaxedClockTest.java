package test.dr.evomodel.branchratemodel;

import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.DiscretizedBranchRates;
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
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.evomodelxml.sitemodel.GammaSiteModelParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.evomodelxml.tree.RateStatisticParser;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.inference.distribution.ExponentialDistributionModel;
import dr.inference.distribution.LogNormalDistributionModel;
import dr.inference.distribution.ParametricDistributionModel;
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
import dr.inferencexml.distribution.DistributionModelParser;
import dr.inferencexml.distribution.LogNormalDistributionModelParser;
import dr.inferencexml.model.CompoundLikelihoodParser;
import dr.math.MathUtils;
import junit.framework.Test;
import junit.framework.TestSuite;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Walter Xie
 * convert testUncorrelatedRelaxedClock.xml in the folder /example
 */
public class UncorrelatedRelaxedClockTest extends TraceCorrelationAssert {

    private Parameter meanParam;
    private Parameter stdevParam;

    public UncorrelatedRelaxedClockTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        MathUtils.setSeed(666);

        createAlignment(DENGUE4_TAXON_SEQUENCE, Nucleotides.INSTANCE);
    }

    public void testLogNormal() throws Exception {
        meanParam = new Parameter.Default(LogNormalDistributionModelParser.MEAN, 2.3E-5, 0, 100.0);
        stdevParam = new Parameter.Default(LogNormalDistributionModelParser.STDEV, 0.1, 0, 10.0);
        ParametricDistributionModel distributionModel = new LogNormalDistributionModel(meanParam, stdevParam, 0.0, true, false); // meanInRealSpace="true"

        ArrayTraceList traceList = UncorrelatedRelaxedClock(distributionModel);

//        <expectation name="posterior" value="-3927.81"/>
//        <expectation name="ucld.mean" value="8.28472E-4"/>
//        <expectation name="ucld.stdev" value="0.17435"/>
//        <expectation name="meanRate" value="8.09909E-4"/>
//        <expectation name="coefficientOfVariation" value="0.15982"/>
//        <expectation name="covariance" value="-3.81803E-2"/>
//        <expectation name="constant.popSize" value="37.3524"/>
//        <expectation name="hky.kappa" value="18.3053"/>
//        <expectation name="treeModel.rootHeight" value="69.2953"/>
//        <expectation name="treeLikelihood" value="-3855.78"/>
//        <expectation name="skyline" value="-72.0313"/>   ???

        TraceCorrelation likelihoodStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(CompoundLikelihoodParser.POSTERIOR));
        assertExpectation(CompoundLikelihoodParser.POSTERIOR, likelihoodStats, -3927.81);

        likelihoodStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(TreeLikelihoodParser.TREE_LIKELIHOOD));
        assertExpectation(TreeLikelihoodParser.TREE_LIKELIHOOD, likelihoodStats, -3855.78);

        TraceCorrelation treeHeightStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(TREE_HEIGHT));
        assertExpectation(TREE_HEIGHT, treeHeightStats, 69.2953);

        TraceCorrelation kappaStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(HKYParser.KAPPA));
        assertExpectation(HKYParser.KAPPA, kappaStats, 18.06518);

        TraceCorrelation ucldStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(LogNormalDistributionModelParser.MEAN));
        assertExpectation(LogNormalDistributionModelParser.MEAN, ucldStats, 8.0591451486E-4);

        ucldStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(LogNormalDistributionModelParser.STDEV));
        assertExpectation(LogNormalDistributionModelParser.STDEV, ucldStats, 0.16846023066431434);

        TraceCorrelation rateStats = traceList.getCorrelationStatistics(traceList.getTraceIndex("meanRate"));
        assertExpectation("meanRate", rateStats, 8.010906E-4);

        TraceCorrelation coefficientOfVariationStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(RateStatisticParser.COEFFICIENT_OF_VARIATION));
        assertExpectation(RateStatisticParser.COEFFICIENT_OF_VARIATION, coefficientOfVariationStats, 0.15982);

        TraceCorrelation covarianceStats = traceList.getCorrelationStatistics(traceList.getTraceIndex("covariance"));
        assertExpectation("covariance", covarianceStats, -0.0260333026);

        TraceCorrelation popStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(ConstantPopulationModelParser.POPULATION_SIZE));
        assertExpectation(ConstantPopulationModelParser.POPULATION_SIZE, popStats, 37.3524);

        TraceCorrelation coalescentStats = traceList.getCorrelationStatistics(traceList.getTraceIndex("coalescent"));
        assertExpectation("coalescent", coalescentStats, -72.0313);
    }

    public void testExponential() throws Exception {
        meanParam = new Parameter.Default(1.0);
        meanParam.setId(DistributionModelParser.MEAN); 
        stdevParam = null;
        ParametricDistributionModel distributionModel = new ExponentialDistributionModel(meanParam); // offset = 0

        ArrayTraceList traceList = UncorrelatedRelaxedClock(distributionModel);

        TraceCorrelation likelihoodStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(CompoundLikelihoodParser.POSTERIOR));
        assertExpectation(CompoundLikelihoodParser.POSTERIOR, likelihoodStats, -3958.7409);
//        System.out.println("likelihoodStats = " + likelihoodStats.getMean());

        likelihoodStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(TreeLikelihoodParser.TREE_LIKELIHOOD));
        assertExpectation(TreeLikelihoodParser.TREE_LIKELIHOOD, likelihoodStats, -3885.26939);
//        System.out.println("treelikelihoodStats = " + likelihoodStats.getMean());

        TraceCorrelation treeHeightStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(TREE_HEIGHT));
        assertExpectation(TREE_HEIGHT, treeHeightStats, 84.3529526);

        TraceCorrelation kappaStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(HKYParser.KAPPA));
        assertExpectation(HKYParser.KAPPA, kappaStats, 18.38065);

        TraceCorrelation ucedStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(DistributionModelParser.MEAN));
        assertExpectation(DistributionModelParser.MEAN, ucedStats, 0.0019344134887784579);
//        System.out.println("ucedStats = " + ucedStats.getMean());

        TraceCorrelation rateStats = traceList.getCorrelationStatistics(traceList.getTraceIndex("meanRate"));
        assertExpectation("meanRate", rateStats, 0.0020538802366337084);
//        System.out.println("rateStats = " + rateStats.getMean());

        TraceCorrelation coefficientOfVariationStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(RateStatisticParser.COEFFICIENT_OF_VARIATION));
        assertExpectation(RateStatisticParser.COEFFICIENT_OF_VARIATION, coefficientOfVariationStats, 0.773609960455);
//        System.out.println("coefficientOfVariationStats = " + coefficientOfVariationStats.getMean());

        TraceCorrelation covarianceStats = traceList.getCorrelationStatistics(traceList.getTraceIndex("covariance"));
        assertExpectation("covariance", covarianceStats, -0.07042030641301375);
//        System.out.println("covarianceStats = " + covarianceStats.getMean());

        TraceCorrelation popStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(ConstantPopulationModelParser.POPULATION_SIZE));
        assertExpectation(ConstantPopulationModelParser.POPULATION_SIZE, popStats, 43.4478);
//        System.out.println("popStats = " + popStats.getMean());

        TraceCorrelation coalescentStats = traceList.getCorrelationStatistics(traceList.getTraceIndex("coalescent"));
        assertExpectation("coalescent", coalescentStats, -73.4715);
//        System.out.println("coalescentStats = " + coalescentStats.getMean());
    }

    private ArrayTraceList UncorrelatedRelaxedClock(ParametricDistributionModel distributionModel) throws Exception {
        Parameter popSize = new Parameter.Default(ConstantPopulationModelParser.POPULATION_SIZE, 380.0, 0, 38000.0);
        ConstantPopulationModel constantModel = createRandomInitialTree(popSize);

        CoalescentLikelihood coalescent = new CoalescentLikelihood(treeModel, null, new ArrayList<TaxonList>(), constantModel);
        coalescent.setId("coalescent");

        // clock model        
        Parameter rateCategoryParameter = new Parameter.Default(32);
        rateCategoryParameter.setId(DiscretizedBranchRates.BRANCH_RATES); 

        DiscretizedBranchRates branchRateModel = new DiscretizedBranchRates(treeModel, rateCategoryParameter, 
                distributionModel, 1, false, Double.NaN, false, false);

        RateStatistic meanRate = new RateStatistic("meanRate", treeModel, branchRateModel, true, true, RateStatisticParser.MEAN);
        RateStatistic coefficientOfVariation = new RateStatistic(RateStatisticParser.COEFFICIENT_OF_VARIATION, treeModel, branchRateModel,
                true, true, RateStatisticParser.COEFFICIENT_OF_VARIATION);
        RateCovarianceStatistic covariance = new RateCovarianceStatistic("covariance", treeModel, branchRateModel);

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

        operator = new ScaleOperator(meanParam, 0.75);
        operator.setWeight(3.0);
        schedule.addOperator(operator);

        if (stdevParam != null) {
            operator = new ScaleOperator(stdevParam, 0.75);
            operator.setWeight(3.0);
            schedule.addOperator(operator);
        }

        Parameter allInternalHeights = treeModel.createNodeHeightsParameter(true, true, false);
        operator = new UpDownOperator(new Scalable[]{new Scalable.Default(meanParam)},
                new Scalable[] {new Scalable.Default(allInternalHeights)}, 0.75, 3.0, CoercionMode.COERCION_ON);
        schedule.addOperator(operator);

        operator = new SwapOperator(rateCategoryParameter, 10);
        operator.setWeight(1.0);
        schedule.addOperator(operator);

        operator = new RandomWalkIntegerOperator(rateCategoryParameter, 1, 10.0);
        schedule.addOperator(operator);

        operator = new UniformIntegerOperator(rateCategoryParameter, (int) (double)rateCategoryParameter.getBounds().getLowerLimit(0),
                (int) (double)rateCategoryParameter.getBounds().getUpperLimit(0), 10.0);
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

        operator = new SubtreeSlideOperator(treeModel, 15.0, 38.0, true, false, false, false, CoercionMode.COERCION_ON);
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
        loggers[0] = new MCLogger(formatter, 10000, false);
        loggers[0].add(posterior);
        loggers[0].add(treeLikelihood);
        loggers[0].add(rootHeight);
        loggers[0].add(meanParam);
        if (stdevParam != null) loggers[0].add(stdevParam);
        loggers[0].add(meanRate);
        loggers[0].add(coefficientOfVariation);
        loggers[0].add(covariance);
        loggers[0].add(popSize);
        loggers[0].add(kappa);
        loggers[0].add(coalescent);

        loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), 100000, false);
        loggers[1].add(posterior);
        loggers[1].add(treeLikelihood);
        loggers[1].add(rootHeight);
        loggers[1].add(meanRate);
        loggers[1].add(coalescent);

        // MCMC
        MCMC mcmc = new MCMC("mcmc1");
        MCMCOptions options = new MCMCOptions(10000000);

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

        return traceList;
    }

    public static Test suite() {
        return new TestSuite(UncorrelatedRelaxedClockTest.class);
    }
}

