package test.dr.evomodel.treelikelihood;

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
import dr.evomodel.treelikelihood.SequenceErrorModel;
import dr.evomodel.treelikelihood.TipPartialsModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.evomodelxml.HKYParser;
import dr.inference.loggers.ArrayLogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.OneOnXPrior;
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
 * convert testPMD.xml in the folder /example
 */
public class PMDTestProblem extends TraceCorrelationAssert {

    public PMDTestProblem(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        MathUtils.setSeed(666);

        createAlignment(NUMBER_TAXON_SEQUENCE, Nucleotides.INSTANCE);
    }


    public void testPMD() throws Exception {
        Parameter popSize = new Parameter.Default(ConstantPopulationModel.POPULATION_SIZE, 496432.69917113904, 0, Double.POSITIVE_INFINITY);
        ConstantPopulationModel constantModel = createRandomInitialTree(popSize);

        CoalescentLikelihood coalescent = new CoalescentLikelihood(treeModel, null, new ArrayList<TaxonList>(), constantModel);
        coalescent.setId("coalescent");

        // clock model
        Parameter rateParameter =  new Parameter.Default(StrictClockBranchRates.RATE, 4.0E-7, 0, 100.0);
        StrictClockBranchRates branchRateModel = new StrictClockBranchRates(rateParameter);

        // Sub model
        Parameter freqs = new Parameter.Default(new double[]{0.25, 0.25, 0.25, 0.25});
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, f);

        //siteModel
        GammaSiteModel siteModel = new GammaSiteModel(hky);
        Parameter mu = new Parameter.Default(GammaSiteModel.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        siteModel.setMutationRateParameter(mu);

        // SequenceErrorModel
        Parameter ageRelatedRateParameter = new Parameter.Default(SequenceErrorModel.AGE_RELATED_RATE, 4.0E-7, 0, 100.0);
        
        TipPartialsModel aDNADamageModel =  new SequenceErrorModel(null, null, SequenceErrorModel.ErrorType.TRANSITIONS_ONLY,
                null, ageRelatedRateParameter);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, branchRateModel, aDNADamageModel,
                false, false, true, false, false);
        treeLikelihood.setId(TreeLikelihood.TREE_LIKELIHOOD);

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

        operator = new ScaleOperator(ageRelatedRateParameter, 0.75);
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

        operator = new SubtreeSlideOperator(treeModel, 15.0, 49643.2699171139, true, false, false, false, CoercionMode.COERCION_ON);
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

        operator = new DeltaExchangeOperator(freqs, new int[] {1, 1, 1, 1}, 0.01, 1.0, false, CoercionMode.COERCION_ON); // ??? correct?
        schedule.addOperator(operator);

        //CompoundLikelihood
        OneOnXPrior likelihood1 = new OneOnXPrior();
        likelihood1.addData(popSize);
        OneOnXPrior likelihood2 = new OneOnXPrior();
        likelihood2.addData(kappa);

        List<Likelihood> likelihoods = new ArrayList<Likelihood>();
        likelihoods.add(likelihood1);
        likelihoods.add(likelihood2);
        likelihoods.add(coalescent);
        Likelihood prior = new CompoundLikelihood(0, likelihoods);
        prior.setId(CompoundLikelihood.PRIOR);

        likelihoods.clear();
        likelihoods.add(treeLikelihood);
        Likelihood likelihood = new CompoundLikelihood(-1, likelihoods);

        likelihoods.clear();
        likelihoods.add(prior);
        likelihoods.add(likelihood);
        Likelihood posterior = new CompoundLikelihood(0, likelihoods);
        posterior.setId(CompoundLikelihood.POSTERIOR);

        // Log
        ArrayLogFormatter formatter = new ArrayLogFormatter(false);

        MCLogger[] loggers = new MCLogger[2];
        loggers[0] = new MCLogger(formatter, 1000, false);
        loggers[0].add(posterior);
        loggers[0].add(treeLikelihood);
        loggers[0].add(rootHeight);
        loggers[0].add(rateParameter);
        loggers[0].add(ageRelatedRateParameter);
        loggers[0].add(popSize);
        loggers[0].add(kappa);
        loggers[0].add(coalescent);

        loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), 10000, false);
        loggers[1].add(posterior);
        loggers[1].add(treeLikelihood);
        loggers[1].add(rootHeight);
        loggers[1].add(rateParameter);

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
        ArrayTraceList traceList = new ArrayTraceList("PMDTest", traces, 0);

        for (int i = 1; i < traces.size(); i++) {
            traceList.analyseTrace(i);
        }

//        <expectation name="clock.rate" value="1.5E-7"/>
//        <expectation name="errorModel.ageRate" value="0.7E-7"/>
//        <expectation name="hky.kappa" value="10"/>

        TraceCorrelation kappaStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(HKYParser.KAPPA));
        assertExpectation(HKYParser.KAPPA, kappaStats, 10);

        TraceCorrelation rateStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(StrictClockBranchRates.RATE));
        assertExpectation(StrictClockBranchRates.RATE, rateStats, 1.5E-7);

        TraceCorrelation ageRateStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(SequenceErrorModel.AGE_RELATED_RATE));
        assertExpectation(SequenceErrorModel.AGE_RELATED_RATE, ageRateStats, 0.7E-7);
    }

    public static Test suite() {
        return new TestSuite(PMDTestProblem.class);
    }
}
