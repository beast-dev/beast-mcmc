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
import dr.evomodel.treelikelihood.TipStatesModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.evomodelxml.sitemodel.GammaSiteModelParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.evomodelxml.treelikelihood.SequenceErrorModelParser;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;
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
import dr.inferencexml.model.CompoundLikelihoodParser;
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
        Parameter popSize = new Parameter.Default(ConstantPopulationModelParser.POPULATION_SIZE, 496432.69917113904, 0, Double.POSITIVE_INFINITY);
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
        Parameter mu = new Parameter.Default(GammaSiteModelParser.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        siteModel.setMutationRateParameter(mu);

        // SequenceErrorModel
        Parameter ageRelatedRateParameter = new Parameter.Default(SequenceErrorModelParser.AGE_RELATED_RATE, 4.0E-7, 0, 100.0);

        TipStatesModel aDNADamageModel =  new SequenceErrorModel(null, null, SequenceErrorModel.ErrorType.TRANSITIONS_ONLY,
                null, ageRelatedRateParameter, null);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, branchRateModel, aDNADamageModel,
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

        TraceCorrelation ageRateStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(SequenceErrorModelParser.AGE_RELATED_RATE));
        assertExpectation(SequenceErrorModelParser.AGE_RELATED_RATE, ageRateStats, 0.7E-7);
    }

    public static Test suite() {
        return new TestSuite(PMDTestProblem.class);
    }
}

/*
****************** BEAST result:

Operator analysis
Operator                                          Tuning  Count      Time     Time/Op  Pr(accept)  Performance suggestion
scale(clock.rate)                                 0.573   358210     685526   1.91     0.2564      good
up:clock.rate down:nodeHeights(treeModel)         0.996   357898     429449   1.2      0.1839      good
scale(errorModel.ageRate)                         0.615   357844     684979   1.91     0.2616      good
scale(constant.popSize)                           0.631   357294     7663     0.02     0.2435      good
scale(treeModel.rootHeight)                       0.644   358109     45734    0.13     0.2387      good
uniform(nodeHeights(treeModel))                           3579577    932944   0.26     0.7178      high
subtreeSlide(treeModel)                           12118.431787556    280311   0.16     0.1591      good
Narrow Exchange(treeModel)                                1789693    248656   0.14     0.3825      good
Wide Exchange(treeModel)                                  357376     13339    0.04     0.0063      low
wilsonBalding(treeModel)                                  357423     104544   0.29     0.018       slightly low
scale(kappa)                                      0.339   119624     228792   1.91     0.2783      good
frequencies                                       0.081   119396     228543   1.91     0.2664      good



2.227766388888889 hours

burnIn   = -1
maxState = 10000000

statistic                mean          hpdLower      hpdUpper      ESS
posterior                -5642.41      -5676.88      -5579.19      156.182
prior                    -2369.96      -2411.87      -2316.7       86.1602       *
likelihood               -3272.46      -3281.91      -3247.88      249.224
clock.rate               1.42656E-7    1.06654E-7    1.82082E-7    227.675
errorModel.ageRate       7.12553E-8    5.91987E-8    8.39501E-8    1273.19
treeModel.rootHeight     1.35663E5     97420.0       1.75097E5     138.118
constant.popSize         94500.1       63793.2       1.15728E5     126.347
kappa                    11.0776       6.40144       16.5056       6722.02
treeLikelihood           -3272.46      -3281.91      -3247.88      249.224
coalescent               -2356.18      -2397.88      -2303.19      86.1716       *

 * WARNING: The results of this MCMC analysis may be invalid as
            one or more statistics had very low effective sample sizes (ESS)
E[clock.rate]=1.5E-7
WARNING: 1.42656E-7     +- 1.3489E-9
E[errorModel.ageRate]=7E-8
WARNING: 7.12553E-8     +- 2.16332E-10


****************************** Java Result:

Operator analysis
Operator                                          Tuning   Count      Time     Time/Op  Pr(accept)  Performance suggestion
scale(kappa)                                      0.384   12009      59970    4.99     0.3266      good
scale(rate)                                       0.599   36245      181150   5.0      0.2856      good
up:rate down:nodeHeights(null)                    0.993   36147      101812   2.82     0.0962      slightly low	Try setting scaleFactor to about 0.9965
scale(populationSize)                             0.627   36263      1572     0.04     0.2424      good
scale(ageRelatedErrorRate)                        0.635   36062      180121   4.99     0.2823      good
scale(treeModel.rootHeight)                       0.679   36251      8469     0.23     0.2907      good
uniform(nodeHeights(null))                                361046     212919   0.59     0.7182      high
subtreeSlide(null)                                24252.175181239     58974    0.33     0.0902      slightly low	Try decreasing size to about 12126.087614484664
Narrow Exchange(null)                                     180508     56685    0.31     0.3857      good
Wide Exchange(null)                                       35991      3409     0.09     0.0086      low
wilsonBalding(null)                                       36104      24678    0.68     0.0222      slightly low

java.lang.NullPointerException
	at dr.util.NumberFormatter.formatToFieldWidth(NumberFormatter.java:92)
	at dr.inference.mcmc.MCMC.formattedOperatorName(MCMC.java:313)
	at dr.inference.mcmc.MCMC.showOperatorAnalysis(MCMC.java:299)
	at dr.inference.mcmc.MCMC.access$300(MCMC.java:50)
	at dr.inference.mcmc.MCMC$1.finished(MCMC.java:241)
	at dr.inference.markovchain.MarkovChain.fireFinished(MarkovChain.java:533)
	at dr.inference.markovchain.MarkovChain.terminateChain(MarkovChain.java:350)
	at dr.inference.mcmc.MCMC.chain(MCMC.java:194)
	at dr.inference.mcmc.MCMC.run(MCMC.java:152)
	at test.dr.evomodel.treelikelihood.PMDTest.testPMD(PMDTest.java:201)
	at sun.reflect.NativeMethodAccessorImpl.invoke0(Native Method)
	at sun.reflect.NativeMethodAccessorImpl.invoke(NativeMethodAccessorImpl.java:39)
	at sun.reflect.DelegatingMethodAccessorImpl.invoke(DelegatingMethodAccessorImpl.java:25)
	at com.intellij.junit3.JUnit3IdeaTestRunner.doRun(JUnit3IdeaTestRunner.java:108)
	at com.intellij.junit3.JUnit3IdeaTestRunner.startRunnerWithArgs(JUnit3IdeaTestRunner.java:42)
	at com.intellij.rt.execution.junit.JUnitStarter.prepareStreamsAndStart(JUnitStarter.java:165)
	at com.intellij.rt.execution.junit.JUnitStarter.main(JUnitStarter.java:60)
	at sun.reflect.NativeMethodAccessorImpl.invoke0(Native Method)
	at sun.reflect.NativeMethodAccessorImpl.invoke(NativeMethodAccessorImpl.java:39)
	at com.intellij.rt.execution.application.AppMain.main(AppMain.java:110)


Process finished with exit code -1


*/