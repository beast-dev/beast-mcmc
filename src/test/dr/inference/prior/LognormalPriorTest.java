package test.dr.inference.prior;

import dr.evolution.util.Units;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.loggers.ArrayLogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.DummyLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.ScaleOperator;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.inference.trace.ArrayTraceList;
import dr.inference.trace.Trace;
import dr.inference.trace.TraceCorrelation;
import dr.math.MathUtils;
import dr.math.distributions.LogNormalDistribution;
import junit.framework.Test;
import junit.framework.TestSuite;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Walter Xie
 * convert testLognormalPrior.xml in the folder /example
 */

public class LognormalPriorTest extends TraceCorrelationAssert {

    public LognormalPriorTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        MathUtils.setSeed(666);
    }


    public void testLognormalPrior() {
//        ConstantPopulation constant = new ConstantPopulation(Units.Type.YEARS);
//        constant.setN0(popSize); // popSize
        Parameter popSize = new Parameter.Default(6.0);
        popSize.setId(ConstantPopulationModelParser.POPULATION_SIZE);
        ConstantPopulationModel demo = new ConstantPopulationModel(popSize, Units.Type.YEARS);

        //Likelihood
        Likelihood dummyLikelihood = new DummyLikelihood(demo);

        // Operators
        OperatorSchedule schedule = new SimpleOperatorSchedule();

        MCMCOperator operator = new ScaleOperator(popSize, 0.75);
        operator.setWeight(1.0);
        schedule.addOperator(operator);

        // Log
        ArrayLogFormatter formatter = new ArrayLogFormatter(false);

        MCLogger[] loggers = new MCLogger[2];
        loggers[0] = new MCLogger(formatter, 1000, false);
//        loggers[0].add(treeLikelihood);
        loggers[0].add(popSize);

        loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), 100000, false);
//        loggers[1].add(treeLikelihood);
        loggers[1].add(popSize);

        // MCMC
        MCMC mcmc = new MCMC("mcmc1");
        MCMCOptions options = new MCMCOptions();
        options.setChainLength(1000000);
        options.setUseCoercion(true); // autoOptimize = true
        options.setCoercionDelay(100);
        options.setTemperature(1.0);
        options.setFullEvaluationCount(2000);

        DistributionLikelihood logNormalLikelihood = new DistributionLikelihood(new LogNormalDistribution(1.0, 1.0), 0); // meanInRealSpace="false"
        logNormalLikelihood.addData(popSize);

        List<Likelihood> likelihoods = new ArrayList<Likelihood>();
        likelihoods.add(logNormalLikelihood);
        Likelihood prior = new CompoundLikelihood(0, likelihoods);

        likelihoods.clear();
        likelihoods.add(dummyLikelihood);
        Likelihood likelihood = new CompoundLikelihood(-1, likelihoods);

        likelihoods.clear();
        likelihoods.add(prior);
        likelihoods.add(likelihood);
        Likelihood posterior = new CompoundLikelihood(0, likelihoods);

        mcmc.setShowOperatorAnalysis(true);
        mcmc.init(options, posterior, schedule, loggers);
        mcmc.run();

        // time
        System.out.println(mcmc.getTimer().toString());

        // Tracer
        List<Trace> traces = formatter.getTraces();
        ArrayTraceList traceList = new ArrayTraceList("LognormalPriorTest", traces, 0);

        for (int i = 1; i < traces.size(); i++) {
            traceList.analyseTrace(i);
        }

//      <expectation name="param" value="4.48168907"/>

        TraceCorrelation popSizeStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(ConstantPopulationModelParser.POPULATION_SIZE));

        System.out.println("Expectation of Log-Normal(1,1) is e^(M+S^2/2) = e^(1.5) = " + Math.exp(1.5));
        assertExpectation(ConstantPopulationModelParser.POPULATION_SIZE, popSizeStats, Math.exp(1.5));
    }

    public static Test suite() {
        return new TestSuite(LognormalPriorTest.class);
    }
}
