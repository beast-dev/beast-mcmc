package dr.evomodel.coalescent;

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodelxml.coalescent.BNPRSamplingLikelihoodParser;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.MaxIterationsExceededException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.integration.RombergIntegrator;

import java.util.Arrays;

/**
 * Created by mkarcher on 3/9/17.
 */
public class BNPRSamplingLikelihood extends AbstractModelLikelihood {
    private boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown = false;
    private double logLikelihood;
    private double storedLogLikelihood;

    private boolean samplingTimesKnown = false;
    private boolean storedSamplingTimesKnown = false;
    private double[] samplingTimes;
    private double[] storedSamplingTimes;
    private double[] logPopSizes;
    private double[] storedLogPopSizes;

    private Tree tree;
    private int numSamples;
    private Parameter betas;
    private DemographicModel population;
    private UnivariateRealFunction[] covariates;
    private DemographicWrapper demo;
    private RombergIntegrator integrator;

    public BNPRSamplingLikelihood(Tree tree, Parameter betas, DemographicModel population) {
        this(tree, betas, population, null);
    }

    public BNPRSamplingLikelihood(Tree tree, Parameter betas, DemographicModel population, UnivariateRealFunction[] covariates) {
        super(BNPRSamplingLikelihoodParser.SAMPLING_LIKELIHOOD);

        this.tree = tree;
        this.betas = betas;
        this.population = population;
        this.covariates = covariates;

        this.likelihoodKnown = false;
        this.samplingTimesKnown = false;

        this.demo = new DemographicWrapper(population.getDemographicFunction());
        this.integrator = new RombergIntegrator();
    }

    private void setupSamplingTimes() {
        this.numSamples = tree.getExternalNodeCount();
        samplingTimes = new double[numSamples];
        logPopSizes = new double[numSamples];

        for (int i = 0; i < numSamples; i++) {
            NodeRef node = tree.getExternalNode(i);
            samplingTimes[i] = tree.getNodeHeight(node);
        }

        Arrays.sort(samplingTimes);

        samplingTimesKnown = true;
    }

    public double[] getSamplingTimes() {
        if (!samplingTimesKnown) {
            setupSamplingTimes();
        }

        return samplingTimes;
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        if (!samplingTimesKnown) {
            setupSamplingTimes();
        }

        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }

        return logLikelihood;
    }

    public double calculateLogLikelihood() {
        double minSample = samplingTimes[0];
        double maxSample = samplingTimes[numSamples-1];

        double beta0 = betas.getParameterValue(0);
        double beta1 = betas.getParameterValue(1);

        double logLik = numSamples * beta0; // Start with beta0 part of likelihood.
        System.out.println(logLik);
        DemographicFunction f = population.getDemographicFunction();

        demo.setF(f); // Possibly unnecessary if the DemographicFunction reference stays the same after model changes
        PowerFunctional pf = new PowerFunctional(demo, beta1);

        // Gather log-population sizes to facilitate likelihood calculation.
        // Separated to make cut-paste into a new setupLogPopSizes() function easy.
        for (int i = 0; i < numSamples; i++) {
            logPopSizes[i] = f.getLogDemographic(samplingTimes[i]);
        }

        double eventLogLik = 0;
        // Calculate the event component of the inhomogeneous Poisson process log-likelihood
        for (int i = 0; i < numSamples; i++) {
            eventLogLik += beta1 * logPopSizes[i];
        }

        System.out.println(eventLogLik);
        logLik += eventLogLik;

        double integralLogLik = 0;
        // Calculate the integral component of the inhomogeneous Poisson process log-likelihood
        try {
            integralLogLik -= Math.exp(beta0) * integrator.integrate(pf, minSample, maxSample);
        } catch (MaxIterationsExceededException e) {
            throw new RuntimeException(e);
        } catch (FunctionEvaluationException e) {
            throw new RuntimeException(e);
        }

        System.out.println(integralLogLik);
        logLik += integralLogLik;

        return logLik;
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
        samplingTimesKnown = false;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        makeDirty();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        makeDirty();
    }

    @Override
    protected void storeState() {
        // copy the sampling times into the storedSamplingTimes
        System.arraycopy(samplingTimes, 0, storedSamplingTimes, 0, samplingTimes.length);

        storedSamplingTimesKnown = samplingTimesKnown;
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;

        System.arraycopy(logPopSizes, 0, storedLogPopSizes, 0, logPopSizes.length);
    }

    @Override
    protected void restoreState() {
        // copy the sampling times into the storedSamplingTimes
        System.arraycopy(storedSamplingTimes, 0, samplingTimes, 0, storedSamplingTimes.length);

        samplingTimesKnown = storedSamplingTimesKnown;
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;

        System.arraycopy(storedLogPopSizes, 0, logPopSizes, 0, storedLogPopSizes.length);
    }

    @Override
    protected void acceptState() {
        // nothing to do
    }

    public String toString() {
        return Double.toString(logLikelihood);
    }

    private class PowerFunctional implements UnivariateRealFunction {
        private double power;
        private UnivariateRealFunction f;

        public PowerFunctional(UnivariateRealFunction f, double power) {
            this.f = f;
            this.power = power;
        }

        @Override
        public double value(double t) throws FunctionEvaluationException {
            return Math.pow(f.value(t), power);
        }

        public double getPower() {
            return power;
        }

        public void setPower(double power) {
            this.power = power;
        }

        public UnivariateRealFunction getF() {
            return f;
        }

        public void setF(UnivariateRealFunction f) {
            this.f = f;
        }
    }

    private class DemographicWrapper implements UnivariateRealFunction {
        private DemographicFunction f;

        public DemographicWrapper(DemographicFunction f) {
            this.f = f;
        }

        public double value(double t) {
            return f.getDemographic(t);
        }

        public void setF(DemographicFunction f) {
            this.f = f;
        }
    }

    public static void main(String[] args) throws Exception {
        NewickImporter importer = new NewickImporter("((((5:0.5,1:0.2):0.5,0:1):0.2,2:0.8):0.2,3:1.4)");
        Tree tree = importer.importNextTree();
        Parameter betas = new Parameter.Default(new double[] {3.0, 1.0});
        Parameter N0 = new Parameter.Default(new double[] {3.5, 8.9, 1.0, 3.4, 1.0});
        double[] epochLengths = new double[] {0.35, 0.35, 0.35, 0.35};

        int numGridPoints = 5;

        DemographicModel Ne = new PiecewisePopulationModel("Ne(t)", N0, epochLengths, false, Units.Type.DAYS);

        System.out.println("Got here A");

        BNPRSamplingLikelihood samplingLikelihood = new BNPRSamplingLikelihood(tree, betas, Ne);

        System.out.println("Sampling times: " + Arrays.toString(samplingLikelihood.getSamplingTimes()));

        System.out.println("Sampling likelihood: " + samplingLikelihood.getLogLikelihood());
    }
}
