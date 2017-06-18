package dr.evomodel.coalescent;

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeModel;
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

    private double[] epochWidths = null;
    private double[] midpoints = null;

    private Tree tree;
    private int numSamples;
    private Parameter betas;
    private DemographicModel population;
    private UnivariateRealFunction[] covariates;
    private DemographicWrapper demo;
    private DemographicModel demoModel = null;
    private RombergIntegrator integrator;

    public BNPRSamplingLikelihood(Tree tree, Parameter betas, DemographicModel population) {
        this(tree, betas, population, null, null);
    }

    public BNPRSamplingLikelihood(Tree tree, Parameter betas, DemographicModel population, double[] epochWidths) {
        this(tree, betas, population, epochWidths, null);
    }

    public BNPRSamplingLikelihood(Tree tree, Parameter betas, DemographicModel population, double[] epochWidths,
                                  UnivariateRealFunction[] covariates) {
        super(BNPRSamplingLikelihoodParser.SAMPLING_LIKELIHOOD);

        this.tree = tree;
        this.betas = betas;
        this.population = population;
        this.covariates = covariates;

        this.likelihoodKnown = false;
        this.samplingTimesKnown = false;

        this.demo = new DemographicWrapper(population.getDemographicFunction());
        this.demoModel = population;
        this.integrator = new RombergIntegrator();

        this.numSamples = tree.getExternalNodeCount();
        samplingTimes = new double[numSamples];
        storedSamplingTimes = new double[numSamples];
        logPopSizes = new double[numSamples];
        storedLogPopSizes = new double[numSamples];

        if (epochWidths != null) {
            setEpochs(epochWidths);
        }

        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }

        addModel(demoModel);

        addVariable(this.betas);

        setupSamplingTimes();
    }

    public void setEpochs(double[] epochWidths) {
        this.epochWidths = epochWidths;

        if (epochWidths == null) {
            this.midpoints = null;
        } else {
            this.midpoints = new double[epochWidths.length + 1];
            double start = 0;
            for (int i = 0; i < epochWidths.length; i++) {
                this.midpoints[i] = start + epochWidths[i] / 2.0;
                start += epochWidths[i];
            }
            this.midpoints[epochWidths.length] = start + epochWidths[epochWidths.length - 1] / 2.0;
        }
    }

    private void setupSamplingTimes() {
//        this.numSamples = tree.getExternalNodeCount();
//        samplingTimes = new double[numSamples];
//        logPopSizes = new double[numSamples];

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

//        System.out.printf("beta0 part = %f\n", logLik);

        DemographicFunction f = population.getDemographicFunction();

        // Gather log-population sizes to facilitate likelihood calculation.
        // Separated to make cut-paste into a new setupLogPopSizes() function easy.
        for (int i = 0; i < numSamples; i++) {
            logPopSizes[i] = f.getLogDemographic(samplingTimes[i]);
        }

        double eventLogLik = 0.0;
        // Calculate the event component of the inhomogeneous Poisson process log-likelihood
        for (int i = 0; i < numSamples; i++) {
            eventLogLik += beta1 * logPopSizes[i];
        }

//        System.out.printf("event part = %f\n", eventLogLik);

        logLik += eventLogLik;

        double integralLogLik = 0.0;

        if (epochWidths == null) {
            demo.setF(f); // Possibly unnecessary if the DemographicFunction reference stays the same after model changes
            PowerFunctional pf = new PowerFunctional(demo, beta1);

            // Calculate the integral component of the inhomogeneous Poisson process log-likelihood
            try {
                integralLogLik -= Math.exp(beta0) * integrator.integrate(pf, minSample, maxSample);
            } catch (MaxIterationsExceededException e) {
                throw new RuntimeException(e);
            } catch (FunctionEvaluationException e) {
                throw new RuntimeException(e);
            }
        } else {
            double[] fBeta = new double[this.midpoints.length];
            for (int i = 0; i < fBeta.length; i++) {
                fBeta[i] = Math.exp(f.getLogDemographic(this.midpoints[i]) * beta1);
            }
            integralLogLik -= Math.exp(beta0) * integratePiecewise(this.epochWidths, fBeta, minSample, maxSample);
        }

//        System.out.printf("integral part = %f\n", integralLogLik);
        logLik += integralLogLik;

        return logLik;
    }

    private static double integratePiecewise(double[] widths, double[] heights, double start, double end) {
        double[] widths2 = new double[widths.length + 1];
        System.arraycopy(widths, 0, widths2, 0, widths.length);
        widths2[widths.length] = Double.POSITIVE_INFINITY;

        double result = 0.0;
        int i = 0;
        while (start > widths2[i]) {
            start -= widths2[i];
            end -= widths2[i];
            i++;
        }

        if (end < widths2[i]) {
            result = (end - start) * heights[i];
        } else {
            result += (widths2[i] - start) * heights[i];
            end -= widths2[i];
            i++;

            while (end > widths2[i]) {
                result += widths2[i] * heights[i];
                end -= widths2[i];
                i++;
            }

            if (end > 0.0) {
                result += end * heights[i];
            }
        }

        return result;
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
//        DistributionLikelihood lik = new DistributionLikelihood(new NormalDistribution(0, 1));
//        Parameter par = new Parameter.Default(new double[]{0.0, 1.0});
//        lik.addData(par);
//        System.out.println("Likelihood = " + lik.getLogLikelihood());

        NewickImporter importer = new NewickImporter("((((5:0.5,1:0.2):0.5,0:1):0.2,2:0.8):0.2,3:1.4)");
        Tree tree = importer.importNextTree();
        Parameter betas = new Parameter.Default(new double[] {3.0, 1.0});
//        Parameter N0 = new Parameter.Default(new double[] {3.0, 3.59, 2.70, 10.0, 100.0});
//        double[] epochLengths = new double[] {0.35, 0.35, 0.35, 0.35};
        Parameter N0 = new Parameter.Default(new double[] {1.0, 0.5, 1.0, 1.0, 10.0});
        double[] epochLengths = new double[] {0.35, 0.35, 0.35, 0.35};

//        int numGridPoints = 2;

        DemographicModel Ne = new PiecewisePopulationModel("Ne(t)", N0, epochLengths, false, Units.Type.DAYS);

        System.out.println("Integration: " + integratePiecewise(epochLengths, N0.getParameterValues(), 3.0, 14.0));

        System.out.println("Got here A");

        BNPRSamplingLikelihood samplingLikelihood1 = new BNPRSamplingLikelihood(tree, betas, Ne);
        BNPRSamplingLikelihood samplingLikelihood2 = new BNPRSamplingLikelihood(tree, betas, Ne, epochLengths);

//        System.out.println("Sampling times: " + Arrays.toString(samplingLikelihood.getSamplingTimes()));
//
        System.out.println("Sampling likelihood 1: " + samplingLikelihood1.getLogLikelihood());
//
//        samplingLikelihood.makeDirty();
//
        System.out.println("Sampling likelihood 2: " + samplingLikelihood2.getLogLikelihood());
//
//        betas.setParameterValue(0, 0.0);
//
//        System.out.println("Sampling likelihood: " + samplingLikelihood.getLogLikelihood());
//
//        samplingLikelihood.makeDirty();
//
//        System.out.println("Sampling likelihood: " + samplingLikelihood.getLogLikelihood());
//
//        System.out.println("Midpoints: " + Arrays.toString(samplingLikelihood.midpoints));
    }
}
