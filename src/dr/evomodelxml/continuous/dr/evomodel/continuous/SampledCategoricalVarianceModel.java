package dr.evomodelxml.continuous.dr.evomodel.continuous;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.math.distributions.NormalDistribution;

public class SampledCategoricalVarianceModel extends AbstractModelLikelihood
        implements GradientWrtParameterProvider {

    private final Parameter traits;
    private final Parameter means;
    private final Parameter variances;
    private final Parameter assignments;

    private boolean likelihoodKnown;
    private boolean storedLikelihoodKnown;
    private double logLikelihood;
    private double storedLogLikelihood;

    private final GradientWrtParameterProvider gradientProvider;

    public SampledCategoricalVarianceModel(String name,
                                           Parameter traits,
                                           Parameter means,
                                           Parameter variances,
                                           Parameter assignments) {
        super(name);

        this.traits = traits;
        this.means = means;
        this.variances = variances;
        this.assignments = assignments;

        addVariable(traits);
        addVariable(means);
        addVariable(variances);
        addVariable(assignments);

        checkParameters();

        this.gradientProvider = makeGradientProvider();
    }

    GradientWrtParameterProvider makeGradientProvider() {

        final Likelihood likelihood = this;
        return new GradientWrtParameterProvider() {
            @Override
            public Likelihood getLikelihood() {
                return likelihood;
            }

            @Override
            public Parameter getParameter() {
                return variances;
            }

            @Override
            public int getDimension() {
                return variances.getDimension();
            }

            @Override
            public double[] getGradientLogDensity() {

                double[] gradient = new double[getDimension()];
                for (int i = 0; i < traits.getDimension(); ++i) {
                    int assignment = getAssignment(i);
                    double sd = getSd(i);

                    gradient[assignment] += NormalDistribution.gradLogPdf(traits.getParameterValue(i),
                            means.getParameterValue(i), sd);
                }

                return gradient;
            }
        };
    }

    void checkParameters() {

        if (traits.getDimension() != means.getDimension() ||
                traits.getDimension() != assignments.getDimension()) {
            throw new IllegalArgumentException("Unequal dimensions");
        }

        double[] assignments = this.assignments.getParameterValues();
        int min = Integer.MAX_VALUE;
        int max = Integer.MIN_VALUE;
        for (double x : assignments) {
            min = Math.min(min, (int) x);
            max = Math.max(max, (int) x);
        }

        if (min < 0 || max > variances.getDimension() - 1) {
            throw new IllegalArgumentException("Bad assignments");
        }
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // Do nothing
    }

    @Override
    protected void storeState() {
        storedLogLikelihood = logLikelihood;
        storedLikelihoodKnown = likelihoodKnown;
    }

    @Override
    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = storedLikelihoodKnown;
    }

    @Override
    protected void acceptState() { }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == traits) {
            likelihoodKnown = false;
        } else if (variable == variances) {
            likelihoodKnown = false;
        } else {
            throw new RuntimeException("Unknown variable");
        }
    }

    @Override
    public Model getModel() { return this; }

    @Override
    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = computeLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    private int getAssignment(int index) {
        return (int) assignments.getParameterValue(index);
    }

    private double getSd(int index) {
        double variance = variances.getParameterValue(getAssignment(index));
        return Math.sqrt(variance);
    }

    private double computeLogLikelihood() {

        double logLikelihood = 0;
        for (int i = 0; i < traits.getDimension(); ++i) {

            double sd = getSd(i);
            logLikelihood += NormalDistribution.logPdf(traits.getParameterValue(i), means.getParameterValue(i), sd);
        }

        return logLikelihood;
    }



    @Override
    public void makeDirty() {
        likelihoodKnown = false;
    }

    @Override
    public Likelihood getLikelihood() {
        return gradientProvider.getLikelihood();
    }

    @Override
    public Parameter getParameter() {
        return gradientProvider.getParameter();
    }

    @Override
    public int getDimension() {
        return gradientProvider.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        return gradientProvider.getGradientLogDensity();
    }
}
