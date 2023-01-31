package dr.inference.operators.hmc;

/**
 * @author Marc A. Suchard
 */
class StepSize {

    final private Options options;
    final private double mu;

    private double stepSize;
    private double logStepSize;
    private double averageLogStepSize;
    private double h;

    StepSize(double initialStepSize) {
        this(initialStepSize, new Options());
    }

    private StepSize(double initialStepSize, Options options) {
        this.options = options;
        this.mu = Math.log(options.muFactor * initialStepSize);
        this.stepSize = initialStepSize;
        this.logStepSize = Math.log(stepSize);
        this.averageLogStepSize = 0;
        this.h = 0;
    }

    void update(long m, double cumAcceptProb, double numAcceptProbStates) {

        if (m <= options.adaptLength) {

            h = (1 - 1 / (m + options.t0)) * h + 1 / (m + options.t0) * (options.targetAcceptRate - (cumAcceptProb / numAcceptProbStates));
            logStepSize = mu - Math.sqrt(m) / options.gamma * h;
            averageLogStepSize = Math.pow(m, -options.kappa) * logStepSize +
                    (1 - Math.pow(m, -options.kappa)) * averageLogStepSize;
            stepSize = Math.exp(logStepSize);
        }
    }

    double getStepSize() {
        return stepSize;
    }

    static class Options { //TODO: these values might be adjusted for dual averaging.
        double kappa = 0.75;
        double t0 = 10.0;
        double gamma = 0.05;
        double targetAcceptRate = 0.85;
        double muFactor = 10.0;
        int adaptLength = 1000;
    }
}
