package dr.inference.glm;

import dr.inference.distribution.LogisticRegression;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

/**
 * Analytic gradient provider for logistic regression log-likelihood with respect to a fixed-effect parameter.
 */
public abstract class AbstractLogisticRegressionGradient implements GradientWrtParameterProvider, Reportable {

    protected final LogisticRegression likelihood;
    protected final Parameter parameter;
    protected final Parameter dependent;

    public AbstractLogisticRegressionGradient(LogisticRegression likelihood, Parameter parameter) {
        this.likelihood = likelihood;
        this.parameter = parameter;
        this.dependent = likelihood.getDependentVariable();
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        return parameter.getDimension();
    }

    protected static double logistic(double x) {
        if (x >= 0) {
            double expNegX = Math.exp(-x);
            return 1.0 / (1.0 + expNegX);
        } else {
            double expX = Math.exp(x);
            return expX / (1.0 + expX);
        }
    }

    @Override
    public String getReport() {
        double tolerance = getTolerance();
        return GradientWrtParameterProvider.getReportAndCheckForError(this, 0, Double.POSITIVE_INFINITY, tolerance);
    }

    double getTolerance() {
        return 1E-3;
    }
}
