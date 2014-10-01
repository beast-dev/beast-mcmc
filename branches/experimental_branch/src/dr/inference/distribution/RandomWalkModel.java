package dr.inference.distribution;

import dr.inference.model.*;

import java.util.logging.Logger;

/**
 * @author Marc A. Suchard
 */
public class RandomWalkModel extends AbstractModelLikelihood {

    public RandomWalkModel(
            ParametricDistributionModel distribution,
            Parameter data, boolean forwardOrder, boolean logScale) {
        super(null);
        this.distribution = distribution;
        this.forwardOrder = forwardOrder;
        this.logScale = logScale;
        this.data = data;
        if (distribution != null) {
            addModel(distribution);
        }

        double lower = Double.NEGATIVE_INFINITY;
        if (logScale)
            lower = 0.0;

        addVariable(data);
        if (data instanceof CompoundParameter) {
            CompoundParameter cp = (CompoundParameter) data;
            for (int i = 0; i < cp.getParameterCount(); i++) {
                Parameter p = cp.getParameter(i);
                p.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, lower, p.getDimension()));
            }
        } else
            data.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, lower, data.getDimension()));

        Logger.getLogger("dr.inference").info("Setting up a first-order random walk:");
        Logger.getLogger("dr.inference").info("\tData parameter: " + data.getId());
        Logger.getLogger("dr.inference").info("\tOn scale: " + (logScale ? "log" : "real"));
        Logger.getLogger("dr.inference").info("\tDistribution: " + distribution.getId());
        Logger.getLogger("dr.inference").info("\tIf you publish results using this model, please cite Suchard and Lemey (in preparation)\n");
    }

    protected double calculateLogLikelihood() {

        final int dim = data.getDimension();

        double logLikelihood = 0;
        double previous = data.getParameterValue(0);
        if (logScale)
            previous = Math.log(previous);
        for (int i = 1; i < dim; i++) {
            double current = data.getParameterValue(i);
            if (logScale)
                current = Math.log(current);
            logLikelihood += distribution.logPdf(current - previous);
            if (logScale)
                logLikelihood -= current;
            previous = current;
        }

        return logLikelihood;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        likelihoodKnown = false;
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        likelihoodKnown = false;
    }

    protected void storeState() {
    }

    protected void restoreState() {
    }

    protected void acceptState() {
    }

    private final ParametricDistributionModel distribution;
    private final boolean logScale;
    private Parameter data;
    protected boolean likelihoodKnown;
    private boolean forwardOrder;

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        return calculateLogLikelihood();
    }

    public void makeDirty() {
    }
}
