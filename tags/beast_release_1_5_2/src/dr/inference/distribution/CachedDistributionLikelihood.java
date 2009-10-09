package dr.inference.distribution;

import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * @author Marc Suchard
 */
public class CachedDistributionLikelihood extends AbstractModelLikelihood {

    public CachedDistributionLikelihood(String name, AbstractDistributionLikelihood likelihood, Variable variable) {
        super(name);
        this.likelihood = likelihood;
        addVariable(variable);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        likelihoodKnown = false;
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        likelihoodKnown = false;
    }

    protected void storeState() {
        storedLikelihoodKnow = likelihoodKnown;
        storedLogLikelihood = logLikelihood;

    }

    protected void restoreState() {
        likelihoodKnown = storedLikelihoodKnow;
        logLikelihood = storedLogLikelihood;

    }

    protected void acceptState() {

    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    private double calculateLogLikelihood() {
        return likelihood.calculateLogLikelihood();
    }

    public void makeDirty() {
        likelihoodKnown = false;
    }

    private boolean likelihoodKnown;
    private boolean storedLikelihoodKnow;
    private double logLikelihood;
    private double storedLogLikelihood;

    private final AbstractDistributionLikelihood likelihood;
}
