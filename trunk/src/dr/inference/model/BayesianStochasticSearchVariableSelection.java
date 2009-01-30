package dr.inference.model;

/**
 * @author Marc Suchard
 */

public interface BayesianStochasticSearchVariableSelection {

    public Parameter getIndicators();

    public boolean validState();

}
