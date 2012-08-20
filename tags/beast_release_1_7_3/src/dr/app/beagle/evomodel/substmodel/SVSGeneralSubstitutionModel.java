package dr.app.beagle.evomodel.substmodel;

import dr.inference.model.*;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.evolution.datatype.DataType;

/**
 * @author Marc Suchard
 */

public class SVSGeneralSubstitutionModel extends GeneralSubstitutionModel implements Likelihood,
        BayesianStochasticSearchVariableSelection {

    public SVSGeneralSubstitutionModel(String name, DataType dataType, FrequencyModel freqModel,
                                       Parameter ratesParameter, Parameter indicatorsParameter) {
        super(name, dataType, freqModel, ratesParameter, -1);

        if (indicatorsParameter == null) {
            this.indicatorsParameter = new Parameter.Default(ratesParameter.getDimension(), 1.0);
        } else {
            this.indicatorsParameter  = indicatorsParameter;
            addVariable(indicatorsParameter);
        }

    }

    @Override
    protected void setupRelativeRates(double[] rates) {
        for (int i = 0; i < rates.length; i++) {
            rates[i] = ratesParameter.getParameterValue(i) * indicatorsParameter.getParameterValue(i);
        }
    }

    public Parameter getIndicators() {
        return indicatorsParameter;
    }

    public boolean validState() {
        return !updateMatrix || Utils.connectedAndWellConditioned(probability,this);
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == ratesParameter && indicatorsParameter.getParameterValue(index) == 0)
            return; // Does not affect likelihood
        super.handleVariableChangedEvent(variable,index,type);
    }

    /**
     * Get the model.
     *
     * @return the model.
     */
    public Model getModel() {
        return this;
    }

    /**
     * Get the log likelihood.
     *
     * @return the log likelihood.
     */
    public double getLogLikelihood() {
        if (updateMatrix) {
            if (!Utils.connectedAndWellConditioned(probability,this)) {
                return Double.NEGATIVE_INFINITY;
            }
        }
        return 0;
    }

    /**
     * Needs to be evaluated before the corresponding data likelihood.
     * @return
     */
    public boolean evaluateEarly() {
        return true;
    }

    /**
     * Forces a complete recalculation of the likelihood next time getLikelihood is called
     */
    public void makeDirty() {
       updateMatrix = true;
    }

    /**
     * @return A detailed name of likelihood for debugging.
     */
    public String prettyName() {
        return "SVSGeneralSubstitutionModel-connectedness";
    }

    @Override
    public boolean isUsed() {
        return super.isUsed() && isUsed;
    }

    public void setUsed() {
        isUsed = true;
    }

    private boolean isUsed = false;

    // **************************************************************
    // Loggable IMPLEMENTATION
    // **************************************************************

    public LogColumn[] getColumns() {
        return new LogColumn[]{
                new LikelihoodColumn(getId())
        };
    }

    protected class LikelihoodColumn extends NumberColumn {
        public LikelihoodColumn(String label) {
            super(label);
        }

        public double getDoubleValue() {
            return getLogLikelihood();
        }
    }

    private double[] probability = null;

    private final Parameter indicatorsParameter;


}
