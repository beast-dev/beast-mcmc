package dr.app.beagle.evomodel.substmodel;

import dr.inference.model.BayesianStochasticSearchVariableSelection;
import dr.inference.model.Parameter;
import dr.inference.model.Model;
import dr.inference.model.Likelihood;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.evolution.datatype.DataType;

/**
 * @author Marc Suchard
 */

public class SVSGeneralSubstitutionModel extends GeneralSubstitutionModel implements Likelihood,
        BayesianStochasticSearchVariableSelection {

    public SVSGeneralSubstitutionModel(DataType dataType, FrequencyModel freqModel,
                                       Parameter ratesParameter, Parameter indicatorsParameter) {
        super(dataType, freqModel, ratesParameter, -1);
        this.indicatorsParameter  = indicatorsParameter;
        addVariable(indicatorsParameter);
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
        return !updateMatrix || checkFullyConnected();
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
            if (!checkFullyConnected()) {
                System.err.println("SVSGeneralSubstitutionModel is not fully connected.");
                return Double.NEGATIVE_INFINITY;
            }
        }
        return 0;
    }

    private boolean checkFullyConnected() {
        if (probability == null)
            probability = new double[stateCount*stateCount];

        getTransitionProbabilities(1.0,probability);
        final int length = stateCount*stateCount;
        for(int i=0; i<length; i++) {
            if(probability[i] == 0)
                return false;
        }
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
