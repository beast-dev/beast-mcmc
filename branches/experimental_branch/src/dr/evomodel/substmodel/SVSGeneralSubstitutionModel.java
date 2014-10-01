package dr.evomodel.substmodel;

import dr.evolution.datatype.*;
import dr.evomodelxml.substmodel.GeneralSubstitutionModelParser;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.*;
import dr.xml.*;

/**
 * <b>A general model of sequence substitution with stochastic variable selection</b>. A general reversible class for any
 * data type.
 *
 * @author Marc Suchard
 * @version $Id: SVSGeneralSubstitutionModel.java,v 1.37 2006/05/05 03:05:10 msuchard Exp $
 */

public class SVSGeneralSubstitutionModel extends GeneralSubstitutionModel implements Likelihood,
        BayesianStochasticSearchVariableSelection {


    public SVSGeneralSubstitutionModel(DataType dataType, FrequencyModel freqModel, Parameter parameter,
                                       Parameter indicator) {
        super(dataType, freqModel, parameter, 1);

        if (indicator != null) {
            rateIndicator = indicator;
            addVariable(rateIndicator);
        } else {
            rateIndicator = new Parameter.Default(parameter.getDimension(), 1.0);
        }
    }

    protected SVSGeneralSubstitutionModel(String name, DataType dataType, FrequencyModel freqModel, int relativeTo) {
        super(name, dataType, freqModel, relativeTo);
    }

    public Parameter getIndicators() {
        return rateIndicator;
    }

    public boolean validState() {
        return !updateMatrix || Utils.connectedAndWellConditioned(probability,this);
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == ratesParameter && rateIndicator.getParameterValue(index) == 0)
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

    protected void setupRelativeRates() {

        for (int i = 0; i < relativeRates.length; i++) {
            relativeRates[i] = ratesParameter.getParameterValue(i) * rateIndicator.getParameterValue(i);
        }
    }

    void normalize(double[][] matrix, double[] pi) {
        double subst = 0.0;
        int dimension = pi.length;

        //final int dim = rateIndicator.getDimension();
        //int sum = 0;
        //for (int i = 0; i < dim; i++)
        //	sum += rateIndicator.getParameterValue(i);


        for (int i = 0; i < dimension; i++)
            subst += -matrix[i][i] * pi[i];

        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                matrix[i][j] = matrix[i][j] / subst; // / sum;
            }
        }
    }
    
    @Override
    public boolean isUsed() {
        return super.isUsed() && isUsed;
    }

    public void setUsed() {
        isUsed = true;
    }

    private boolean isUsed = false;

    private Parameter rateIndicator;
}
