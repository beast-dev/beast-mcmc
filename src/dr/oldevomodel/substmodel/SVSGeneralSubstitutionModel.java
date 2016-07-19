/*
 * SVSGeneralSubstitutionModel.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.oldevomodel.substmodel;

import dr.evolution.datatype.*;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.*;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.*;

/**
 * <b>A general model of sequence substitution with stochastic variable selection</b>. A general reversible class for any
 * data type.
 *
 * @author Marc Suchard
 * @version $Id: SVSGeneralSubstitutionModel.java,v 1.37 2006/05/05 03:05:10 msuchard Exp $
 */

public class SVSGeneralSubstitutionModel extends GeneralSubstitutionModel implements Likelihood,
        BayesianStochasticSearchVariableSelection, Citable {


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
        return !updateMatrix || BayesianStochasticSearchVariableSelection.Utils.connectedAndWellConditioned(probability,this);
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
            if (!BayesianStochasticSearchVariableSelection.Utils.connectedAndWellConditioned(probability,this)) {
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
    public Set<Likelihood> getLikelihoodSet() {
        return new HashSet<Likelihood>(Arrays.asList(this));
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

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.SUBSTITUTION_MODELS;
    }

    @Override
    public String getDescription() {
        return "Stochastic search variable selection, reversible substitution model";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CommonCitations.LEMEY_2009_BAYESIAN);
    }
}
