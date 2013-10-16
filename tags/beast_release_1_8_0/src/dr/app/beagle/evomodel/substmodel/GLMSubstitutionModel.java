/*
 * GLMSubstitutionModel.java
 *
 * Copyright (C) 2002-2012 Alexei Drummond, Andrew Rambaut & Marc A. Suchard
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

package dr.app.beagle.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.inference.distribution.LogLinearModel;
import dr.inference.loggers.LogColumn;
import dr.inference.model.BayesianStochasticSearchVariableSelection;
import dr.inference.model.Model;

/**
 * @author Marc A. Suchard
 */
public class GLMSubstitutionModel extends ComplexSubstitutionModel {

    public GLMSubstitutionModel(String name, DataType dataType, FrequencyModel rootFreqModel,
                                LogLinearModel glm) {

        super(name, dataType, rootFreqModel, null);
        this.glm = glm;
        addModel(glm);
        testProbabilities = new double[stateCount * stateCount];

    }

    protected void setupRelativeRates(double[] rates) {
        System.arraycopy(glm.getXBeta(),0,rates,0,rates.length);       
    }


    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == glm) {
            updateMatrix = true;
            fireModelChanged();
        } else
            super.handleModelChangedEvent(model, object, index);
    }

    public LogColumn[] getColumns() {
        return glm.getColumns();
    }

    public double getLogLikelihood() {
        double logL = super.getLogLikelihood();
        if (logL == 0 &&
                BayesianStochasticSearchVariableSelection.Utils.connectedAndWellConditioned(testProbabilities, this)) {
            // Also check that graph is connected
            return 0;
        }
        return Double.NEGATIVE_INFINITY;
    }


    private LogLinearModel glm;
    private double[] testProbabilities;

}
