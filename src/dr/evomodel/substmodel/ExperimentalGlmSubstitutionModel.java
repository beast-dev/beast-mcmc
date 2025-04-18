/*
 * ExperimentalGlmSubstitutionModel.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.inference.glm.ExperimentalGeneralizedLinearModel;
import dr.inference.loggers.LogColumn;
import dr.inference.model.BayesianStochasticSearchVariableSelection;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.*;

/**
 * @author Marc A. Suchard
 */
public class ExperimentalGlmSubstitutionModel extends ComplexSubstitutionModel {

    public ExperimentalGlmSubstitutionModel(String name, DataType dataType, FrequencyModel rootFreqModel,
                                            ExperimentalGeneralizedLinearModel glm) {

        super(name, dataType, rootFreqModel, null);
        this.glm = glm;
        addModel(glm);
        testProbabilities = new double[stateCount*stateCount];
            
    }

    public double[] getRates() {
        return glm.getXBeta();
    }


    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == glm) {
            updateMatrix = true;
            fireModelChanged();
        } else
            super.handleModelChangedEvent(model, object, index);
    }

    // This info can be gotten from the GLM
//    public LogColumn[] getColumns() {
//        return glm.getColumns();
//    }

    public LogColumn[] getColumns() {
        //Aggregate columns from ComplexSubstitutionModel with glm.columns
        LogColumn[] aggregated = new LogColumn[glm.getColumns().length + 2];
        int index = 0;
        for (LogColumn col : glm.getColumns()) {
            aggregated[index] = col;
            index++;
        }
        aggregated[index++] = new LikelihoodColumn(getId() + ".L");
        aggregated[index++] = new NormalizationColumn(getId() + ".Norm");

        return aggregated;
        //return glm.getColumns();
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

    public ExperimentalGeneralizedLinearModel getGeneralizedLinearModel() { return glm; }

    @Override
    public String getDescription() {
        return "Generalized linear (model, GLM) substitution model"; // TODO Horrible; fix
    }
    
    @Override
    public Set<Likelihood> getLikelihoodSet() {
        return new HashSet<Likelihood>(Arrays.asList(this, glm));
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CommonCitations.LEMEY_2014_UNIFYING);
    }

    private ExperimentalGeneralizedLinearModel glm;
    private double[] testProbabilities;    
}
