/*
 * LogRateSubstitutionModel.java
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
import dr.inference.loggers.LogColumn;
import dr.inference.model.BayesianStochasticSearchVariableSelection;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.*;

/**
 * @author Filippo Monti
 * @author Marc A. Suchard
 */

public class LogRateSubstitutionModel extends ComplexSubstitutionModel { // implements ParameterReplaceableSubstitutionModel, DifferentiableSubstitutionModel

    public LogRateSubstitutionModel(String name, DataType dataType, FrequencyModel rootFreqModel,
                                    LogAdditiveCtmcRateProvider lrm) {

        super(name, dataType, rootFreqModel, null);
        this.lrm = lrm;
        addModel(lrm);
        testProbabilities = new double[stateCount*stateCount];

    }

    @Override
    public LogAdditiveCtmcRateProvider getRateProvider() {
        return lrm;
    }

    protected void setupRelativeRates(double[] rates) {
        double[] transformedRates = lrm.getRates();
        System.arraycopy(transformedRates,0,rates,0,rates.length);
    }

    @Override
    public Set<Likelihood> getLikelihoodSet() {
        return new HashSet<>(Arrays.asList(this, lrm));
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == lrm) {
            updateMatrix = true;
            fireModelChanged();
        } else {
            super.handleModelChangedEvent(model, object, index);
        }
    }

    public LogColumn[] getColumns() {
        //Aggregate columns from ComplexSubstitutionModel with glm.columns
        LogColumn[] aggregated = new LogColumn[lrm.getColumns().length + 2];
        int index = 0;
        for (LogColumn col : lrm.getColumns()) {
            aggregated[index] = col;
            index++;
        }
        aggregated[index++] = new LikelihoodColumn(getId() + ".L");
        aggregated[index] = new NormalizationColumn(getId() + ".Norm");

        return aggregated;
        //return lrm.getColumns();
    }

    public double getLogLikelihood() {
        double logL = super.getLogLikelihood();
        if (logL == 0 &&
                BayesianStochasticSearchVariableSelection.Utils.connectedAndWellConditioned(testProbabilities,this)) { // Also check that graph is connected
            return 0;
        }
        return Double.NEGATIVE_INFINITY;
    }

    @Override
    public String getDescription() { return "Log Rate Substitution Model"; }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CommonCitations.LEMEY_2014_UNIFYING); // TODO change citation
    }

    private final LogAdditiveCtmcRateProvider lrm;
    private final double[] testProbabilities;
}
