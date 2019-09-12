/*
 * DuplicationBranchSubstitutionModel.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.substmodel.geneconversion;

import dr.evomodel.substmodel.BaseSubstitutionModel;
import dr.inference.model.Parameter;

import java.util.List;

/**
 * @author Xiang Ji
 * @author Jeff Thorne
 * @author Marc A. Suchard
 */
public class DuplicationBranchSubstitutionModel extends BaseSubstitutionModel{

    private final List<PairedParalogGeneConversionSubstitutionModel> substitutionModels;
    private final Parameter timeUntilDuplicationProportion;

    public DuplicationBranchSubstitutionModel(String name,
                                              List<PairedParalogGeneConversionSubstitutionModel> substitutionModels,
                                              Parameter timeUntilDuplicationProportion) {
        super(name, substitutionModels.get(0).getDataType(), substitutionModels.get(0).getFrequencyModel());
        this.substitutionModels = substitutionModels;
        this.timeUntilDuplicationProportion = timeUntilDuplicationProportion;
    }

    @Override
    protected void frequenciesChanged() {

    }

    @Override
    protected void ratesChanged() {

    }

    @Override
    protected void setupRelativeRates(double[] rates) {

    }

    @Override
    protected double setupMatrix() {
        return 1.0;
    }

    @Override
    public void getTransitionProbabilities(double distance, double[] matrix) {
        final double singleParalogDistance = distance * timeUntilDuplicationProportion.getParameterValue(0);
        final double pairedParalogDistance = distance - singleParalogDistance;

        double[] singleParalogMatrix = new double[stateCount * stateCount];
        double[] pairedParalogMatrix = new double[stateCount * stateCount];
        substitutionModels.get(0).getTransitionProbabilities(singleParalogDistance, singleParalogMatrix);
        substitutionModels.get(1).getTransitionProbabilities(pairedParalogDistance, pairedParalogMatrix);

        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                double tmp = 0.0;
                for (int k = 0; k < stateCount; k++) {
                    tmp += singleParalogMatrix[i * stateCount + k] * pairedParalogMatrix[k * stateCount + j];
                }
                matrix[i * stateCount + j] = tmp;
            }
        }
    }

}
