/*
 * EmpiricalAminoAcidModel.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.substmodel.aminoacid;

import dr.evomodel.substmodel.BaseSubstitutionModel;
import dr.evomodel.substmodel.EmpiricalRateMatrix;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evolution.datatype.AminoAcids;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class EmpiricalAminoAcidModel extends BaseSubstitutionModel {

    public EmpiricalAminoAcidModel(EmpiricalRateMatrix rateMatrix, FrequencyModel freqModel) {

        super(rateMatrix.getName(), rateMatrix.getDataType(), freqModel, null);

        if (freqModel == null) {
            double[] freqs = rateMatrix.getEmpiricalFrequencies();
            this.freqModel = new FrequencyModel(AminoAcids.INSTANCE, new Parameter.Default(freqs));
        }

        this.rateMatrix = rateMatrix;
    }

    protected void frequenciesChanged() {
        // Nothing to precalculate
    }

    protected void ratesChanged() {
        // Nothing to precalculate
    }

    protected void setupRelativeRates(double[] rates) {
        double[] empiricalRates = rateMatrix.getEmpiricalRates();
        System.arraycopy(empiricalRates, 0, rates, 0, empiricalRates.length);

    }

    private EmpiricalRateMatrix rateMatrix;
}
