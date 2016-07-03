/*
 * EmpiricalAminoAcidModel.java
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

import dr.evolution.datatype.AminoAcids;
import dr.inference.model.Parameter;

/**
 * <b>A general model of sequence substitution</b>. A general reversible class for any
 * data type.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: EmpiricalAminoAcidModel.java,v 1.15 2005/05/24 20:25:58 rambaut Exp $
 */
public class EmpiricalAminoAcidModel extends AbstractAminoAcidModel {
    /**
     * constructor
     *
     * @param rateMatrix      EmpiricalRateMatrix
     * @param freqModel       the frequency model
     */
    public EmpiricalAminoAcidModel(EmpiricalRateMatrix rateMatrix, FrequencyModel freqModel) {

        super(rateMatrix.getName(), freqModel);

        if (freqModel == null) {
            areFrequenciesConstant = true;

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

    protected void setupRelativeRates() {
        double[] rates = rateMatrix.getEmpiricalRates();
        System.arraycopy(rates, 0, relativeRates, 0, relativeRates.length);
    }

    // *****************************************************************
    // Interface Model
    // *****************************************************************


    protected void storeState() {
    } // nothing to do

    protected void restoreState() {
        updateMatrix = !areFrequenciesConstant;
    }

    protected void acceptState() {
    } // nothing to do

    private EmpiricalRateMatrix rateMatrix;
    private boolean areFrequenciesConstant = false;
}
