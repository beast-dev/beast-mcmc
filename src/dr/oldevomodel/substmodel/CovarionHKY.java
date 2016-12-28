/*
 * CovarionHKY.java
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

import dr.evolution.datatype.OldHiddenNucleotides;
import dr.oldevomodelxml.substmodel.CovarionHKYParser;
import dr.inference.model.Parameter;

/**
 * A model with hidden states that represent different rates.
 *
 * @author Alexei Drummond
 * @version $Id: CovarionHKY.java,v 1.4 2005/05/24 20:25:58 rambaut Exp $
 */
public class CovarionHKY extends AbstractCovarionDNAModel {
    
    /**
     * kappa
     */
    private Parameter kappaParameter;

    /**
     * @param dataType         the datatype to be used
     * @param kappaParameter   the rate of transitions versus transversions
     * @param hiddenClassRates the relative rates of the hidden categories
     *                         (first hidden category has rate 1.0 so this parameter
     *                         has dimension one less than number of hidden categories.
     *                         each hidden category.
     * @param switchingRates   rate of switching between hidden categories
     * @param freqModel        the frequencies
     */
    public CovarionHKY(OldHiddenNucleotides dataType, Parameter kappaParameter, Parameter hiddenClassRates, Parameter switchingRates, FrequencyModel freqModel) {

        super(CovarionHKYParser.COVARION_HKY, dataType, hiddenClassRates, switchingRates, freqModel);

        this.kappaParameter = kappaParameter;
        addVariable(kappaParameter);
        setupRelativeRates();
    }

    double[] getRelativeDNARates() {
        double kappa = kappaParameter.getParameterValue(0);
        return new double[]{1.0, kappa, 1.0, 1.0, kappa, 1.0};
    }

    /**
     * set kappa
     *
     * @param kappa the new value of kappa to set
     */
    public void setKappa(double kappa) {
        kappaParameter.setParameterValue(0, kappa);
    }

    /**
     * @return kappa
     */
    public double getKappa() {
        return kappaParameter.getParameterValue(0);
    }

    public String toString() {
        StringBuilder builder = new StringBuilder();
        builder.append("Covarion HKY model with ");
        builder.append(getHiddenClassCount()).append(" rate classes.\n");
        builder.append("Relative rates: \n");
        builder.append(SubstitutionModelUtils.toString(relativeRates, dataType, true, 2));
        return builder.toString();

    }

}