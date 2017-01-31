/*
 * CovarionGTR.java
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
import dr.oldevomodelxml.substmodel.CovarionGTRParser;
import dr.inference.model.Parameter;

/**
 * A general time reversible model of nucleotide evolution with
 * covarion hidden rate categories.
 *
 * @author Alexei Drummond
 * @version $Id$
 */
public class CovarionGTR extends AbstractCovarionDNAModel {
    /**
     * constructor
     *
     * @param dataType        the data type
     * @param freqModel       the equlibrium frequencies
     * @param hiddenRates     the relative rates of the hidden categories
     *                        (first hidden category has rate 1.0 so this parameter
     *                        has dimension one less than number of hidden categories.
     * @param rateACParameter the relative rate of A<->C substitutions
     * @param rateAGParameter the relative rate of A<->G substitutions
     * @param rateATParameter the relative rate of A<->T substitutions
     * @param rateCGParameter the relative rate of C<->G substitutions
     * @param rateCTParameter the relative rate of C<->T substitutions
     * @param rateGTParameter the relative rate of G<->T substitutions
     *                        each hidden category.
     * @param switchingRates  rates of switching between hidden rate classes
     */
    public CovarionGTR(
            OldHiddenNucleotides dataType,
            Parameter hiddenRates,
            Parameter switchingRates,
            Parameter rateACParameter,
            Parameter rateAGParameter,
            Parameter rateATParameter,
            Parameter rateCGParameter,
            Parameter rateCTParameter,
            Parameter rateGTParameter,
            FrequencyModel freqModel) {

        super(CovarionGTRParser.GTR_COVARION_MODEL, dataType, hiddenRates, switchingRates, freqModel);
        if (rateACParameter != null) {
            addVariable(rateACParameter);
            rateACParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateACParameter = rateACParameter;
        }

        if (rateAGParameter != null) {
            addVariable(rateAGParameter);
            rateAGParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateAGParameter = rateAGParameter;
        }

        if (rateATParameter != null) {
            addVariable(rateATParameter);
            rateATParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateATParameter = rateATParameter;
        }

        if (rateCGParameter != null) {
            addVariable(rateCGParameter);
            rateCGParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateCGParameter = rateCGParameter;
        }

        if (rateCTParameter != null) {
            addVariable(rateCTParameter);
            rateCTParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateCTParameter = rateCTParameter;
        }

        if (rateGTParameter != null) {
            addVariable(rateGTParameter);
            rateGTParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateGTParameter = rateGTParameter;
        }
    }

    double[] getRelativeDNARates() {
        double[] rr = new double[6];
        for (int i = 0; i < rr.length; i++) {
            rr[i] = 1.0;
        }

        if (rateACParameter != null) {
            rr[0] = rateACParameter.getParameterValue(0);
        }
        if (rateAGParameter != null) {
            rr[1] = rateAGParameter.getParameterValue(0);
        }
        if (rateATParameter != null) {
            rr[2] = rateATParameter.getParameterValue(0);
        }
        if (rateCGParameter != null) {
            rr[3] = rateCGParameter.getParameterValue(0);
        }
        if (rateCTParameter != null) {
            rr[4] = rateCTParameter.getParameterValue(0);
        }
        if (rateGTParameter != null) {
            rr[5] = rateGTParameter.getParameterValue(0);
        }
        return rr;
    }

    private Parameter rateACParameter = null;
    private Parameter rateAGParameter = null;
    private Parameter rateATParameter = null;
    private Parameter rateCGParameter = null;
    private Parameter rateCTParameter = null;
    private Parameter rateGTParameter = null;
}
