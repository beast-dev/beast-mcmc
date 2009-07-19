/*
 * GTR.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.beagle.substmodel;

import dr.evolution.datatype.Nucleotides;
import dr.inference.model.Parameter;

/**
 * General Time Reversible model of nucleotide evolution
 * This is really just a place-holder because all the implementation
 * already exists in NucleotideModel and GeneralModel, its base classes.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: GTR.java,v 1.19 2005/05/24 20:25:58 rambaut Exp $
 */
public class GTR extends BaseSubstitutionModel {

    private Parameter rateACParameter = null;
    private Parameter rateAGParameter = null;
    private Parameter rateATParameter = null;
    private Parameter rateCGParameter = null;
    private Parameter rateCTParameter = null;
    private Parameter rateGTParameter = null;

    /**
     * @param rateACParameter rate of A<->C substitutions
     * @param rateAGParameter rate of A<->G substitutions
     * @param rateATParameter rate of A<->T substitutions
     * @param rateCGParameter rate of C<->G substitutions
     * @param rateCTParameter rate of C<->T substitutions
     * @param rateGTParameter rate of G<->T substitutions
     * @param freqModel       frequencies
     */
    public GTR(
            Parameter rateACParameter,
            Parameter rateAGParameter,
            Parameter rateATParameter,
            Parameter rateCGParameter,
            Parameter rateCTParameter,
            Parameter rateGTParameter,
            FrequencyModel freqModel) {

        super("GTR", Nucleotides.INSTANCE, freqModel);

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

    public void setAbsoluteRates(double[] rates, int relativeTo) {
        for (int i = 0; i < relativeRates.length; i++) {
            relativeRates[i] = rates[i] / rates[relativeTo];
        }
        updateMatrix = true;
        fireModelChanged();
    }

    public void setRelativeRates(double[] rates) {
        System.arraycopy(rates, 0, relativeRates, 0, relativeRates.length);
        updateMatrix = true;
        fireModelChanged();
    }

    protected void frequenciesChanged() {
        // nothing to do...
    }

    protected void ratesChanged() {
        // nothing to do...
    }

    protected void setupRelativeRates(double[] rates) {
        if (rateACParameter != null) {
            rates[0] = rateACParameter.getParameterValue(0);
        }
        if (rateAGParameter != null) {
            rates[1] = rateAGParameter.getParameterValue(0);
        }
        if (rateATParameter != null) {
            rates[2] = rateATParameter.getParameterValue(0);
        }
        if (rateCGParameter != null) {
            rates[3] = rateCGParameter.getParameterValue(0);
        }
        if (rateCTParameter != null) {
            rates[4] = rateCTParameter.getParameterValue(0);
        }
        if (rateGTParameter != null) {
            rates[5] = rateGTParameter.getParameterValue(0);
        }
    }

    // **************************************************************
    // XHTMLable IMPLEMENTATION
    // **************************************************************

    public String toXHTML() {
        StringBuffer buffer = new StringBuffer();

        buffer.append("<em>GTR Model</em> Instantaneous Rate Matrix = <table><tr><td></td><td>A</td><td>C</td><td>G</td><td>T</td></tr>");
        buffer.append("<tr><td>A</td><td></td><td>");
        buffer.append(relativeRates[0]);
        buffer.append("</td><td>");
        buffer.append(relativeRates[1]);
        buffer.append("</td><td>");
        buffer.append(relativeRates[2]);
        buffer.append("</td></tr>");

        buffer.append("<tr><td>C</td><td></td><td></td><td>");
        buffer.append(relativeRates[3]);
        buffer.append("</td><td>");
        buffer.append(relativeRates[4]);
        buffer.append("</td></tr>");

        buffer.append("<tr><td>G</td><td></td><td></td><td></td><td>");
        buffer.append(relativeRates[5]);
        buffer.append("</td></tr>");

        buffer.append("<tr><td>G</td><td></td><td></td><td></td><td></td></tr></table>");

        return buffer.toString();
    }

}