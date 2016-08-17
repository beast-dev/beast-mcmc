/*
 * GTR.java
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

import dr.oldevomodelxml.substmodel.GTRParser;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Collections;
import java.util.List;

/**
 * General Time Reversible model of nucleotide evolution
 * This is really just a place-holder because all the implementation
 * already exists in NucleotideModel and GeneralModel, its base classes.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: GTR.java,v 1.19 2005/05/24 20:25:58 rambaut Exp $
 */
public class GTR extends AbstractNucleotideModel implements Citable {

    private Variable<Double> rateACValue = null;
    private Variable<Double> rateAGValue = null;
    private Variable<Double> rateATValue = null;
    private Variable<Double> rateCGValue = null;
    private Variable<Double> rateCTValue = null;
    private Variable<Double> rateGTValue = null;

    /**
     * @param rateACValue rate of A<->C substitutions
     * @param rateAGValue rate of A<->G substitutions
     * @param rateATValue rate of A<->T substitutions
     * @param rateCGValue rate of C<->G substitutions
     * @param rateCTValue rate of C<->T substitutions
     * @param rateGTValue rate of G<->T substitutions
     * @param freqModel       frequencies
     */
    public GTR(
            Variable rateACValue,
            Variable rateAGValue,
            Variable rateATValue,
            Variable rateCGValue,
            Variable rateCTValue,
            Variable rateGTValue,
            FrequencyModel freqModel) {

        super(GTRParser.GTR_MODEL, freqModel);

        if (rateACValue != null) {
            addVariable(rateACValue);
            rateACValue.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateACValue = rateACValue;
        }

        if (rateAGValue != null) {
            addVariable(rateAGValue);
            rateAGValue.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateAGValue = rateAGValue;
        }

        if (rateATValue != null) {
            addVariable(rateATValue);
            rateATValue.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateATValue = rateATValue;
        }

        if (rateCGValue != null) {
            addVariable(rateCGValue);
            rateCGValue.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateCGValue = rateCGValue;
        }

        if (rateCTValue != null) {
            addVariable(rateCTValue);
            rateCTValue.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateCTValue = rateCTValue;
        }

        if (rateGTValue != null) {
            addVariable(rateGTValue);
            rateGTValue.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateGTValue = rateGTValue;
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

    protected void setupRelativeRates() {

        if (rateACValue != null) {
            relativeRates[0] = rateACValue.getValue(0);
        }
        if (rateAGValue != null) {
            relativeRates[1] = rateAGValue.getValue(0);
        }
        if (rateATValue != null) {
            relativeRates[2] = rateATValue.getValue(0);
        }
        if (rateCGValue != null) {
            relativeRates[3] = rateCGValue.getValue(0);
        }
        if (rateCTValue != null) {
            relativeRates[4] = rateCTValue.getValue(0);
        }
        if (rateGTValue != null) {
            relativeRates[5] = rateGTValue.getValue(0);
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

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.SUBSTITUTION_MODELS;
    }

    @Override
    public String getDescription() {
        return "GTR nucleotide substitution model";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }


    public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("S", "Tavaré")
            },
            "Some probabilistic and statistical problems in the analysis of DNA sequences",
            1985,
            "In: Miura R. M., editor. Lectures on mathematics in the life sciences.",
            17, 57, 86
    );

}
