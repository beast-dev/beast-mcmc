/*
 * GTR.java
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

package dr.evomodel.substmodel.nucleotide;

import dr.evomodel.substmodel.BaseSubstitutionModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evolution.datatype.Nucleotides;
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
public class GTR extends BaseSubstitutionModel implements Citable {

    public static final String A_TO_C = "AC";
    public static final String A_TO_G = "AG";
    public static final String A_TO_T = "AT";
    public static final String C_TO_G = "CG";
    public static final String C_TO_T = "CT";
    public static final String G_TO_T = "GT";

    private Parameter ratesParameter = null;

    private Variable<Double> rateACVariable = null;
    private Variable<Double> rateAGVariable = null;
    private Variable<Double> rateATVariable = null;
    private Variable<Double> rateCGVariable = null;
    private Variable<Double> rateCTVariable = null;
    private Variable<Double> rateGTVariable = null;

    /**
     * @param rates           vector of all 6 rates
     * @param freqModel       frequencies
     */
    public GTR(
            Parameter rates,
            FrequencyModel freqModel) {
        super("GTR", Nucleotides.INSTANCE, freqModel);
        this.ratesParameter = rates;
        this.ratesParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, rates.getDimension()));
        addVariable(rates);
    }

    /**
     * @param rateACVariable rate of A<->C substitutions
     * @param rateAGVariable rate of A<->G substitutions
     * @param rateATVariable rate of A<->T substitutions
     * @param rateCGVariable rate of C<->G substitutions
     * @param rateCTVariable rate of C<->T substitutions
     * @param rateGTVariable rate of G<->T substitutions
     * @param freqModel       frequencies
     */
    public GTR(
            Variable<Double> rateACVariable,
            Variable<Double> rateAGVariable,
            Variable<Double> rateATVariable,
            Variable<Double> rateCGVariable,
            Variable<Double> rateCTVariable,
            Variable<Double> rateGTVariable,
            FrequencyModel freqModel) {

        super("GTR", Nucleotides.INSTANCE, freqModel);

        if (rateACVariable != null) {
            addVariable(rateACVariable);
            rateACVariable.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateACVariable = rateACVariable;
        }

        if (rateAGVariable != null) {
            addVariable(rateAGVariable);
            rateAGVariable.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateAGVariable = rateAGVariable;
        }

        if (rateATVariable != null) {
            addVariable(rateATVariable);
            rateATVariable.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateATVariable = rateATVariable;
        }

        if (rateCGVariable != null) {
            addVariable(rateCGVariable);
            rateCGVariable.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateCGVariable = rateCGVariable;
        }

        if (rateCTVariable != null) {
            addVariable(rateCTVariable);
            rateCTVariable.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateCTVariable = rateCTVariable;
        }

        if (rateGTVariable != null) {
            addVariable(rateGTVariable);
            rateGTVariable.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateGTVariable = rateGTVariable;
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
        if (ratesParameter == null) {
            if (rateACVariable != null) {
                rates[0] = rateACVariable.getValue(0);
            }
            if (rateAGVariable != null) {
                rates[1] = rateAGVariable.getValue(0);
            }
            if (rateATVariable != null) {
                rates[2] = rateATVariable.getValue(0);
            }
            if (rateCGVariable != null) {
                rates[3] = rateCGVariable.getValue(0);
            }
            if (rateCTVariable != null) {
                rates[4] = rateCTVariable.getValue(0);
            }
            if (rateGTVariable != null) {
                rates[5] = rateGTVariable.getValue(0);
            }
        } else {
            System.arraycopy(this.ratesParameter.getParameterValues(), 0, rates, 0, rates.length);
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
            "In: Miura R. M., editor. Lectures on Mathematics in the Life Sciences",
            17, 57, 86
    );

}