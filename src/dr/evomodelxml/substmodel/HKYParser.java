/*
 * HKYParser.java
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

package dr.evomodelxml.substmodel;

import dr.inference.model.Parameter;
import dr.xml.*;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.HKY;

import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class HKYParser extends AbstractXMLObjectParser {

    public static final String HKY_MODEL = "hkyModel";
    public static final String KAPPA = "kappa";
    public static final String FREQUENCIES = "frequencies";

    public String getParserName() {
        return HKY_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter kappaParam = (Parameter) xo.getElementFirstChild(KAPPA);
        FrequencyModel freqModel = (FrequencyModel) xo.getElementFirstChild(FrequencyModelParser.FREQUENCIES);

        if (kappaParam.getId() == null && kappaParam.getParameterValue(0) == 1.0) {
            Logger.getLogger("dr.evomodel").info("\nCreating JC substitution model.");
        } else {
            Logger.getLogger("dr.evomodel").info("\nCreating HKY substitution model. Initial kappa = " +
                    kappaParam.getParameterValue(0));
        }
        return new HKY(kappaParam, freqModel);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an instance of the HKY85 " +
                "(Hasegawa, Kishino & Yano, 1985) model of nucleotide evolution.";
    }

    public Class getReturnType() {
        return HKY.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(FrequencyModelParser.FREQUENCIES,
                    new XMLSyntaxRule[]{new ElementRule(FrequencyModel.class)}),
            new ElementRule(KAPPA,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
    };

}