/*
 * TN93Parser.java
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

package dr.oldevomodelxml.substmodel;

import dr.oldevomodel.substmodel.FrequencyModel;
import dr.oldevomodel.substmodel.NucModelType;
import dr.oldevomodel.substmodel.TN93;
import dr.inference.model.Variable;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * Parses an element from an DOM document into a DemographicModel. Recognises
 * ConstantPopulation and ExponentialGrowth.
 */
public class TN93Parser extends AbstractXMLObjectParser {
    public static final String TN93_MODEL = NucModelType.TN93.getXMLName();
    public static final String KAPPA1 = "kappa1";
    public static final String KAPPA2 = "kappa2";
    public static final String FREQUENCIES = "frequencies";

    public String getParserName() {
        return TN93_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Variable kappa1Param = (Variable) xo.getElementFirstChild(KAPPA1);
        Variable kappa2Param = (Variable) xo.getElementFirstChild(KAPPA2);
        FrequencyModel freqModel = (FrequencyModel) xo.getElementFirstChild(FREQUENCIES);

        Logger.getLogger("dr.evomodel").info("Creating TN93 substitution model. Initial kappa = "
                + kappa1Param.getValue(0) + "," + kappa2Param.getValue(0));

        return new TN93(kappa1Param, kappa2Param, freqModel);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an instance of the TN93 (Tamura and Nei 1993) model of nucleotide evolution.";
    }

    public Class getReturnType() {
        return TN93.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules;{
        rules = new XMLSyntaxRule[]{
                new ElementRule(FREQUENCIES,
                        new XMLSyntaxRule[]{new ElementRule(FrequencyModel.class)}),
                new ElementRule(KAPPA1,
                        new XMLSyntaxRule[]{new ElementRule(Variable.class)}),
                new ElementRule(KAPPA2,
                        new XMLSyntaxRule[]{new ElementRule(Variable.class)})
        };
    }


}
