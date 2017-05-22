/*
 * IstvanOperatorParser.java
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

package dr.oldevomodelxml.indel;

import dr.oldevomodel.indel.IstvanOperator;
import dr.oldevomodel.indel.TKF91Likelihood;
import dr.xml.*;

/**
 *
 */
public class IstvanOperatorParser extends AbstractXMLObjectParser {

    public static final String ALIGMENT_CHUNK_OPERATOR = "alignmentChunkOperator";

    public String getParserName() {
        return ALIGMENT_CHUNK_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TKF91Likelihood likelihood = (TKF91Likelihood) xo.getChild(TKF91Likelihood.class);
        double weight = xo.getDoubleAttribute("weight");
        double iP = xo.getDoubleAttribute("iP");
        double exponent = xo.getDoubleAttribute("exponent");
        double gapPenalty = xo.getDoubleAttribute("gapPenalty");


        return new IstvanOperator(iP, exponent, gapPenalty, weight, likelihood);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an operator that re-aligns a small chunk of an alignment.";
    }

    public Class getReturnType() {
        return IstvanOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule("weight"),
            AttributeRule.newDoubleRule("iP", false, "tuning probability, values near zero resamples entire alignment and near 1.0 resamples single columns."),
            AttributeRule.newDoubleRule("exponent", false, "tuning parameter, value of 1.0 samples random alignments, large values (e.g. 2.7) sample alignment peaked around 'best' alignment"),
            AttributeRule.newDoubleRule("gapPenalty", false, "tuning parameter, must be negative, large values penalize gaps more in the alignment proposal."),
            new ElementRule(TKF91Likelihood.class)
    };
}
