/*
 * TruncatedDistributionLikelihoodParser.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.inferencexml.distribution;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.Parameter;
import dr.inference.model.TruncatedDistributionLikelihood;
import dr.xml.*;

/**
 * @author Max Tolkoff
 */
public class TruncatedDistributionLikelihoodParser extends AbstractXMLObjectParser {
    public static final String TRUNCATED_DISTRIBUTION_LIKELIHOOD = "truncatedDistributionLikelihood";
    public static final String LOW = "low";
    public static final String HIGH = "high";

    @Override
    public String getParserName() {
        return TRUNCATED_DISTRIBUTION_LIKELIHOOD;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        DistributionLikelihood likelihood = (DistributionLikelihood) xo.getChild(DistributionLikelihood.class);
        Parameter low;
        Parameter high;
        if(xo.getChild(LOW) != null){
            low = (Parameter) xo.getChild(LOW).getChild(Parameter.class);
        }
        else{
            low = new Parameter.Default("low", Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
        }
        if(xo.getChild(HIGH) != null){
            high = (Parameter) xo.getChild(HIGH).getChild(Parameter.class);
        }
        else{
            high = new Parameter.Default("high", Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
        }

        return new TruncatedDistributionLikelihood(likelihood, low, high);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "Produces a truncated distribution likelihood";
    }

    @Override
    public Class getReturnType() {
        return TruncatedDistributionLikelihood.class;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(DistributionLikelihood.class),
            new ElementRule(LOW, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
            new ElementRule(HIGH, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
    };
}
