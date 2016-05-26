/*
 * DistributionLikelihoodParser.java
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

package dr.inferencexml.distribution;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Statistic;
import dr.math.distributions.RandomGenerator;
import dr.xml.*;

/**
 *
 */
public class DistributionLikelihoodParser extends AbstractXMLObjectParser {

    public static final String DISTRIBUTION = "distribution";
    public static final String DATA = "data";
    public static final String FROM = "from";
    public static final String TO = "to";


    public String getParserName() {
        return DistributionLikelihood.DISTRIBUTION_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final XMLObject cxo = xo.getChild(DISTRIBUTION);
        ParametricDistributionModel model = (ParametricDistributionModel) cxo.getChild(ParametricDistributionModel.class);

        DistributionLikelihood likelihood = new DistributionLikelihood(model);

        XMLObject cxo1 = xo.getChild(DATA);
        final int from = cxo1.getAttribute(FROM, -1);
        int to = cxo1.getAttribute(TO, -1);
        if (from >= 0 || to >= 0) {
            if (to < 0) {
                to = Integer.MAX_VALUE;
            }
            if (!(from >= 0 && to >= 0 && from < to)) {
                throw new XMLParseException("ill formed from-to");
            }
            likelihood.setRange(from, to);
        }

        for (int j = 0; j < cxo1.getChildCount(); j++) {
            if (cxo1.getChild(j) instanceof Statistic) {

                likelihood.addData((Statistic) cxo1.getChild(j));
            } else {
                throw new XMLParseException("illegal element in " + cxo1.getName() + " element");
            }
        }

        return likelihood;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(DISTRIBUTION,
                    new XMLSyntaxRule[]{new ElementRule(ParametricDistributionModel.class)}),
            new ElementRule(DATA, new XMLSyntaxRule[]{
                    AttributeRule.newIntegerRule(FROM, true),
                    AttributeRule.newIntegerRule(TO, true),
                    new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
            })
    };

    public String getParserDescription() {
        return "Calculates the likelihood of some data given some parametric or empirical distribution.";
    }

    public Class getReturnType() {
        return DistributionLikelihood.class;
    }
}