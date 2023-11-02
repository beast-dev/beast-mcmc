/*
 * RandomFieldParser.java
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

import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.distribution.OldGaussianMarkovRandomFieldModel;
import dr.inference.distribution.RandomField;
import dr.inference.model.Parameter;
import dr.math.distributions.RandomFieldDistribution;
import dr.xml.*;

public class RandomFieldParser extends AbstractXMLObjectParser {

    private static final String NORMAL_DISTRIBUTION_MODEL = "randomField";
    private static final String DATA = "data";
    private static final String DISTRIBUTION = "distribution";

    public String getParserName() {
        return NORMAL_DISTRIBUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter field = (Parameter) xo.getElementFirstChild(MultivariateDistributionLikelihood.DATA);

        RandomFieldDistribution distribution = (RandomFieldDistribution)
                xo.getElementFirstChild(DISTRIBUTION);

        if (field.getDimension() != distribution.getDimension()) {
            throw new XMLParseException("Field dimension (" + field.getDimension() +
                    ") != distribution dimension (" + distribution.getDimension() + ")");
        }

        String id = xo.hasId() ? xo.getId() : null;

        return new RandomField(id, field, distribution);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(DATA, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class),
            }),
            new ElementRule(DISTRIBUTION, new XMLSyntaxRule[] {
                    new ElementRule(RandomFieldDistribution.class),
            }),
    };

    public String getParserDescription() {
        return "Describes a normal distribution with a given mean and precision " +
                "that can be used in a distributionLikelihood element";
    }

    public Class getReturnType() {
        return OldGaussianMarkovRandomFieldModel.class;
    }

}
