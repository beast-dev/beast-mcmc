/*
 * RandomWalkModelParser.java
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
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.distribution.RandomWalkModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 *
 */
public class RandomWalkModelParser extends AbstractXMLObjectParser {

    public static final String RANDOM_WALK = "randomWalk";
    public static final String LOG_SCALE = "logScale";

    public String getParserName() {
        return RANDOM_WALK;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter data = (Parameter) xo.getChild(Parameter.class);
        ParametricDistributionModel distribution = (ParametricDistributionModel) xo.getChild(ParametricDistributionModel.class);

        boolean logScale = false;
        if (xo.hasAttribute(LOG_SCALE))
            logScale = xo.getBooleanAttribute(LOG_SCALE);

        return new RandomWalkModel(distribution, data, false, logScale);

    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newBooleanRule(LOG_SCALE, true),
            new ElementRule(Parameter.class),
            new XORRule(
                    new ElementRule(ParametricDistributionModel.class),
                    new ElementRule(DistributionLikelihood.class)
            )
    };

    public String getParserDescription() {
        return "Describes a first-order random walk. No prior is assumed on the first data element";
    }

    public Class getReturnType() {
        return RandomWalkModel.class;
    }

}
