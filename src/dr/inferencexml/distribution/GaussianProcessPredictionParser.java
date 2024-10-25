/*
 * GaussianProcessPredictionParser.java
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

import dr.inference.distribution.RandomField;
import dr.inference.model.DesignMatrix;
import dr.inference.model.Parameter;
import dr.math.distributions.gp.AdditiveGaussianProcessDistribution;
import dr.math.distributions.gp.GaussianProcessPrediction;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

public class GaussianProcessPredictionParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "gaussianProcessPrediction";
    private static final String BASES = "bases";

    public String getParserName() { return PARSER_NAME; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        AdditiveGaussianProcessDistribution gp = (AdditiveGaussianProcessDistribution)
                xo.getChild(AdditiveGaussianProcessDistribution.class);

        Parameter realizedValues = (Parameter) xo.getChild(Parameter.class);
        List<DesignMatrix> predictiveDesigns = new ArrayList<>(xo.getChild(BASES).getAllChildren(DesignMatrix.class));
        //List<DesignMatrix> predictiveDesigns = new ArrayList<>(xo.getAllChildren(DesignMatrix.class));

        return new GaussianProcessPrediction(gp, realizedValues, predictiveDesigns);
    }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(AdditiveGaussianProcessDistribution.class),
            new ElementRule(Parameter.class),
            new ElementRule(BASES, new XMLSyntaxRule[] {
                    new ElementRule(DesignMatrix.class, 1, Integer.MAX_VALUE)
            }),
    };

    public String getParserDescription() { // TODO update
        return "Describes a normal distribution with a given mean and precision " +
                "that can be used in a distributionLikelihood element";
    }

    public Class getReturnType() { return RandomField.class; }
}
