/*
 * MultivariateNormalDistributionModelParser.java
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

package dr.inferencexml.distribution.shrinkage;

import dr.inference.distribution.shrinkage.*;
import dr.inference.model.Parameter;
import dr.inference.model.DuplicatedParameter;
import dr.inference.model.ParameterParser;
import dr.xml.*;

import static dr.inferencexml.distribution.shrinkage.BayesianBridgeLikelihoodParser.*;

public class BayesianBridgeDistributionModelParser extends AbstractXMLObjectParser {

    private static final String BAYESIAN_BRIDGE_DISTRIBUTION = "bayesianBridgeDistribution";
    private static final String DIMENSION = "dim";

    public String getParserName() {
        return BAYESIAN_BRIDGE_DISTRIBUTION;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject globalXo = xo.getChild(GLOBAL_SCALE);
        Parameter globalScale = (Parameter) globalXo.getChild(Parameter.class);

        Parameter localScale = null;
        int dim;
        if (xo.hasChildNamed(LOCAL_SCALE)) {
            XMLObject localXo = xo.getChild(LOCAL_SCALE);
            localScale = (Parameter) localXo.getChild(Parameter.class);
            if(localScale instanceof DuplicatedParameter) {
              throw new XMLParseException("Local scale cannot be a duplicated parameter");
            }
            dim = localScale.getDimension();

            if (xo.hasAttribute(DIMENSION) && (xo.getIntegerAttribute(DIMENSION) != dim)) {
                throw new XMLParseException("Invalid dimensions");
            }
        } else {
            dim = xo.getAttribute(DIMENSION, 1);
        }

        XMLObject exponentXo = xo.getChild(EXPONENT);
        Parameter exponent = (Parameter) exponentXo.getChild(Parameter.class);

        Parameter slabWidth = ParameterParser.getOptionalParameter(xo, SLAB_WIDTH);

        if (localScale == null && slabWidth != null) {
            throw new XMLParseException("Slab-regularization is only available under the joint Bayesian bridge");
        }

        if (localScale != null) {
            return new JointBayesianBridgeDistributionModel(globalScale, localScale, exponent, slabWidth, dim);
        } else {
            return new MarginalBayesianBridgeDistributionModel(globalScale, exponent, dim);
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(GLOBAL_SCALE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(EXPONENT,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(LOCAL_SCALE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(SLAB_WIDTH,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),               
            AttributeRule.newIntegerRule(DIMENSION, true),
    };

    public String getParserDescription() {
        return "Describes a scaled mixture of normals distribution with a given global and local scale " +
                "that can be used in a distributionLikelihood element";
    }

    public Class getReturnType() {
        return BayesianBridgeDistributionModel.class;
    }
}
