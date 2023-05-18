/*
 * DirectionalBayesianBridgeMarkovRandomFieldLikelihoodParser.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.inference.distribution.BayesianBridgeMarkovRandomFieldLikelihood;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.distribution.shrinkage.*;
import dr.inference.model.Parameter;
import dr.inference.model.ParameterParser;
import dr.inferencexml.distribution.shrinkage.BayesianBridgeDistributionModelParser;
import dr.util.FirstOrderFiniteDifferenceTransform;
import dr.util.InverseFirstOrderFiniteDifferenceTransform;
import dr.util.Transform;
import dr.xml.*;

public class BayesianBridgeMarkovRandomFieldLikelihoodParser extends AbstractXMLObjectParser {

    public static final String NAME = "bayesianBridgeMarkovRandomField";
    public static final String BAYESIAN_BRIDGE_DISTRIBUTION = "incrementDistribution";
    public static final String FIRST_ELEMENT_DISTRIBUTION = "firstElementDistribution";
    public static final String INCREMENT_TRANSFORM = "type";

    public String getParserName() {
        return NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter variables = (Parameter) xo.getChild(Parameter.class);

//        BayesianBridgeDistributionModelParser bridgeParser = new BayesianBridgeDistributionModelParser();
//        BayesianBridgeDistributionModel bridge = (BayesianBridgeDistributionModel) bridgeParser.parseXMLObject(xo.getChild(BAYESIAN_BRIDGE_DISTRIBUTION).getChild(0));
        BayesianBridgeDistributionModel bridge = (BayesianBridgeDistributionModel) xo.getChild(BAYESIAN_BRIDGE_DISTRIBUTION).getChild(0);

        ParametricDistributionModel firstElementDistribution = (ParametricDistributionModel) xo.getChild(FIRST_ELEMENT_DISTRIBUTION).getChild(0);

        double upper = xo.getAttribute("upper", 1.0);
        double lower = xo.getAttribute("lower", 0.0);

        String ttype = (String) xo.getAttribute(INCREMENT_TRANSFORM, "none");
        Transform.UnivariableTransform incrementTransform;
        if ( ttype.equalsIgnoreCase("none") ) {
            incrementTransform = new Transform.NoTransform();
        } else if ( ttype.equalsIgnoreCase("log") ) {
            incrementTransform = new Transform.LogTransform();
        } else if ( ttype.equalsIgnoreCase("logit") ) {
            incrementTransform = new Transform.ScaledLogitTransform(lower, upper);
        } else {
            throw new RuntimeException("Invalid option for "+ INCREMENT_TRANSFORM);
        }

        int dim = variables.getDimension();

        InverseFirstOrderFiniteDifferenceTransform transform = new InverseFirstOrderFiniteDifferenceTransform(dim, incrementTransform);

        return new BayesianBridgeMarkovRandomFieldLikelihood(variables, bridge, firstElementDistribution, transform);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Parameter.class),
            new StringAttributeRule(INCREMENT_TRANSFORM, "The scale on which this model lives."),
            new ElementRule(BAYESIAN_BRIDGE_DISTRIBUTION,
                    new XMLSyntaxRule[]{new ElementRule(BayesianBridgeDistributionModel.class)}),
            new ElementRule(FIRST_ELEMENT_DISTRIBUTION,
                    new XMLSyntaxRule[]{new ElementRule(ParametricDistributionModel.class)})
    };

    public String getParserDescription() {
        return "Describes a scaled mixture of normals distribution with a given global and local scale " +
                "that can be used in a distributionLikelihood element";
    }

    public Class getReturnType() {
        return BayesianBridgeLikelihood.class;
    }
}
