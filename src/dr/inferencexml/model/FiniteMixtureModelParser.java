/*
 * FiniteMixtureParameterParser.java
 *
 * Copyright © 2002-2026 the BEAST Development Team
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

package dr.inferencexml.model;

import dr.inference.model.FiniteMixtureModel;
import dr.inference.model.FiniteMixtureParameter;
import dr.math.distributions.Distribution;
import dr.math.distributions.GeometricDistribution;
import dr.math.distributions.PoissonDistribution;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class FiniteMixtureModelParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "finiteMixtureModel";
    private static final String CONSTRAINT = "constraint";
    private static final String PRIOR = "prior";
    private static final String HYPER_PARAMETER = "hyperParameter";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        FiniteMixtureParameter parameter = (FiniteMixtureParameter) xo.getChild(FiniteMixtureParameter.class);
        Distribution prior = parsePrior(xo);
        FiniteMixtureModel.Constraint constraint = parseConstraint(xo);

        return new FiniteMixtureModel(xo.getId(),parameter, prior, constraint);
    }

    private FiniteMixtureModel.Constraint parseConstraint(XMLObject xo) throws XMLParseException {
        String label = xo.getStringAttribute(CONSTRAINT, FiniteMixtureModel.Constraint.NONE.getName());
        FiniteMixtureModel.Constraint constraint = FiniteMixtureModel.Constraint.parse(label);
        if (constraint == null) {
            throw new XMLParseException("Unknown constraint");
        }
        return constraint;
    }

    private Distribution parsePrior(XMLObject xo) throws XMLParseException {
        String label = xo.getStringAttribute(PRIOR);
        double hyper = xo.getDoubleAttribute(HYPER_PARAMETER);

        if (label.equalsIgnoreCase("geometric")) {
            if (hyper <= 0.0 || hyper > 1.0) {
                throw new XMLParseException("Invalid hyper-parameter value");
            }
            return new GeometricDistribution(hyper);
        } else if (label.equalsIgnoreCase("poisson")) {
            if (hyper <= 0.0) {
                throw new XMLParseException("Invalid hyper-parameter value");
            }
            return new PoissonDistribution(hyper);
        } else {
            throw new XMLParseException("Unknown prior type");
        }
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(FiniteMixtureParameter.class),
            AttributeRule.newStringRule(CONSTRAINT, true),
            AttributeRule.newStringRule(PRIOR),
            AttributeRule.newDoubleRule(HYPER_PARAMETER),
    };

    public String getParserDescription() {
        return "A finite mixture parameter.";
    }

    public Class getReturnType() {
        return FiniteMixtureParameter.class;
    }

    public String getParserName() {
        return PARSER_NAME;
    }
}
