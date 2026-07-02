/*
 * LikelihoodBasedDiscontinuousPotentialProviderParser.java
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

package dr.inferencexml.hmc;

import dr.inference.hmc.DiscontinuousPotentialProvider;
import dr.inference.hmc.LikelihoodBasedDiscontinuousPotentialProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Filippo Monti (powered by OpenAI)
 */
public class LikelihoodBasedDiscontinuousPotentialProviderParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "likelihoodBasedDiscontinuousPotential";
    public static final String QUIET = "quiet";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final Likelihood likelihood = (Likelihood) xo.getChild(Likelihood.class);
        final Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        final boolean quiet = xo.getAttribute(QUIET, false);

        if (parameter == null) {
            throw new XMLParseException("A parameter must be supplied to the likelihood-based discontinuous potential provider");
        }

        return new LikelihoodBasedDiscontinuousPotentialProvider(likelihood, parameter, quiet);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                AttributeRule.newBooleanRule(QUIET, true,
                        "If true, probe moves use setParameterValueQuietly and only dirty the wrapped likelihood. " +
                                "Default false fires normal parameter-change events and is production-correct but slower."),
                new ElementRule(Likelihood.class),
                new ElementRule(Parameter.class),
        };
    }

    @Override
    public String getParserDescription() {
        return "Wraps a likelihood and parameter as a discontinuous potential provider by reevaluating the likelihood after one-coordinate moves.";
    }

    @Override
    public Class getReturnType() {
        return DiscontinuousPotentialProvider.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }
}
