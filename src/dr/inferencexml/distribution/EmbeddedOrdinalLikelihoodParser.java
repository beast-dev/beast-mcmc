/*
 * EmbeddedOrdinalLikelihoodParser.java
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

import dr.inference.distribution.EmbeddedOrdinalLikelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Filippo Monti (powered by OpenAI)
 */
public class EmbeddedOrdinalLikelihoodParser extends AbstractXMLObjectParser {

    public static final String LATENT_PARAMETER = "latent";
    public static final String LOG_WEIGHTS = "logWeights";
    public static final String CUTS = "cuts";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final Parameter latent = (Parameter) xo.getElementFirstChild(LATENT_PARAMETER);
        final Parameter logWeights = (Parameter) xo.getElementFirstChild(LOG_WEIGHTS);
        final Parameter cuts = (Parameter) xo.getElementFirstChild(CUTS);

        return new EmbeddedOrdinalLikelihood(latent, logWeights, cuts);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                new ElementRule(LATENT_PARAMETER, new XMLSyntaxRule[] {
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(LOG_WEIGHTS, new XMLSyntaxRule[] {
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(CUTS, new XMLSyntaxRule[] {
                        new ElementRule(Parameter.class)
                })
        };
    }

    @Override
    public String getParserDescription() {
        return "Likelihood for an embedded ordinal latent parameter with piecewise-constant log weights.";
    }

    @Override
    public Class getReturnType() {
        return EmbeddedOrdinalLikelihood.class;
    }

    @Override
    public String getParserName() {
        return EmbeddedOrdinalLikelihood.EMBEDDED_ORDINAL_LIKELIHOOD;
    }
}
