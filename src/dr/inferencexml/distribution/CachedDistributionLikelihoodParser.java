/*
 * CachedDistributionLikelihoodParser.java
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

import dr.inference.distribution.AbstractDistributionLikelihood;
import dr.inference.distribution.CachedDistributionLikelihood;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.Variable;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Marc Suchard
 */
public class CachedDistributionLikelihoodParser extends AbstractXMLObjectParser {

    public static final String CACHED_PRIOR = "cachedPrior";

   // public static final String RATE_BLOCK = "rates";
    //public static final String INDICATOR_BLOCK = "indicators";

    public String getParserName() {
        return CACHED_PRIOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final String name = xo.hasId() ? xo.getId() : CACHED_PRIOR;

        final AbstractDistributionLikelihood likelihood = (AbstractDistributionLikelihood) xo.getChild(AbstractDistributionLikelihood.class);
        final Variable variable = (Variable) xo.getChild(Variable.class);

        final Logger logger = Logger.getLogger("dr.inference");
        logger.info("Constructing a cache around likelihood '" + likelihood.getId() + "', signal = " + variable.getVariableName());

        return new CachedDistributionLikelihood(name, likelihood, variable);

    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(AbstractDistributionLikelihood.class),
            new ElementRule(Variable.class),
    };

    public String getParserDescription() {
        return "Calculates a cached likelihood of some data given some parametric or empirical distribution.";
    }

    public Class getReturnType() {
        return DistributionLikelihood.class;
    }
}
