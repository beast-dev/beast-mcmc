/*
 * HamiltonianMonteCarloOperatorParser.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inferencexml.operators.hmc;

import dr.inference.hmc.PrecisionColumnProvider;
import dr.inference.hmc.ReversibleHMCProvider;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.hmc.NoUTurnOperator;
import dr.xml.*;

import static dr.evomodelxml.continuous.hmc.TaskPoolParser.THREAD_COUNT;

/**
 * @author Zhenyu Zhang
 */

public class NoUTurnOperatorParser extends AbstractXMLObjectParser {

    private final static String NUTS = "NoUTurnOperator";
    private final static String ADAPTIVE_STEPSIZE_FLG = "adaptiveStepsize";
    private final static String ADAPTIVE_STEPSIZE_DELAY = "adaptiveDelay";

    @Override
    public String getParserName() {
        return NUTS;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        ReversibleHMCProvider reversibleHMCprovider = (ReversibleHMCProvider) xo.getChild(ReversibleHMCProvider.class);
        boolean adaptiveStepsize = xo.getAttribute(ADAPTIVE_STEPSIZE_FLG, true);
        int adaptiveDelay = xo.getAttribute(ADAPTIVE_STEPSIZE_DELAY, -1);
        return new NoUTurnOperator(reversibleHMCprovider, adaptiveStepsize, adaptiveDelay, weight);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    final static XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT)
    };

    private final XMLSyntaxRule[] additionalRules = {
            new ElementRule(PrecisionColumnProvider.class),
            AttributeRule.newIntegerRule(THREAD_COUNT, true),
    };

    @Override
    public String getParserDescription() {
        return "Returns a NUTS transition kernel";
    }

    @Override
    public Class getReturnType() {
        return NoUTurnOperator.class;
    }
}
