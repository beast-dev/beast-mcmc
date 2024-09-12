/*
 * MLOptimizerParser.java
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

package dr.inferencexml;

import dr.inference.loggers.Logger;
import dr.inference.ml.MLOptimizer;
import dr.inference.model.Likelihood;
import dr.inference.operators.OperatorSchedule;
import dr.xml.*;

import java.util.ArrayList;

/**
 *
 */
public class MLOptimizerParser extends AbstractXMLObjectParser {

    public static final String CHAIN_LENGTH = "chainLength";
    public static final String OPTIMIZER = "optimizer";    

    public String getParserName() { return OPTIMIZER; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        int chainLength = xo.getIntegerAttribute(CHAIN_LENGTH);

        OperatorSchedule opsched = null;
        dr.inference.model.Likelihood likelihood = null;
        ArrayList<Logger> loggers = new ArrayList<Logger>();

        for (int i = 0; i < xo.getChildCount(); i++) {
            Object child = xo.getChild(i);
            if (child instanceof dr.inference.model.Likelihood) {
                likelihood = (dr.inference.model.Likelihood)child;
            } else if (child instanceof OperatorSchedule) {
                opsched = (OperatorSchedule)child;
            } else if (child instanceof Logger) {
                loggers.add((Logger)child);
            } else {
                throw new XMLParseException("Unrecognized element found in optimizer element:" + child);
            }
        }

        Logger[] loggerArray = new Logger[loggers.size()];
        loggers.toArray(loggerArray);

        return new MLOptimizer("optimizer1", chainLength, likelihood, opsched, loggerArray);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************


    public String getParserDescription() {
        return "This element returns a maximum likelihood heuristic optimizer and runs the optimization as a side effect.";
    }

    public Class getReturnType() { return MLOptimizer.class; }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
        AttributeRule.newIntegerRule(CHAIN_LENGTH),
        new ElementRule(OperatorSchedule.class ),
        new ElementRule(Likelihood.class ),
        new ElementRule(Logger.class, 1, Integer.MAX_VALUE )
    };

}
