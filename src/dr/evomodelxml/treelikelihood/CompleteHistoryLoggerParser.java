/*
 * CompleteHistoryLoggerParser.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.treelikelihood;

import dr.evomodel.treelikelihood.MarkovJumpsTraitProvider;
import dr.evomodel.treelikelihood.utilities.CompleteHistoryLogger;
import dr.evomodel.treelikelihood.utilities.HistoryFilter;
import dr.inference.loggers.Logger;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */
public class CompleteHistoryLoggerParser extends AbstractXMLObjectParser {

    public static final String NAME = "completeHistoryLogger";
    public static final String EXTERNAL = "external";
    public static final String INTERNAL = "internal";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        boolean logExternal = xo.getAttribute(EXTERNAL, true);
        boolean logInternal = xo.getAttribute(INTERNAL, true);
        MarkovJumpsTraitProvider treeLikelihood =
                (MarkovJumpsTraitProvider) xo.getChild(MarkovJumpsTraitProvider.class);

        HistoryFilter filter = (HistoryFilter) xo.getChild(HistoryFilter.class);

        return new CompleteHistoryLogger(treeLikelihood, filter, logInternal, logExternal);
    }

    public String getParserName() {
        return NAME;
    }

    public String getParserDescription() {
        return "A logger to record all transitions in the complete history.";
    }

    public Class getReturnType() {
        return Logger.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(EXTERNAL, true),
            AttributeRule.newBooleanRule(INTERNAL, true),
            new ElementRule(MarkovJumpsTraitProvider.class),
            new ElementRule(HistoryFilter.class, true),
    };
}
