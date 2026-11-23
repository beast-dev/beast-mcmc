/*
 * ChainDiagnosticParser.java
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

package dr.inferencexml.loggers;

import dr.inference.loggers.ChainDiagnosticColumn;
import dr.xml.*;

/**
 * Parses the convergence diagnostics that can be put in the screen log. They are shown when
 * several independent chains are run (-chains) and are left out otherwise.
 *
 * @author Andrew Rambaut
 */
public class ChainDiagnosticParser extends AbstractXMLObjectParser {

    public static final String MIN_ESS = "minESS";
    public static final String MAX_RHAT = "maxRhat";
    public static final String ASDSF = "asdsf";
    public static final String COMMON_BURNIN = "commonBurnin";
    public static final String SHOW_NAME = "showName";

    public String getParserName() {
        return MIN_ESS;
    }

    @Override
    public String[] getParserNames() {
        return new String[]{MIN_ESS, MAX_RHAT, ASDSF, COMMON_BURNIN};
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        ChainDiagnosticColumn.Type type;
        String name = xo.getName();
        if (name.equals(MIN_ESS)) {
            type = ChainDiagnosticColumn.Type.MIN_ESS;
        } else if (name.equals(MAX_RHAT)) {
            type = ChainDiagnosticColumn.Type.MAX_RHAT;
        } else if (name.equals(ASDSF)) {
            type = ChainDiagnosticColumn.Type.ASDSF;
        } else {
            type = ChainDiagnosticColumn.Type.BURNIN;
        }
        boolean showName = xo.getAttribute(SHOW_NAME, false);
        if (showName && type != ChainDiagnosticColumn.Type.MIN_ESS && type != ChainDiagnosticColumn.Type.MAX_RHAT) {
            throw new XMLParseException("The " + SHOW_NAME + " attribute is only for " + MIN_ESS + " and " + MAX_RHAT);
        }
        return new ChainDiagnosticColumn(type, showName);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(SHOW_NAME, true,
                    "Follow the value with the name of the logged value that gives it (" + MIN_ESS + " and " + MAX_RHAT + ")"),
    };

    public String getParserDescription() {
        return "A convergence diagnostic of the independent chains of a run, for the screen log: " +
                MIN_ESS + " is the smallest effective sample size (summed over the chains) of the values in the main log, " +
                MAX_RHAT + " the largest split R-hat of those values, " +
                ASDSF + " the average standard deviation of clade frequencies and " +
                COMMON_BURNIN + " the burn-in (in states) after which the chains' posteriors agree. " +
                "All are computed after that burn-in.";
    }

    public Class getReturnType() {
        return ChainDiagnosticColumn.class;
    }
}
