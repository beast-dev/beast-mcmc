/*
 * LogFileTraceExporterParser.java
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

package dr.inferencexml.trace;

import dr.inference.trace.LogFileTraceExporter;
import dr.util.FileHelpers;
import dr.xml.*;

import java.io.File;

/**
 *
 */
public class LogFileTraceExporterParser extends AbstractXMLObjectParser {
    public static final String LOG_FILE_TRACE = "logFileTrace";

    private static final String FILENAME = "fileName";
    private static final String BURN_IN = "burnIn";

    public String getParserName() {
        return LOG_FILE_TRACE;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final File file = FileHelpers.getFile(xo.getStringAttribute(FILENAME));
        final int burnIn = xo.getAttribute(BURN_IN, -1);

        try {
            return new LogFileTraceExporter(file, burnIn);
        } catch (Exception e) {
            throw new XMLParseException(e.getMessage());
        }
    }

    public String getParserDescription() {
        return "reconstruct population graph from variable dimension run.";
    }

    public Class getReturnType() {
        return LogFileTraceExporter.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(FILENAME, false, "trace log."),
            AttributeRule.newIntegerRule(BURN_IN, true,
                    "The number of states (not sampled states, but actual states) that are discarded from the" +
                            " beginning of the trace before doing the analysis"),
    };
}
