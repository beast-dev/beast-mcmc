/*
 * HiddenLinkageLoggerParser.java
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

package dr.evomodelxml.tree;

import java.io.PrintWriter;

import dr.evomodel.tree.HiddenLinkageLogger;
import dr.evomodel.tree.HiddenLinkageModel;
import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inferencexml.loggers.LoggerParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.StringAttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Aaron Darling
 */
public class HiddenLinkageLoggerParser extends LoggerParser {
    public static final String LOG_HIDDEN_LINKAGE = "logHiddenLinkage";

    public String getParserName() {
        return LOG_HIDDEN_LINKAGE;
    }

    /**
     * @return an object based on the XML element it was passed.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
    	HiddenLinkageModel hlm = (HiddenLinkageModel)xo.getChild(HiddenLinkageModel.class);
        // logEvery of zero only displays at the end
        final int logEvery = xo.getAttribute(LOG_EVERY, 0);
        final PrintWriter pw = getLogFile(xo, getParserName());
        final LogFormatter formatter = new TabDelimitedFormatter(pw);

        return new HiddenLinkageLogger(hlm, formatter, logEvery);
    }

    public String getParserDescription() {
        return "Logs a linkage groups for metagenomic reads to a file";
    }

    public Class getReturnType() {
        return HiddenLinkageLogger.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(LOG_EVERY, true),
            new StringAttributeRule(FILE_NAME,
                    "The name of the file to send log output to. " +
                            "If no file name is specified then log is sent to standard output", true),
            new ElementRule(HiddenLinkageModel.class, "The linkage model which is to be logged"),
    };
}
