/*
 * PosteriorPredictiveSequenceLoggerParser.java
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

package dr.evomodelxml.treedatalikelihood.pps;

import dr.evomodel.treedatalikelihood.pps.PosteriorPredictiveSequenceLogger;
import dr.evomodel.treedatalikelihood.pps.PredictiveDataGenerator;
import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inferencexml.loggers.LoggerParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.StringAttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

import java.io.PrintWriter;

/**
 * Parses a logPosteriorPredictive element.
 */
public class PosteriorPredictiveSequenceLoggerParser extends LoggerParser {

    public static final String LOG_POSTERIOR_PREDICTIVE = "logPosteriorPredictive";

    public String getParserName() {
        return LOG_POSTERIOR_PREDICTIVE;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final long logEvery = xo.getLongIntegerAttribute(LOG_EVERY);

        final PrintWriter pw = getLogFile(xo, getParserName());
        final LogFormatter formatter = new TabDelimitedFormatter(pw);

        final PredictiveDataGenerator generator = (PredictiveDataGenerator) xo.getChild(PredictiveDataGenerator.class);

        return new PosteriorPredictiveSequenceLogger(generator, formatter, logEvery);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newLongIntegerRule(LOG_EVERY),
            AttributeRule.newBooleanRule(ALLOW_OVERWRITE_LOG, true),
            new StringAttributeRule(FILE_NAME,
                    "The name of the file to log posterior predictive datasets to", true),
            new ElementRule(PredictiveDataGenerator.class)
    };

    public String getParserDescription() {
        return "Logs posterior predictive sequence datasets, simulated from the current model state, " +
                "as a single growing native-BEAST-XML file.";
    }

    public Class getReturnType() {
        return PosteriorPredictiveSequenceLogger.class;
    }
}
