/*
 * MLLoggerParser.java
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

import dr.inference.loggers.*;
import dr.inference.model.Likelihood;
import dr.util.Identifiable;
import dr.xml.*;

import java.io.PrintWriter;

/**
 *
 */
public class MLLoggerParser extends LoggerParser {

    public static final String LOG_ML = "logML";
    public static final String LIKELIHOOD = "ml";

    public String getParserName() {
        return LOG_ML;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Likelihood likelihood = (Likelihood) xo.getElementFirstChild(LIKELIHOOD);

        // logEvery of zero only displays at the end
        int logEvery = xo.getAttribute(LOG_EVERY, 0);

        final PrintWriter pw = getLogFile(xo, getParserName());

        LogFormatter formatter = new TabDelimitedFormatter(pw);

        MLLogger logger = new MLLogger(likelihood, formatter, logEvery);

        if (xo.hasAttribute(TITLE)) {
            logger.setTitle(xo.getStringAttribute(TITLE));
        }

        for (int i = 0; i < xo.getChildCount(); i++) {
            Object child = xo.getChild(i);

            if (child instanceof Columns) {

                logger.addColumns(((Columns) child).getColumns());

            } else if (child instanceof Loggable) {

                logger.add((Loggable) child);

            } else if (child instanceof Identifiable) {

                logger.addColumn(new LogColumn.Default(((Identifiable) child).getId(), child));

            } else {

                logger.addColumn(new LogColumn.Default(child.getClass().toString(), child));
            }
        }

        return logger;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(LOG_EVERY, true),
            new ElementRule(LIKELIHOOD,
                    new XMLSyntaxRule[]{new ElementRule(Likelihood.class)}),
            new OrRule(
                    new ElementRule(Columns.class, 1, Integer.MAX_VALUE),
                    new ElementRule(Loggable.class, 1, Integer.MAX_VALUE)
            )
    };

    public String getParserDescription() {
        return "Logs one or more items every time the given likelihood improves";
    }

    public Class getReturnType() {
        return MLLogger.class;
    }
}
