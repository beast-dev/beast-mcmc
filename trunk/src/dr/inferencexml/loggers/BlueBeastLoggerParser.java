/*
* BlueBeastLoggerParser.java
*
* Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
* BEAST is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with BEAST; if not, write to the
* Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
* Boston, MA  02110-1301  USA
*/

package dr.inferencexml.loggers;

//import beast.inference.loggers.BlueBeastLogger;
import beast.loggers.BlueBeastLogger;
import dr.app.beast.BeastVersion;
import dr.inference.loggers.*;
import dr.math.MathUtils;
import dr.util.Identifiable;
import dr.util.Property;
import dr.xml.*;

import java.io.PrintWriter;
import java.util.Date;

/**
* @author Wai Lok Li
*/
public class BlueBeastLoggerParser extends AbstractXMLObjectParser {

    public static final String LOG = "blueBeastLog";
    public static final String ECHO = "echo";
    public static final String ECHO_EVERY = "echoEvery";
    public static final String TITLE = "title";
//    public static final String FILE_NAME = FileHelpers.FILE_NAME;
    public static final String FORMAT = "format";
    public static final String TAB = "tab";
    public static final String HTML = "html";
    public static final String PRETTY = "pretty";
//    public static final String INITIAL_CHECK_INTERVAL = "initalLogInterval";
    public static final String LOG_EVERY = "logEvery";
//    public static final String ALLOW_OVERWRITE_LOG = "overwrite";

    public static final String COLUMNS = "columns";
    public static final String COLUMN = "column";
    public static final String LABEL = "label";
    public static final String SIGNIFICANT_FIGURES = "sf";
    public static final String DECIMAL_PLACES = "dp";
    public static final String WIDTH = "width";

    public String getParserName() {
        return LOG;
    }

    /**
     * @return an object based on the XML element it was passed.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        // You must say how often you want to log
        final int logEvery = xo.getIntegerAttribute(LOG_EVERY);
        // And also how often you want to check the logged values
//        int initialLogInterval = 1000;
//        if (xo.hasAttribute(INITIAL_CHECK_INTERVAL)) {
//            initialLogInterval = xo.getIntegerAttribute(INITIAL_CHECK_INTERVAL);
//            System.out.println("Starting log interval: " + initialLogInterval);
//        }

        final PrintWriter pw = getLogFile(xo, getParserName());

        final LogFormatter formatter = new TabDelimitedFormatter(pw);

        boolean performanceReport = false;

        // added a performance measurement delay to avoid the full evaluation period.
        //final BlueBeastLogger logger = new BlueBeastLogger(formatter, logEvery, initialLogInterval, performanceReport, 10000);

        final BlueBeastLogger logger = new BlueBeastLogger(formatter, logEvery, performanceReport, 10000);

        if (xo.hasAttribute(TITLE)) {
            logger.setTitle(xo.getStringAttribute(TITLE));
        } else {

            final BeastVersion version = new BeastVersion();

            final String title = "BEAST " + version.getVersionString() +
                    ", " + version.getBuildString() + "\n" +

                    "Generated " + (new Date()).toString() + " [seed=" + MathUtils.getSeed() + "]";
            logger.setTitle(title);
        }

        for (int i = 0; i < xo.getChildCount(); i++) {

            final Object child = xo.getChild(i);

//            if (child instanceof Columns) {
//
//                logger.addColumns(((Columns) child).getColumns());
//
//            } else
            if (child instanceof Loggable) {

                logger.add((Loggable) child);

            }
            else if (child instanceof Identifiable) {
                throw new RuntimeException("Inputted an Identifiable to the logger, not valid in BlueBeast");
//                logger.addColumn(new LogColumn.Default(((Identifiable) child).getId(), child));

            } else if (child instanceof Property) {
                throw new RuntimeException("Inputted a Property to the logger, not valid in BlueBeast");
//                logger.addColumn(new LogColumn.Default(((Property) child).getAttributeName(), child));
            } else {
                throw new RuntimeException("Inputted a wrong object type to the logger, not valid in BlueBeast");
//                logger.addColumn(new LogColumn.Default(child.getClass().toString(), child));
            }
        }

        return logger;
    }

    public static PrintWriter getLogFile(XMLObject xo, String parserName) throws XMLParseException {
        return XMLParser.getFilePrintWriter(xo, parserName);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
//            AttributeRule.newIntegerRule(INITIAL_CHECK_INTERVAL, true),
            AttributeRule.newIntegerRule(LOG_EVERY),
            new StringAttributeRule(TITLE,
                    "The title of the log", true),
            new OrRule(
                    new XMLSyntaxRule[]{
                            //new ElementRule(Columns.class, 1, Integer.MAX_VALUE),
                            new ElementRule(Loggable.class, 1, Integer.MAX_VALUE),
                            new ElementRule(Object.class, 1, Integer.MAX_VALUE)
                    }
            )
    };

    public String getParserDescription() {
        return "Logs one or more items at a given frequency to the screen or to a file";
    }

    public Class getReturnType() {
        return BlueBeastLogger.class;
    }
}

