/*
 * Report.java
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

package dr.xml;

import java.io.PrintWriter;
import java.util.ArrayList;

/**
 * A generates a text report from the elements within it.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Marc A. Suchard
 * @version $Id: Report.java,v 1.15 2005/05/24 20:26:01 rambaut Exp $
 */
public class Report {

    public static final String REPORT = "report";
    public static final String FILENAME = "fileName";

    protected String title = "";
    protected ArrayList<Object> objects = new ArrayList<Object>();
    private PrintWriter writer;

    public void createReport() {
    	
		if (!title.equalsIgnoreCase("")) {
			writer.println(getTitle());
			writer.println();
		}
 
        for (Object object : objects) {
            final String item;

            if (object instanceof Reportable) {
                item = ((Reportable) object).getReport();
            } else {
                item = object.toString();
            }
            writer.print(item.trim());
            writer.print("\n");
        }
        writer.println();
        writer.flush();
    }

    public void setOutput(PrintWriter writer) {
        this.writer = writer;
    }

    /**
     * @param title the title of the report.
     */
    public void setTitle(String title) {
        this.title = title;
    }

    /**
     * @return the report's title.
     */
    public String getTitle() {
        return title;
    }

    /**
     * @param object the object to be added
     */
    public void add(Object object) {
        objects.add(object);
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return REPORT;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Report report;

            if (xo.hasAttribute("type")) {
                if (xo.getAttribute("type").equals("TEXT"))
                    report = new Report();
                else if (xo.getAttribute("type").equals("XHTML"))
                    report = new XHTMLReport();
                else
                    throw new XMLParseException("unknown document type, " + xo.getAttribute("type") + ", for report");
            } else
                report = new Report();

            report.setTitle(xo.getAttribute("title", ""));

            for (int i = 0; i < xo.getChildCount(); i++) {
                Object child = xo.getChild(i);
                report.add(child);
            }

            report.setOutput(XMLParser.getFilePrintWriter(xo, getParserName()));
            report.createReport();

            return report;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Generates a report using the given text and elements";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new StringAttributeRule("type", "The format of the report", new String[]{"TEXT", "XHTML"}, true),
                new StringAttributeRule("title", "The title of the report", "Report", true),
                new ElementRule(Object.class, "An arbitrary mixture of text and elements to report", 1, Integer.MAX_VALUE),
                AttributeRule.newStringRule(FILENAME, true),
        };

        public Class getReturnType() {
            return Report.class;
        }
    };

}
