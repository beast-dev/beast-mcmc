/*
 * Report.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

import java.util.ArrayList;


/**
 * A generates a text report from the elements within it.
 *
 * @version $Id: Report.java,v 1.15 2005/05/24 20:26:01 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class Report {
	
	public static final String REPORT = "report";

	protected String title = "";
	protected ArrayList<Object> objects = new ArrayList<Object>();
	
	public void createReport() {
		System.out.println(getTitle());
		System.out.println();

        for (Object object : objects) {
            final String item = object.toString();
            System.out.print(item.trim());
            System.out.print(" ");
        }
        System.out.println();
	}

	/** set the report's title. */
	public void setTitle(String title) { this.title = title; }

	/** @return the report's title. */
	public String getTitle() { return title; }
	

	/** add an Object element. */
	public void add(Object object) {
		objects.add(object); 
	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
	
		public String getParserName() { return REPORT; }
	
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
			
			if (xo.hasAttribute("title")) {
				report.setTitle(xo.getStringAttribute("title"));
			}
			
			for (int i = 0; i < xo.getChildCount(); i++) {
				Object child = xo.getChild(i);
				report.add(child);
			}
						
			report.createReport();
			
			return report;
		}
		
		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "Generates a report using the given text and elements";
		}
		
		public XMLSyntaxRule[] getSyntaxRules() { return rules; }
		
		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new StringAttributeRule("type", "The format of the report", new String[] { "TEXT", "XHTML" }, true),
			new StringAttributeRule("title", "The title of the report", "Report", true),
			new ElementRule(Object.class, "An arbitrary mixture of text and elements to report", 1, Integer.MAX_VALUE)
		};
	
		public Class getReturnType() { return Report.class; }
	};

}
