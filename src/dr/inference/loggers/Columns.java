/*
 * Columns.java
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

package dr.inference.loggers;

import dr.util.Identifiable;
import dr.xml.*;

import java.util.ArrayList;

/**
 * A class which parses column elements
 *
 * @version $Id: Columns.java,v 1.10 2005/05/24 20:25:59 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class Columns {

	public static final String COLUMN = "column";
	public static final String LABEL = "label";
	public static final String SIGNIFICANT_FIGURES = "sf";
	public static final String DECIMAL_PLACES = "dp";
	public static final String WIDTH = "width";
    public static final String FORMAT = "format";
    public static final String PERCENT = "percent";
    public static final String BOOL = "boolean";

    public Columns(LogColumn[] columns) { this.columns = columns; }
	public LogColumn[] getColumns() { return columns; }

	/**
	 * @return an object based on the XML element it was passed.
	 */
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return COLUMN; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			String label = null;
//			int sf = -1;
//			int dp = -1;
//			int width = -1;

			if (xo.hasAttribute(LABEL)) {
				label = xo.getStringAttribute(LABEL);
			}

            final int sf = xo.getAttribute(SIGNIFICANT_FIGURES, -1);
            final int dp = xo.getAttribute(DECIMAL_PLACES, -1);
            final int width = xo.getAttribute(WIDTH, -1);

//            if (xo.hasAttribute(SIGNIFICANT_FIGURES)) {
//				sf = xo.getIntegerAttribute(SIGNIFICANT_FIGURES);
//			}
//
//			if (xo.hasAttribute(DECIMAL_PLACES)) {
//				dp = xo.getIntegerAttribute(DECIMAL_PLACES);
//			}
//
//			if (xo.hasAttribute(WIDTH)) {
//				width = xo.getIntegerAttribute(WIDTH);
//			}

            String format = xo.getAttribute(FORMAT, "");

            ArrayList colList = new ArrayList();

			for (int i = 0; i < xo.getChildCount(); i++) {

				Object child = xo.getChild(i);
				LogColumn[] cols;

				if (child instanceof Loggable) {
					cols = ((Loggable)child).getColumns();
				} else if (child instanceof Identifiable) {
					cols = new LogColumn[] { new LogColumn.Default(((Identifiable)child).getId(), child) };
				} else {
					cols = new LogColumn[] { new LogColumn.Default(child.getClass().toString(), child) };
				}

                if( format.equals(PERCENT) ) {
                    for(int k = 0; k < cols.length; ++k) {
                        if( cols[k] instanceof NumberColumn ) {
                            cols[k] = new PercentColumn((NumberColumn)cols[k]);
                        }
                    }
                } else if( format.equals(BOOL) ) {
                    for(int k = 0; k < cols.length; ++k) {
                        if( cols[k] instanceof NumberColumn ) {
                            cols[k] = new BooleanColumn((NumberColumn)cols[k]);
                        }
                    }
                }

                for (int j = 0; j < cols.length; j++) {

					if (label != null) {
						if (cols.length > 1) {
							cols[j].setLabel(label + Integer.toString(j+1));
						} else {
							cols[j].setLabel(label);
						}
					}

					if (cols[j] instanceof NumberColumn) {
						if (sf != -1) {
							((NumberColumn)cols[j]).setSignificantFigures(sf);
						}
						if (dp != -1) {
							((NumberColumn)cols[j]).setDecimalPlaces(dp);
						}
					}

					if (width > 0) {
						cols[j].setMinimumWidth(width);
					}

					colList.add(cols[j]);
				}
			}

			LogColumn[] columns = new LogColumn[colList.size()];
			colList.toArray(columns);

			return new Columns(columns);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "Specifies formating options for one or more columns in a log file.";
		}

		public Class getReturnType() { return Columns.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
                new StringAttributeRule(LABEL,
                        "The label of the column. " +
                                "If this is specified and more than one statistic is in this column, " +
                                "then the label will be appended by the index of the statistic to create individual column names", true),
                AttributeRule.newIntegerRule(SIGNIFICANT_FIGURES, true),
                AttributeRule.newIntegerRule(DECIMAL_PLACES, true),
                AttributeRule.newIntegerRule(WIDTH, true),
                AttributeRule.newStringRule(FORMAT, true),
                // Anything goes???
                new ElementRule(Object.class, 1, Integer.MAX_VALUE),
        };

    };

	private LogColumn[] columns = null;
}
