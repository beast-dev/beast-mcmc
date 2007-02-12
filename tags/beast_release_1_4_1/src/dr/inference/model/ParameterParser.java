/*
 * ParameterParser.java
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

package dr.inference.model;

import dr.xml.AttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Parses a multi-dimensional continuous parameter.
 * 
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 * @version $Id: ParameterParser.java,v 1.12 2005/05/24 20:26:00 rambaut Exp $
 */
public class ParameterParser extends dr.xml.AbstractXMLObjectParser {
	
	public static final String UPPER = "upper";
	public static final String LOWER = "lower";
	public static final String DIMENSION = "dimension";
	public static final String VALUE = "value";	
	public static final String PARAMETER = "parameter";	
		
	public String getParserName() { return PARAMETER; }
	
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		
		double[] values = null;
		double[] uppers = null;
		double[] lowers = null;
		
		if (xo.hasAttribute(DIMENSION)) {
				values = new double[xo.getIntegerAttribute(DIMENSION)];
		} 
		
		if (xo.hasAttribute(VALUE)) {
			if (values == null) {
				values = xo.getDoubleArrayAttribute(VALUE);
			} else {
				double[] v = xo.getDoubleArrayAttribute(VALUE);
				if (v.length == values.length) {
					System.arraycopy(v, 0, values, 0, v.length);
				} else if (v.length == 1) {
					for (int i = 0; i < values.length; i++) {
						values[i] = v[0];
					}
				} else {
					throw new XMLParseException("value string must have 1 value or dimension values");
				}
			}
		} else {
			if (xo.hasAttribute(DIMENSION)) {
				values = new double[xo.getIntegerAttribute(DIMENSION)];
			} else {
				// parameter dimension will get set correctly by TreeModel presumably.
				return new Parameter.Default(1);
			}
		}

        uppers = new double[values.length];
        for (int i =0; i < values.length; i++) {
            uppers[i] = Double.POSITIVE_INFINITY;
		}
		
        lowers = new double[values.length];
        for (int i =0; i < values.length; i++) {
            lowers[i] = Double.NEGATIVE_INFINITY;
        }

		if (xo.hasAttribute(UPPER)) {
            double[] v = xo.getDoubleArrayAttribute(UPPER);
            if (v.length == uppers.length) {
                System.arraycopy(v, 0, uppers, 0, v.length);
            } else if (v.length == 1) {
                for (int i = 0; i < uppers.length; i++) {
                    uppers[i] = v[0];
                }
            } else {
                throw new XMLParseException("uppers string must have 1 value or dimension values");
            }
		}
		
		if (xo.hasAttribute(LOWER)) {
            double[] v = xo.getDoubleArrayAttribute(LOWER);
            if (v.length == lowers.length) {
                System.arraycopy(v, 0, lowers, 0, v.length);
            } else if (v.length == 1) {
                for (int i = 0; i < lowers.length; i++) {
                    lowers[i] = v[0];
                }
            } else {
                throw new XMLParseException("lowers string must have 1 value or dimension values");
            }
		}

		if (uppers != null && (uppers.length != values.length)) {
			throw new XMLParseException("value and upper limit strings have different dimension, in parameter");
		}

		if (lowers != null && (lowers.length != values.length)) {
			throw new XMLParseException("value and lower limit strings have different dimension, in parameter");
		}
		
		// check if uppers and lowers are consistent
		for (int i =0; i < values.length; i++) {
			if (uppers[i] < lowers[i]) {
				throw new XMLParseException("upper is lower than lower, in parameter");
			}
		}
		
		// make values consistent with bounds
		for (int i =0; i < values.length; i++) {
			if (uppers[i] < values[i]) values[i] = uppers[i];
		}
					
		Parameter param = new Parameter.Default(values.length);
		for (int i =0; i < values.length; i++) {
			param.setParameterValue(i, values[i]); 
		}
		param.addBounds(new Parameter.DefaultBounds(uppers, lowers));
		return param;
	}
	
			public XMLSyntaxRule[] getSyntaxRules() { return rules; }
			
			private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
				AttributeRule.newDoubleArrayRule(VALUE, true),
				AttributeRule.newIntegerRule(DIMENSION, true),
				AttributeRule.newDoubleArrayRule(UPPER, true),
				AttributeRule.newDoubleArrayRule(LOWER, true)
			};


	public String getParserDescription() {
		return "A real-valued parameter of one or more dimensions.";
	}
	
	public Class getReturnType() { return Parameter.class; }
}
