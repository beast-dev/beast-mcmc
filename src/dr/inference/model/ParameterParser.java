/*
 * ParameterParser.java
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

package dr.inference.model;

import java.io.File;
import java.io.IOException;
import java.util.List;

import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceException;
import dr.math.MathUtils;
import dr.math.distributions.Distribution;
import dr.xml.*;

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
    public static final String RANDOMIZE = "randomize";
    public static final String FILENAME = "fileName";
    public static final String BURNIN = "burnin";
    public static final String PARAMETERCOLUMN = "parameterColumn";

    public String getParserName() {
        return PARAMETER;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double[] values = null;
        double[] uppers;
        double[] lowers;

        if( xo.hasAttribute(DIMENSION) ) {
            values = new double[xo.getIntegerAttribute(DIMENSION)];
        }
        
        if ( xo.hasAttribute(FILENAME)) {
        	//read samples from a file (used for a reference prior) and calculate the mean
        	String fileName = xo.getStringAttribute(FILENAME);
        	File file = new File(fileName);
        	fileName = file.getName();
        	String parent = file.getParent();
        	if (!file.isAbsolute()) {
				parent = System.getProperty("user.dir");
			}
        	file = new File(parent, fileName);
        	fileName = file.getAbsolutePath();
        	String columnName = "";
        	if ( xo.hasAttribute(PARAMETERCOLUMN)) {
        		columnName = xo.getStringAttribute(PARAMETERCOLUMN);
        	} else {
        		throw new XMLParseException("when providing a file name you must provide a parameter column as well");
        	}
        	//if a number for the burnin is not given, use 0 ...
        	int burnin;
        	if ( xo.hasAttribute(BURNIN)) {
        		burnin = xo.getIntegerAttribute(BURNIN);
        	} else {
        		burnin = 0;
        	}
        	LogFileTraces traces = new LogFileTraces(fileName, file);
        	List parameterSamples = null;
			try {
				traces.loadTraces();
				traces.setBurnIn(burnin);
				int traceIndexParameter = -1;
				for (int i = 0; i < traces.getTraceCount(); i++) {
					String traceName = traces.getTraceName(i);
					if (traceName.trim().equals(columnName)) {
						traceIndexParameter = i;
					}
				}
				if (traceIndexParameter == -1) {
					throw new XMLParseException("Column '" + columnName + "' can not be found for " + getParserName() + " element.");
				}
				parameterSamples = traces.getValues(traceIndexParameter);
			} catch (TraceException e) {
				throw new XMLParseException(e.getMessage());
			} catch (IOException ioe) {
				throw new XMLParseException(ioe.getMessage());
			}
        	values = new double[1];
        	for (int i = 0, stop = parameterSamples.size(); i < stop; i++) {
        		values[0] += ((Double)parameterSamples.get(i))/((double)stop);
        	}
        	System.out.println("Number of samples: " + parameterSamples.size());
        	System.out.println("Parameter mean: " + values[0]);
        } else if( xo.hasAttribute(VALUE) ) {
            if( values == null ) {
                values = xo.getDoubleArrayAttribute(VALUE);
            } else {
                double[] v = xo.getDoubleArrayAttribute(VALUE);
                if( v.length == values.length ) {
                    System.arraycopy(v, 0, values, 0, v.length);
                } else if( v.length == 1 ) {
                    for(int i = 0; i < values.length; i++) {
                        values[i] = v[0];
                    }
                } else {
                    throw new XMLParseException("value string must have 1 value or dimension values");
                }
            }
        } else {
            if( xo.hasAttribute(DIMENSION) ) {
                values = new double[xo.getIntegerAttribute(DIMENSION)];
            } else {
                // parameter dimension will get set correctly by TreeModel presumably.
                if (!xo.hasChildNamed(RANDOMIZE)) {
                    return new Parameter.Default(1);
                }
                values = new double[1];
                values[0] = 1.0;
            }
        }

        uppers = new double[values.length];
        for(int i = 0; i < values.length; i++) {
            uppers[i] = Double.POSITIVE_INFINITY;
        }

        lowers = new double[values.length];
        for(int i = 0; i < values.length; i++) {
            lowers[i] = Double.NEGATIVE_INFINITY;
        }

        if( xo.hasAttribute(UPPER) ) {
            double[] v = xo.getDoubleArrayAttribute(UPPER);
            if( v.length == uppers.length ) {
                System.arraycopy(v, 0, uppers, 0, v.length);
            } else if( v.length == 1 ) {
                for(int i = 0; i < uppers.length; i++) {
                    uppers[i] = v[0];
                }
            } else {
                throw new XMLParseException("uppers string must have 1 value or dimension values");
            }
        }

        if( xo.hasAttribute(LOWER) ) {
            double[] v = xo.getDoubleArrayAttribute(LOWER);
            if( v.length == lowers.length ) {
                System.arraycopy(v, 0, lowers, 0, v.length);
            } else if( v.length == 1 ) {
                for(int i = 0; i < lowers.length; i++) {
                    lowers[i] = v[0];
                }
            } else {
                throw new XMLParseException("lowers string must have 1 value or dimension values");
            }
        }

        //  assert uppers != null && lowers != null;

        if( (uppers.length != values.length) ) {
            throw new XMLParseException("value and upper limit strings have different dimension, in parameter");
        }

        if( (lowers.length != values.length) ) {
            throw new XMLParseException("value and lower limit strings have different dimension, in parameter");
        }

        // check if uppers and lowers are consistent
        for(int i = 0; i < values.length; i++) {
            if( uppers[i] < lowers[i] ) {
                throw new XMLParseException("upper is lower than lower, in parameter");
            }
        }

        if (xo.hasChildNamed(RANDOMIZE)) {
                     
            Distribution distribution = (Distribution) xo.getChild(RANDOMIZE).getChild(Distribution.class);
            for (int i = 0; i < values.length; i++) {
                do {
                    // Not an efficient way to draw random variables, but this is currently the only general interface
                    values[i] = distribution.quantile(MathUtils.nextDouble());
                } while (values[i] < lowers[i] || values[i] > uppers[i]);
            }

        } else {

            // make values consistent with bounds
            for(int i = 0; i < values.length; i++) {
                if( uppers[i] < values[i] ) values[i] = uppers[i];
            }

            for(int i = 0; i < values.length; i++) {
                if (lowers[i] > values[i]) values[i] = lowers[i];
            }
        }

        Parameter param = new Parameter.Default(values.length);
        for(int i = 0; i < values.length; i++) {
            param.setParameterValue(i, values[i]);
        }
        param.addBounds(new Parameter.DefaultBounds(uppers, lowers));
        return param;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleArrayRule(VALUE, true),
            AttributeRule.newIntegerRule(DIMENSION, true),
            AttributeRule.newStringRule(FILENAME, true),
            AttributeRule.newStringRule(PARAMETERCOLUMN, true),
            AttributeRule.newIntegerRule(BURNIN, true),
            AttributeRule.newDoubleArrayRule(UPPER, true),
            AttributeRule.newDoubleArrayRule(LOWER, true),
            new ElementRule(RANDOMIZE, new XMLSyntaxRule[] {
                    new ElementRule(Distribution.class),
            },true),
    };


    public String getParserDescription() {
        return "A real-valued parameter of one or more dimensions.";
    }

    public Class getReturnType() {
        return Parameter.class;
    }

    static public void replaceParameter(XMLObject xo, Parameter newParam) throws XMLParseException {

        for (int i = 0; i < xo.getChildCount(); i++) {

            if (xo.getChild(i) instanceof Parameter) {

                XMLObject rxo;
                Object obj = xo.getRawChild(i);

                if (obj instanceof Reference ) {
                    rxo = ((Reference) obj).getReferenceObject();
                } else if (obj instanceof XMLObject) {
                    rxo = (XMLObject) obj;
                } else {
                    throw new XMLParseException("object reference not available");
                }

                if (rxo.getChildCount() > 0) {
                    throw new XMLParseException("No child elements allowed in parameter element.");
                }

                if (rxo.hasAttribute(XMLParser.IDREF)) {
                    throw new XMLParseException("References to " + xo.getName() + " parameters are not allowed in treeModel.");
                }

                if (rxo.hasAttribute(VALUE)) {
                    throw new XMLParseException("Parameters in " + xo.getName() + " have values set automatically.");
                }

                if (rxo.hasAttribute(UPPER)) {
                    throw new XMLParseException("Parameters in " + xo.getName() + " have bounds set automatically.");
                }

                if (rxo.hasAttribute(LOWER)) {
                    throw new XMLParseException("Parameters in " + xo.getName() + " have bounds set automatically.");
                }

                if (rxo.hasAttribute(XMLParser.ID)) {
                    newParam.setId(rxo.getStringAttribute(XMLParser.ID));
                }

                rxo.setNativeObject(newParam);

                return;
            }
        }
    }

//    static public Parameter getParameter(XMLObject xo) throws XMLParseException {
//
//        int paramCount = 0;
//        Parameter param = null;
//        for (int i = 0; i < xo.getChildCount(); i++) {
//            if (xo.getChild(i) instanceof Parameter) {
//                param = (Parameter) xo.getChild(i);
//                paramCount += 1;
//            }
//        }
//
//        if (paramCount == 0) {
//            throw new XMLParseException("no parameter element in treeModel " + xo.getName() + " element");
//        } else if (paramCount > 1) {
//            throw new XMLParseException("More than one parameter element in treeModel " + xo.getName() + " element");
//        }
//
//        return param;
//    }
}
