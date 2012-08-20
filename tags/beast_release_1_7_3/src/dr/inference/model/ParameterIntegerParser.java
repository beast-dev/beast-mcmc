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
 * @author Walter Xie
 *
 * @version $Id: ParameterParser.java,v 1.12 2005/05/24 20:26:00 rambaut Exp $
 */
public class ParameterIntegerParser extends dr.xml.AbstractXMLObjectParser {

//    public static final String UPPER = "upper";
//    public static final String LOWER = "lower";
    public static final String DIMENSION = "dimension";
    public static final String VALUE = "value";
    public static final String PARAMETER = "integerParameter";
//    public static final String RANDOMIZE = "randomize";

    public String getParserName() {
        return PARAMETER;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        int[] values = new int[0];
//        double[] uppers;
//        double[] lowers;

//        if( xo.hasAttribute(DIMENSION) ) {
//            values = new int[xo.getIntegerAttribute(DIMENSION)];
//        }

        if( xo.hasAttribute(VALUE) ) {
            if( values == null ) {
                values = xo.getIntegerArrayAttribute(VALUE);
            } else {
                int[] v = xo.getIntegerArrayAttribute(VALUE);
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
                values = new int[xo.getIntegerAttribute(DIMENSION)];
            } else {
                // parameter dimension will get set correctly by TreeModel presumably.
//                if (!xo.hasChildNamed(RANDOMIZE)) {
//                    return new Parameter.Default(1);
//                }
                values = new int[1];
                values[0] = 0;
            }
        }

//        uppers = new double[values.length];
//        for(int i = 0; i < values.length; i++) {
//            uppers[i] = Double.POSITIVE_INFINITY;
//        }
//
//        lowers = new double[values.length];
//        for(int i = 0; i < values.length; i++) {
//            lowers[i] = Double.NEGATIVE_INFINITY;
//        }

//        if( xo.hasAttribute(UPPER) ) {
//            double[] v = xo.getDoubleArrayAttribute(UPPER);
//            if( v.length == uppers.length ) {
//                System.arraycopy(v, 0, uppers, 0, v.length);
//            } else if( v.length == 1 ) {
//                for(int i = 0; i < uppers.length; i++) {
//                    uppers[i] = v[0];
//                }
//            } else {
//                throw new XMLParseException("uppers string must have 1 value or dimension values");
//            }
//        }
//
//        if( xo.hasAttribute(LOWER) ) {
//            double[] v = xo.getDoubleArrayAttribute(LOWER);
//            if( v.length == lowers.length ) {
//                System.arraycopy(v, 0, lowers, 0, v.length);
//            } else if( v.length == 1 ) {
//                for(int i = 0; i < lowers.length; i++) {
//                    lowers[i] = v[0];
//                }
//            } else {
//                throw new XMLParseException("lowers string must have 1 value or dimension values");
//            }
//        }

        //  assert uppers != null && lowers != null;

//        if( (uppers.length != values.length) ) {
//            throw new XMLParseException("value and upper limit strings have different dimension, in parameter");
//        }
//
//        if( (lowers.length != values.length) ) {
//            throw new XMLParseException("value and lower limit strings have different dimension, in parameter");
//        }
//
//        // check if uppers and lowers are consistent
//        for(int i = 0; i < values.length; i++) {
//            if( uppers[i] < lowers[i] ) {
//                throw new XMLParseException("upper is lower than lower, in parameter");
//            }
//        }
//
//        if (xo.hasChildNamed(RANDOMIZE)) {
//
//            Distribution distribution = (Distribution) xo.getChild(RANDOMIZE).getChild(Distribution.class);
//            for (int i = 0; i < values.length; i++) {
//                do {
//                    // Not an efficient way to draw random variables, but this is currently the only general interface
//                    values[i] = distribution.quantile(MathUtils.nextDouble());
//                } while (values[i] < lowers[i] || values[i] > uppers[i]);
//            }
//
//        } else {
//
//            // make values consistent with bounds
//            for(int i = 0; i < values.length; i++) {
//                if( uppers[i] < values[i] ) values[i] = uppers[i];
//            }
//
//            for(int i = 0; i < values.length; i++) {
//                if (lowers[i] > values[i]) values[i] = lowers[i];
//            }
//        }

        Variable<Integer> param = new Variable.I(values);

        param.addBounds(new Bounds.Staircase(param));
        return param;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerArrayRule(VALUE, true),
            AttributeRule.newIntegerRule(DIMENSION),
//            AttributeRule.newDoubleArrayRule(UPPER, true),
//            AttributeRule.newDoubleArrayRule(LOWER, true),
//            new ElementRule(RANDOMIZE, new XMLSyntaxRule[] {
//                    new ElementRule(Distribution.class),
//            },true),
    };


    public String getParserDescription() {
        return "An integer-valued parameter only for staircase bound.";
    }

    public Class getReturnType() {
        return Variable.class;
    }

//    static public void replaceParameter(XMLObject xo, Parameter newParam) throws XMLParseException {
//
//        for (int i = 0; i < xo.getChildCount(); i++) {
//
//            if (xo.getChild(i) instanceof Parameter) {
//
//                XMLObject rxo;
//                Object obj = xo.getRawChild(i);
//
//                if (obj instanceof Reference ) {
//                    rxo = ((Reference) obj).getReferenceObject();
//                } else if (obj instanceof XMLObject) {
//                    rxo = (XMLObject) obj;
//                } else {
//                    throw new XMLParseException("object reference not available");
//                }
//
//                if (rxo.getChildCount() > 0) {
//                    throw new XMLParseException("No child elements allowed in parameter element.");
//                }
//
//                if (rxo.hasAttribute(XMLParser.IDREF)) {
//                    throw new XMLParseException("References to " + xo.getName() + " parameters are not allowed in treeModel.");
//                }
//
//                if (rxo.hasAttribute(VALUE)) {
//                    throw new XMLParseException("Parameters in " + xo.getName() + " have values set automatically.");
//                }
//
////                if (rxo.hasAttribute(UPPER)) {
////                    throw new XMLParseException("Parameters in " + xo.getName() + " have bounds set automatically.");
////                }
////
////                if (rxo.hasAttribute(LOWER)) {
////                    throw new XMLParseException("Parameters in " + xo.getName() + " have bounds set automatically.");
////                }
//
//                if (rxo.hasAttribute(XMLParser.ID)) {
//                    newParam.setId(rxo.getStringAttribute(XMLParser.ID));
//                }
//
//                rxo.setNativeObject(newParam);
//
//                return;
//            }
//        }
//    }

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
