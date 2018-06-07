/*
 * TransformParsers.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.util;

import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;

/**
 * @author Marc A. Suchard
 */

public class TransformParsers {

    @SuppressWarnings("unused")
    public static XMLObjectParser COMPOUND_PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Parameter parameter = (Parameter) xo.getChild(Parameter.class);
            return new Transform.Collection(xo.getAllChildren(Transform.ParsedTransform.class),
                    parameter);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[] {
                    new ElementRule(Transform.ParsedTransform.class, 1, Integer.MAX_VALUE),
                    new ElementRule(Parameter.class),
            };
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return Transform.Collection.class;
        }

        @Override
        public String getParserName() {
            return COMPOUND;
        }
    };

    @SuppressWarnings("unused")
    public static XMLObjectParser COMPOSE_PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            XMLObject outerXo = xo.getChild(OUTER);
            XMLObject innerXo = xo.getChild(INNER);

            Transform.ParsedTransform outerPT = (Transform.ParsedTransform)
                    outerXo.getChild(Transform.ParsedTransform.class);

            Transform.ParsedTransform innerPT = (Transform.ParsedTransform)
                    innerXo.getChild(Transform.ParsedTransform.class);

            if (!outerPT.equivalent(innerPT)) {
                throw new XMLParseException("Not equivalent transformations");
            }

            if (outerPT.transform instanceof Transform.UnivariableTransform &&
                    innerPT.transform instanceof Transform.UnivariableTransform) {

                Transform.ParsedTransform composition = outerPT.clone();
                composition.transform = new Transform.Compose((Transform.UnivariableTransform) outerPT.transform,
                        (Transform.UnivariableTransform) innerPT.transform);
                return composition;

            } else {
                throw new XMLParseException("Not composable transform types");
            }
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[] {
                    new ElementRule(OUTER, new XMLSyntaxRule[] {
                            new ElementRule(Transform.ParsedTransform.class),
                    }),
                    new ElementRule(INNER, new XMLSyntaxRule[] {
                            new ElementRule(Transform.ParsedTransform.class),
                    }),
            };
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return Transform.ParsedTransform.class;
        }

        @Override
        public String getParserName() {
            return COMPOSE;
        }
    };

    @SuppressWarnings("unused")
    public static XMLObjectParser INVERSE_PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Transform.ParsedTransform parsedTransform = (Transform.ParsedTransform)
                    xo.getChild(Transform.ParsedTransform.class);

            Transform transform = parsedTransform.transform;

            if (transform instanceof Transform.UnivariableTransform) {

                Transform.ParsedTransform inverseTransform = parsedTransform.clone();
                inverseTransform.transform = new Transform.Inverse((Transform.UnivariableTransform)
                        parsedTransform.transform);

                return inverseTransform;
            } else {
                throw new XMLParseException("Not invertible transform type");
            }
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[] {
                    new ElementRule(Transform.ParsedTransform.class),
            };
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return Transform.ParsedTransform.class;
        }

        @Override
        public String getParserName() {
            return INVERSE;
        }
    };

    @SuppressWarnings("unused")
    public static XMLObjectParser TRANSFORM_PARSER = new AbstractXMLObjectParser() {

         public Object parseXMLObject(XMLObject xo) throws XMLParseException {

             Transform thisTransform;
             String name = (String) xo.getAttribute(TYPE);
//             System.err.println("name: " + name);

             thisTransform = null;
             for (Transform.Type type: Transform.Type.values()) {
                 if (name.equalsIgnoreCase(type.getName())) {
                     thisTransform = type.getTransform();
                 }
             }
             if (thisTransform == null) {
                 throw new XMLParseException("Unrecognized transform type, " + name);
             }

             if (xo.getAttribute(INVERSE, false)) {
                 if (thisTransform instanceof Transform.UnivariableTransform) {
                     thisTransform = new Transform.Inverse((Transform.UnivariableTransform)thisTransform);
                 } else {
                     throw new XMLParseException("Non-invertible transform type, " + name);
                 }
             }

             Transform.ParsedTransform transform = new Transform.ParsedTransform();
             transform.transform = thisTransform;
             if (xo.hasAttribute(START)) {
                 transform.start = xo.getIntegerAttribute(START);
                 transform.end = xo.getAttribute(END, Integer.MAX_VALUE);
                 transform.every = xo.getAttribute(EVERY, 1);
                 // todo: check values are valid
                 transform.start--; // zero-indexed
             } else {
                 if (xo.hasAttribute(SUM)) {
                     transform.fixedSum = xo.getDoubleAttribute(SUM);
                 }
                 transform.parameters = new ArrayList<Parameter>();

                 transform.parameters.addAll(xo.getAllChildren(Parameter.class));
             }

             return transform;
         }

         public XMLSyntaxRule[] getSyntaxRules() {
             return new XMLSyntaxRule[]{
                     AttributeRule.newStringRule(TYPE),
                     AttributeRule.newIntegerRule(START, true),
                     AttributeRule.newIntegerRule(END, true),
                     AttributeRule.newIntegerRule(EVERY, true),
                     AttributeRule.newDoubleRule(SUM, true),
                     AttributeRule.newBooleanRule(INVERSE, true),
                     new ElementRule(Transform.ParsedTransform.class, 0, 1),
                     new ElementRule(Parameter.class, 0, Integer.MAX_VALUE)
             };
         }

         public String getParserDescription() {
             return null;
         }

         public Class getReturnType() {
             return Transform.ParsedTransform.class;
         }

         public String getParserName() {
             return TRANSFORM;
         }
     };


    public static final String TRANSFORM = "transform";
    public static final String TYPE = "type";
    public static final String START = "start";
    public static final String END = "end";
    public static final String SUM = "sum";
    public static final String EVERY = "every";
    private static final String INVERSE = "inverseTransform";
    private static final String COMPOSE = "composedTransform";
    private static final String OUTER = "outer";
    private static final String INNER = "inner";
    private static final String COMPOUND = "compoundTransform";
}
