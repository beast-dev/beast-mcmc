/*
 * Transform.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.xml.*;

/**
 * interface for the one-to-one transform of a continuous variable.
 * A static member Transform.LOG provides an instance of LogTransform
 *
 * @author Andrew Rambaut
 * @author Guy Baele
 * @author Marc Suchard
 * @version $Id: Transform.java,v 1.5 2005/05/24 20:26:01 rambaut Exp $
 */
public interface Transform {
    /**
     * @return the transformed value
     */
    double transform(double value);

    /**
     * @return the inverse transformed value
     */
    double inverse(double value);

    /**
     * @return the transform's name
     */
    String getTransformName();


    public double getLogJacobian(double oldValue, double newValue);


    public static class LogTransform implements Transform {

        public LogTransform() {
        }

        public double transform(double value) {
            return Math.log(value);
        }

        public double inverse(double value) {
            return Math.exp(value);
        }

        public String getTransformName() {
            return "log";
        }

        public double getLogJacobian(double oldValue, double newValue) {
            return Math.log(newValue * (1.0 / oldValue));
        }
    }

    public static class LogitTransform implements Transform {

        public LogitTransform() {
        }

        public double transform(double value) {
            return Math.log(value / (1.0 - value));
        }

        public double inverse(double value) {
            return 1.0 / (1.0 + Math.exp(-value));
        }

        public String getTransformName() {
            return "logit";
        }

        public double getLogJacobian(double oldValue, double newValue) {
            return Math.log((1.0 - newValue) * (1.0 / oldValue - 1.0 / (1.0 - oldValue)));
        }
    }

    public static class NoTransform implements Transform {

        public NoTransform() {
        }

        public double transform(double value) {
            return value;
        }

        public double inverse(double value) {
            return value;
        }

        public String getTransformName() {
            return "none";
        }

        public double getLogJacobian(double oldValue, double newValue) {
            return 0.0;
        }
    }

    public class ParsedTransform {
        public Transform transform;
        public int start; // zero-indexed
        public int end; // zero-indexed
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Transform thisTransform = Transform.NONE;
            String name = (String) xo.getAttribute(TYPE);
            for (Transform type : Transform.transformList) {
                if (name.equals(type.getTransformName())) {
                    thisTransform = type;
                    break;
                }
            }

            ParsedTransform transform = new ParsedTransform();
            transform.transform = thisTransform;
            transform.start = Integer.parseInt((String) xo.getAttribute(START)) - 1;
            transform.end = Integer.parseInt((String) xo.getAttribute(END));
            return transform;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    AttributeRule.newStringRule(TYPE),
                    AttributeRule.newIntegerRule(START),
                    AttributeRule.newIntegerRule(END),
            };
        }

        public String getParserDescription() {
            return null;
        }

        public Class getReturnType() {
            return ParsedTransform.class;
        }

        public String getParserName() {
            return TRANSFORM;
        }
    };


    public static final LogTransform LOG = new LogTransform();
    public static final LogitTransform LOGIT = new LogitTransform();
    public static final NoTransform NONE = new NoTransform();

    public static final Transform[] transformList = {LOG, LOGIT, NONE};

    public static final String TRANSFORM = "transform";
    public static final String TYPE = "type";
    public static final String START = "start";
    public static final String END = "end";
}
