/*
 * Transform.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

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
     * @param value evaluation point
     * @return the transformed value
     */
    double transform(double value);

    /**
     * overloaded transformation that takes and returns an array of doubles
     * @param values evaluation points
     * @param from start transformation at this index
     * @param to end transformation at this index
     * @return the transformed values
     */
    double[] transform(double[] values, int from, int to);

    /**
     * @param value evaluation point
     * @return the inverse transformed value
     */
    double inverse(double value);

    /**
     * overloaded transformation that takes and returns an array of doubles
     * @param values evaluation points
     * @param from start transformation at this index
     * @param to end transformation at this index
     * @return the transformed values
     */
    double[] inverse(double[] values, int from, int to);

    /**
     * @return the transform's name
     */
    String getTransformName();

    /**
     * @param value evaluation point
     * @return the log of the transform's jacobian
     */
    double getLogJacobian(double value);

    /**
     * @param values evaluation points
     * @param from start calculation at this index
     * @param to end calculation at this index
     * @return the log of the transform's jacobian
     */
    double getLogJacobian(double[] values, int from, int to);


    class LogTransform implements Transform {

        public LogTransform() {
        }

        public double transform(double value) {
            return Math.log(value);
        }

        public double[] transform(double[] values, int from, int to) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

        public double inverse(double value) {
            return Math.exp(value);
        }

        public double[] inverse(double[] values, int from, int to) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

        public String getTransformName() {
            return "log";
        }

        public double getLogJacobian(double value) {
            return -Math.log(value);
        }

        public double getLogJacobian(double[] values, int from, int to) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }
    }

    class LogConstrainedSumTransform implements Transform {

        public LogConstrainedSumTransform() {
        }

        public double transform(double value) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

        public double[] transform(double[] values, int from, int to) {
            double[] transformedValues = new double[to - from + 1];
            int counter = 0;
            for (int i = from; i <= to; i++) {
                transformedValues[counter] = Math.log(values[i]);
                counter++;
            }
            return transformedValues;
        }

        public double inverse(double value) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

        //inverse transformation assumes a sum of elements equal to the number of elements
        public double[] inverse(double[] values, int from, int to) {
            double sum = (double)(to - from + 1);
            double[] transformedValues = new double[to - from + 1];
            int counter = 0;
            double newSum = 0.0;
            for (int i = from; i <= to; i++) {
                transformedValues[counter] = Math.exp(values[i]);
                newSum += transformedValues[counter];
                counter++;
            }
            for (int i = 0; i < sum; i++) {
                transformedValues[i] = (transformedValues[i] / newSum) * sum;
            }
            return transformedValues;
        }

        public String getTransformName() {
            return "logConstrainedSum";
        }

        public double getLogJacobian(double value) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

        public double getLogJacobian(double[] values, int from, int to) {
            double sum = 0.0;
            for (int i = from; i <= to; i++) {
                sum -= Math.log(values[i]);
            }
            return sum;
        }

        public static void main(String[] args) {

            //specify starting values
            double[] startValues = {1.5, 0.6, 0.9};
            System.err.print("Starting values: ");
            double startSum = 0.0;
            for (int i = 0; i < startValues.length; i++) {
                System.err.print(startValues[i] + " ");
                startSum += startValues[i];
            }
            System.err.println("\nSum = " + startSum);

            //perform transformation
            double[] transformedValues = LOG_CONSTRAINED_SUM.transform(startValues, 0, startValues.length-1);
            System.err.print("Transformed values: ");
            for (int i = 0; i < transformedValues.length; i++) {
                System.err.print(transformedValues[i] + " ");
            }
            System.err.println();

            //add draw for normal distribution to transformed elements
            for (int i = 0; i < transformedValues.length; i++) {
                transformedValues[i] += 0.20 * Math.random();
            }

            //perform inverse transformation
            transformedValues = LOG_CONSTRAINED_SUM.inverse(transformedValues, 0, transformedValues.length-1);
            System.err.print("New values: ");
            double endSum = 0.0;
            for (int i = 0; i < transformedValues.length; i++) {
                System.err.print(transformedValues[i] + " ");
                endSum += transformedValues[i];
            }
            System.err.println("\nSum = " + endSum);

            if (startSum != endSum) {
                System.err.println("Starting and ending constraints differ!");
            }

        }

    }

    class LogitTransform implements Transform {

        public LogitTransform() {
            range = 1.0;
            lower = 0.0;
        }

        public double transform(double value) {
            return Math.log(value / (1.0 - value));
        }

        public double[] transform(double[] values, int from, int to) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

        public double inverse(double value) {
            return 1.0 / (1.0 + Math.exp(-value));
        }

        public double[] inverse(double[] values, int from, int to) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

        public String getTransformName() {
            return "logit";
        }

        public double getLogJacobian(double value) {
            return -Math.log(1.0 - value) - Math.log(value);
        }

        public double getLogJacobian(double[] values, int from, int to) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

        private final double range;
        private final double lower;
    }

    class FisherZTransform implements Transform {

        public FisherZTransform() {
        }

        public double transform(double value) {
            return 0.5 * (Math.log(1.0 + value) - Math.log(1.0 - value));
        }

        public double[] transform(double[] values, int from, int to) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

        public double inverse(double value) {
            return (Math.exp(2 * value) - 1) / (Math.exp(2 * value) + 1);
        }

        public double[] inverse(double[] values, int from, int to) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

        public String getTransformName() {
            return "fisherz";
        }

        public double getLogJacobian(double value) {
            return -Math.log(1 - value) - Math.log(1 + value);
        }

        public double getLogJacobian(double[] values, int from, int to) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

    }

    class NoTransform implements Transform {

        public NoTransform() {
        }

        public double transform(double value) {
            return value;
        }

        public double[] transform(double[] values, int from, int to) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

        public double inverse(double value) {
            return value;
        }

        public double[] inverse(double[] values, int from, int to) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

        public String getTransformName() {
            return "none";
        }

        public double getLogJacobian(double value) {
            return 0.0;
        }

        public double getLogJacobian(double[] values, int from, int to) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

    }

    class ParsedTransform {
        public Transform transform;
        public int start; // zero-indexed
        public int end; // zero-indexed, i.e, i = start; i < end; ++i
        public int every;
        public List<Parameter> parameters = null;
    }

    XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Transform thisTransform = Transform.NONE;
            String name = (String) xo.getAttribute(TYPE);
            System.err.println("name: " + name);

            thisTransform = Type.valueOf(name).transform;

            ParsedTransform transform = new ParsedTransform();
            transform.transform = thisTransform;
            if (xo.hasAttribute(START)) {
                transform.start = xo.getIntegerAttribute(START);
                transform.end = xo.getAttribute(END, Integer.MAX_VALUE);
                transform.every = xo.getAttribute(EVERY, 1);
                // todo: check values are valid
                transform.start--; // zero-indexed
            } else {
                transform.parameters = new ArrayList<Parameter>();

                for (Parameter param : xo.getAllChildren(Parameter.class)) {
                    transform.parameters.add(param);
                }
            }

            return transform;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    AttributeRule.newStringRule(TYPE),
                    AttributeRule.newIntegerRule(START, true),
                    AttributeRule.newIntegerRule(END, true),
                    AttributeRule.newIntegerRule(EVERY, true),
                    new ElementRule(Parameter.class, 0, Integer.MAX_VALUE)
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

    class Util {
        public static Transform[] getListOfNoTransforms(int size) {
            Transform[] transforms = new Transform[size];
            for (int i = 0; i < size; ++i) {
                transforms[i] = NONE;
            }
            return transforms;
        }
    }

    NoTransform NONE = new NoTransform();
    LogTransform LOG = new LogTransform();
    LogConstrainedSumTransform LOG_CONSTRAINED_SUM = new LogConstrainedSumTransform();
    LogitTransform LOGIT = new LogitTransform();
    FisherZTransform FISHER_Z = new FisherZTransform();

    enum Type {
        NONE("none", new NoTransform()),
        LOG("log", new LogTransform()),
        LOG_CONSTRAINED_SUM("logConstrainedSum", new LogConstrainedSumTransform()),
        LOGIT("logit", new LogitTransform()),
        FISHER_Z("fisherZ",new FisherZTransform());

        Type(String name, Transform transform) {
            this.name = name;
            this.transform = transform;
        }

        public Transform getTransform() {
            return transform;
        }

        public String getName() {
            return name;
        }

        private Transform transform;
        private String name;
    }

    String TRANSFORM = "transform";
    String TYPE = "type";
    String START = "start";
    String END = "end";
    String EVERY = "every";

}
