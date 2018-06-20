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

import dr.inference.model.Parameter;
import dr.math.MathUtils;

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
     * overloaded transformation that takes and returns an array of doubles
     * @param values evaluation points
     * @param from start transformation at this index
     * @param to end transformation at this index
     * @param sum fixed sum of values that needs to be enforced
     * @return the transformed values
     */
    double[] inverse(double[] values, int from, int to, double sum);

    double updateGradientLogDensity(double gradient, double value);

    double[] updateGradientLogDensity(double[] gradient, double[] value, int from, int to);

    double gradientInverse(double value);

    double[] gradientInverse(double[] values, int from, int to);

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

    abstract class UnivariableTransform implements Transform {

        public abstract double transform(double value);

        public double[] transform(double[] values, int from, int to) {
            double[] result = values.clone();
            for (int i = from; i < to; ++i) {
                result[i] = transform(values[i]);
            }
            return result;
        }

        public abstract double inverse(double value);

        public double[] inverse(double[] values, int from, int to) {
            double[] result = values.clone();
            for (int i = from; i < to; ++i) {
                result[i] = inverse(values[i]);
            }
            return result;
        }

        public double[] inverse(double[] values, int from, int to, double sum) {
            throw new RuntimeException("Fixed sum cannot be enforced for a univariate transformation.");
        }

        public abstract double gradientInverse(double value);

        public double[] gradientInverse(double[] values, int from, int to) {
            double[] result = values.clone();
            for (int i = from; i < to; ++i) {
                result[i] = gradientInverse(values[i]);
            }
            return result;
        }

        public abstract double updateGradientLogDensity(double gradient, double value);

        public double[] updateGradientLogDensity(double[] gradient, double[] value , int from, int to) {
            double[] result = value.clone();
            for (int i = from; i < to; ++i) {
                result[i] = updateGradientLogDensity(gradient[i], value[i]);
            }
            return result;
        }

        public abstract double getLogJacobian(double value);

        public double getLogJacobian(double[] values, int from, int to) {
            double sum = 0.0;
            for (int i = from; i < to; ++i) {
                sum += getLogJacobian(values[i]);
            }
            return sum;
        }
    }

    abstract class MultivariableTransform implements Transform {

        public double transform(double value) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

        public double inverse(double value) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

        public double updateGradientLogDensity(double gradient, double value) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

        public double gradientInverse(double value) {
             throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
         }

        public double getLogJacobian(double value) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }
    }

    abstract class MultivariableTransformWithParameter extends MultivariableTransform {
        abstract public Parameter getParameter();
    }

    class LogTransform extends UnivariableTransform {

        public double transform(double value) {
            return Math.log(value);
        }

        public double inverse(double value) {
            return Math.exp(value);
        }

        public double gradientInverse(double value) { return Math.exp(value); }

        public double updateGradientLogDensity(double gradient, double value) {
            // value == gradient of inverse()
            // 1.0 == gradient of log Jacobian of inverse()
            return gradient * value + 1.0;
        }

        public String getTransformName() { return "log"; }

        public double getLogJacobian(double value) { return -Math.log(value); }
    }

    class LogConstrainedSumTransform extends MultivariableTransform {

        //private double fixedSum;

        public LogConstrainedSumTransform() {
        }

        /*public LogConstrainedSumTransform(double fixedSum) {
            this.fixedSum = fixedSum;
        }

        public double getConstrainedSum() {
            return this.fixedSum;
        }*/

        public double[] transform(double[] values, int from, int to) {
            double[] transformedValues = new double[to - from + 1];
            int counter = 0;
            for (int i = from; i <= to; i++) {
                transformedValues[counter] = Math.log(values[i]);
                counter++;
            }
            return transformedValues;
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
            /*for (int i = 0; i < sum; i++) {
                transformedValues[i] = (transformedValues[i] / newSum) * sum;
            }*/
            for (int i = 0; i < transformedValues.length; i++) {
                transformedValues[i] = (transformedValues[i] / newSum) * sum;
            }
            return transformedValues;
        }

        //inverse transformation assumes a given sum provided as an argument
        public double[] inverse(double[] values, int from, int to, double sum) {
            //double sum = (double)(to - from + 1);
            double[] transformedValues = new double[to - from + 1];
            int counter = 0;
            double newSum = 0.0;
            for (int i = from; i <= to; i++) {
                transformedValues[counter] = Math.exp(values[i]);
                newSum += transformedValues[counter];
                counter++;
            }
            /*for (int i = 0; i < sum; i++) {
                transformedValues[i] = (transformedValues[i] / newSum) * sum;
            }*/
            for (int i = 0; i < transformedValues.length; i++) {
                transformedValues[i] = (transformedValues[i] / newSum) * sum;
            }
            return transformedValues;
        }

        public String getTransformName() {
            return "logConstrainedSum";
        }

        public double[] updateGradientLogDensity(double[] gradient, double[] value, int from, int to) {
            throw new RuntimeException("Not yet implemented");
        }

        public double[] gradientInverse(double[] values, int from, int to) {
            throw new RuntimeException("Not yet implemented");
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
                transformedValues[i] += 0.20 * MathUtils.nextDouble();
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

    class LogitTransform extends UnivariableTransform {

        public LogitTransform() {
            range = 1.0;
            lower = 0.0;
        }

        public double transform(double value) {
            return Math.log(value / (1.0 - value));
        }

        public double inverse(double value) {
            return 1.0 / (1.0 + Math.exp(-value));
        }

        public double gradientInverse(double value) {
            throw new RuntimeException("Not yet implemented");
        }

        public double updateGradientLogDensity(double gradient, double value) {
            throw new RuntimeException("Not yet implemented");
        }

        public String getTransformName() {
            return "logit";
        }

        public double getLogJacobian(double value) {
            return -Math.log(1.0 - value) - Math.log(value);
        }

        private final double range;
        private final double lower;
    }

    class FisherZTransform extends UnivariableTransform {

        public double transform(double value) {
            return 0.5 * (Math.log(1.0 + value) - Math.log(1.0 - value));
        }

        public double inverse(double value) {
            return (Math.exp(2 * value) - 1) / (Math.exp(2 * value) + 1);
        }

        public double gradientInverse(double value) {
            throw new RuntimeException("Not yet implemented");
        }

        public double updateGradientLogDensity(double gradient, double value) {
            throw new RuntimeException("Not yet implemented");
        }

        public String getTransformName() {
            return "fisherz";
        }

        public double getLogJacobian(double value) {
            return -Math.log(1 - value) - Math.log(1 + value);
        }
    }

    class NegateTransform extends UnivariableTransform {

        public double transform(double value) {
            return -value;
        }

        public double inverse(double value) {
            return -value;
        }

        public double updateGradientLogDensity(double gradient, double value) {
            // -1 == gradient of inverse()
            // 0.0 == gradient of log Jacobian of inverse()
            return -gradient;
        }

        public double gradientInverse(double value) { return -1.0; }

        public String getTransformName() {
            return "negate";
        }

        public double getLogJacobian(double value) {
            return 0.0;
        }
    }

    class PowerTransform extends UnivariableTransform{
        private double power;

        PowerTransform(){
            this.power = 2;
        }

        PowerTransform(double power){
            this.power = power;
        }

        @Override
        public String getTransformName() {
            return "Power Transform";
        }

        @Override
        public double transform(double value) {
            return Math.pow(value, power);
        }

        @Override
        public double inverse(double value) {
            return Math.pow(value, 1 / power);
        }

        @Override
        public double gradientInverse(double value) {
            throw new RuntimeException("not implemented yet");
//            return 0;
        }

        @Override
        public double updateGradientLogDensity(double gradient, double value) {
            throw new RuntimeException("not implemented yet");
        }

        @Override
        public double getLogJacobian(double value) {
            throw new RuntimeException("not implemented yet");
        }
    }

    class NoTransform extends UnivariableTransform {

        public double transform(double value) {
            return value;
        }

        public double inverse(double value) {
            return value;
        }

        public double updateGradientLogDensity(double gradient, double value) {
            return gradient;
        }

        public double gradientInverse(double value) { return 1.0; }

        public String getTransformName() {
            return "none";
        }

        public double getLogJacobian(double value) {
            return 0.0;
        }
    }

    class Compose extends UnivariableTransform  {

        public Compose(UnivariableTransform outer, UnivariableTransform inner) {
            this.outer = outer;
            this.inner = inner;
        }

        @Override
        public String getTransformName() {
            return "compose." + outer.getTransformName() + "." + inner.getTransformName();
        }

        @Override
        public double transform(double value) {
            final double outerValue = inner.transform(value);
            final double outerTransform = outer.transform(outerValue);

//            System.err.println(value + " " + outerValue + " " + outerTransform);
//            System.exit(-1);

            return outerTransform;
//            return outer.transform(inner.transform(value));
        }

        @Override
        public double inverse(double value) {
            return inner.inverse(outer.inverse(value));
        }

        @Override
        public double gradientInverse(double value) {
            return inner.gradientInverse(value) * outer.gradientInverse(inner.transform(value));
        }

        @Override
        public double updateGradientLogDensity(double gradient, double value) {
//            final double innerGradient = inner.updateGradientLogDensity(gradient, value);
//            final double outerValue = inner.transform(value);
//            final double outerGradient = outer.updateGradientLogDensity(innerGradient, outerValue);
//            return outerGradient;

            return outer.updateGradientLogDensity(inner.updateGradientLogDensity(gradient, value), inner.transform(value));
        }

        @Override
        public double getLogJacobian(double value) {
            return inner.getLogJacobian(value) + outer.getLogJacobian(inner.transform(value));
        }

        private final UnivariableTransform outer;
        private final UnivariableTransform inner;
    }

    class Inverse extends UnivariableTransform {

        public Inverse(UnivariableTransform inner) {
            this.inner = inner;
        }

        @Override
        public String getTransformName() {
            return "inverse." + inner.getTransformName();
        }

        @Override
        public double transform(double value) {
            return inner.inverse(value);  // Purposefully switched

        }

        @Override
        public double updateGradientLogDensity(double gradient, double value) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double inverse(double value) {
            return inner.transform(value); // Purposefully switched
        }

        @Override
        public double gradientInverse(double value) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double getLogJacobian(double value) {
            return -inner.getLogJacobian(value);
        }

        private final UnivariableTransform inner;
    }


    class Array extends MultivariableTransformWithParameter {

          private final List<Transform> array;
          private final Parameter parameter;

          public Array(List<Transform> array, Parameter parameter) {
              this.parameter = parameter;
              this.array = array;

//              if (parameter.getDimension() != array.size()) {
//                  throw new IllegalArgumentException("Dimension mismatch");
//              }
          }

          public Parameter getParameter() { return parameter; }

          @Override
          public double[] transform(double[] values, int from, int to) {

              final double[] result = values.clone();

              for (int i = from; i < to; ++i) {
                  result[i] = array.get(i).transform(values[i]);
              }
              return result;
          }

          @Override
          public double[] inverse(double[] values, int from, int to) {

              final double[] result = values.clone();

              for (int i = from; i < to; ++i) {
                  result[i] = array.get(i).inverse(values[i]);
              }
              return result;
          }

          @Override
          public double[] inverse(double[] values, int from, int to, double sum) {
              throw new RuntimeException("Not yet implemented.");
          }

          @Override
          public double[] gradientInverse(double[] values, int from, int to) {

              final double[] result = values.clone();

              for (int i = from; i < to; ++i) {
                  result[i] = array.get(i).gradientInverse(values[i]);
              }
              return result;
          }

          @Override
          public double[] updateGradientLogDensity(double[] gradient, double[] values, int from, int to) {

              final double[] result = values.clone();

              for (int i = from; i < to; ++i) {
                  result[i] = array.get(i).updateGradientLogDensity(gradient[i], values[i]);
              }
              return result;
          }

          @Override
          public String getTransformName() {
              return "array";
          }

          @Override
          public double getLogJacobian(double[] values, int from, int to) {

              double sum = 0.0;

              for (int i = from; i < to; ++i) {
                  sum += array.get(i).getLogJacobian(values[i]);
              }
              return sum;
          }
    }

    class Collection extends MultivariableTransformWithParameter {

        private final List<ParsedTransform> segments;
        private final Parameter parameter;

        public Collection(List<ParsedTransform> segments, Parameter parameter) {
            this.parameter = parameter;
            this.segments = ensureContiguous(segments);
        }

        public Parameter getParameter() { return parameter; }

        private List<ParsedTransform> ensureContiguous(List<ParsedTransform> segments) {

            final List<ParsedTransform> contiguous = new ArrayList<ParsedTransform>();

            int current = 0;
            for (ParsedTransform segment : segments) {
                if (current < segment.start) {
                    contiguous.add(new ParsedTransform(NONE, current, segment.start));
                }
                contiguous.add(segment);
                current = segment.end;
            }
            if (current < parameter.getDimension()) {
                contiguous.add(new ParsedTransform(NONE, current, parameter.getDimension()));
            }

//            System.err.println("Segments:");
//            for (ParsedTransform transform : contiguous) {
//                System.err.println(transform.transform.getTransformName() + " " + transform.start + " " + transform.end);
//            }
//            System.exit(-1);

            return contiguous;
        }

        @Override
        public double[] transform(double[] values, int from, int to) {

            final double[] result = values.clone();

            for (ParsedTransform segment : segments) {
                if (from < segment.end && to >= segment.start) {
                    final int begin = Math.max(segment.start, from);
                    final int end = Math.min(segment.end, to);
                    for (int i = begin; i < end; ++i) {
                        result[i] = segment.transform.transform(values[i]);
                    }
                }
            }
            return result;
        }

        @Override
        public double[] inverse(double[] values, int from, int to) {

            final double[] result = values.clone();

            for (ParsedTransform segment : segments) {
                if (from < segment.end && to >= segment.start) {
                    final int begin = Math.max(segment.start, from);
                    final int end = Math.min(segment.end, to);
                    for (int i = begin; i < end; ++i) {
                        result[i] = segment.transform.inverse(values[i]);
                    }
                }
            }
            return result;
        }

        @Override
        public double[] inverse(double[] values, int from, int to, double sum) {
            throw new RuntimeException("Not yet implemented.");
        }

        @Override
        public double[] gradientInverse(double[] values, int from, int to) {

            final double[] result = values.clone();

            for (ParsedTransform segment : segments) {
                if (from < segment.end && to >= segment.start) {
                    final int begin = Math.max(segment.start, from);
                    final int end = Math.min(segment.end, to);
                    for (int i = begin; i < end; ++i) {
                        result[i] = segment.transform.gradientInverse(values[i]);
                    }
                }
            }
            return result;
        }

        @Override
        public double[] updateGradientLogDensity(double[] gradient, double[] values, int from, int to) {

            final double[] result = values.clone();

            for (ParsedTransform segment : segments) {
                if (from < segment.end && to >= segment.start) {
                    final int begin = Math.max(segment.start, from);
                    final int end = Math.min(segment.end, to);
                    for (int i = begin; i < end; ++i) {
                        result[i] = segment.transform.updateGradientLogDensity(gradient[i], values[i]);
                    }
                }
            }
            return result;
        }

        @Override
        public String getTransformName() {
            return "collection";
        }

        @Override
        public double getLogJacobian(double[] values, int from, int to) {

            double sum = 0.0;

            for (ParsedTransform segment : segments) {
                if (from < segment.end && to >= segment.start) {
                    final int begin = Math.max(segment.start, from);
                    final int end = Math.min(segment.end, to);
                    for (int i = begin; i < end; ++i) {
                        sum += segment.transform.getLogJacobian(values[i]);
                    }
                }
            }
//            System.err.println("Log: " + sum + " " + segments.size());
            return sum;
        }

//        class Segment {
//
//            public Segment(Transform transform, int start, int end) {
//                this.transform = transform;
//                this.start = start;
//                this.end = end;
//            }
//            public Transform transform;
//            public int start;
//            public int end;
//        }
    }

    class ParsedTransform {
        public Transform transform;
        public int start; // zero-indexed
        public int end; // zero-indexed, i.e, i = start; i < end; ++i
        public int every = 1;
        public double fixedSum = 0.0;
        public List<Parameter> parameters = null;

        public ParsedTransform() {
            
        }

        public ParsedTransform(Transform transform, int start, int end) {
            this.transform = transform;
            this.start = start;
            this.end = end;
        }

        public ParsedTransform clone() {
            ParsedTransform clone = new ParsedTransform();
            clone.transform = transform;
            clone.start = start;
            clone.end = end;
            clone.every = every;
            clone.fixedSum = fixedSum;
            clone.parameters = parameters;
            return clone;
        }

        public boolean equivalent(ParsedTransform other) {
            if (start == other.start && end == other.end && every == other.every && parameters == other.parameters) {
                return true;
            } else {
                return false;
            }
        }
    }

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
    NegateTransform NEGATE = new NegateTransform();
    Compose LOG_NEGATE = new Compose(new LogTransform(), new NegateTransform());
    LogConstrainedSumTransform LOG_CONSTRAINED_SUM = new LogConstrainedSumTransform();
    LogitTransform LOGIT = new LogitTransform();
    FisherZTransform FISHER_Z = new FisherZTransform();

    enum Type {
        NONE("none", new NoTransform()),
        LOG("log", new LogTransform()),
        NEGATE("negate", new NegateTransform()),
        LOG_NEGATE("log-negate", new Compose(new LogTransform(), new NegateTransform())),
        LOG_CONSTRAINED_SUM("logConstrainedSum", new LogConstrainedSumTransform()),
        LOGIT("logit", new LogitTransform()),
        FISHER_Z("fisherZ",new FisherZTransform()),
        POWER("power", new PowerTransform());

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

//    String TRANSFORM = "transform";
//    String TYPE = "type";
//    String START = "start";
//    String END = "end";
//    String EVERY = "every";
//    String INVERSE = "inverse";

}
