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
import dr.math.matrixAlgebra.Matrix;
import dr.xml.XMLObject;
import org.apache.commons.math.util.FastMath;

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

    double updateDiagonalHessianLogDensity(double diagonalHessian, double gradient, double value);

    double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] value, int from, int to);

    double[][] updateHessianLogDensity(double[][] hessian, double[][] transformationHessian, double[] gradient, double[] value, int from, int to);

    double updateOffdiagonalHessianLogDensity(double offdiagonalHessian, double transformationHessian, double gradientI, double gradientJ, double valueI, double valueJ);

    double updateGradientInverseUnWeightedLogDensity(double gradient, double value);

    double[] updateGradientInverseUnWeightedLogDensity(double[] gradient, double[] value, int from, int to);

    double updateGradientUnWeightedLogDensity(double gradient, double value);

    double[] updateGradientUnWeightedLogDensity(double[] gradient, double[] value, int from, int to);

    double gradient(double value);

    double[] gradient(double[] values, int from, int to);

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

    /**
     * @return true if the transform is multivatiate (i.e. components not independents)
     */
    boolean isMultivariate();

    /**
     * @return false if the point is on the frontier of the domain
     */
    boolean isInInteriorDomain(double value);

    boolean isInInteriorDomain(double[] values, int from, int to);

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

        public double updateGradientLogDensity(double gradient, double value) {
            // value : untransformed. TODO:use updateGradientUnWeightedLogDensity()
            return updateGradientInverseUnWeightedLogDensity(gradient, transform(value)) + getGradientLogJacobianInverse(transform(value));
        }

        public double[] updateGradientLogDensity(double[] gradient, double[] value , int from, int to) {
            double[] result = value.clone();
            for (int i = from; i < to; ++i) {
                result[i] = updateGradientLogDensity(gradient[i], value[i]);
            }
            return result;
        }

        public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] value, int from, int to) {
            double[] result = value.clone();
            for (int i = from; i < to; ++i) {
                result[i] = updateDiagonalHessianLogDensity(diagonalHessian[i], gradient[i], value[i]);
            }
            return result;
        }

        public double updateGradientInverseUnWeightedLogDensity(double gradient, double value) {
            // value is transformed
            return gradient * gradientInverse(value);
        }

        public double updateGradientUnWeightedLogDensity(double gradient, double value) {
            // value is unTransformed
            return gradient * gradient(value);
        }

        public double[] updateGradientInverseUnWeightedLogDensity(double[] gradient, double[] value, int from, int to) {
            double[] result = value.clone();
            for (int i = from; i < to; ++i) {
                result[i] = updateGradientInverseUnWeightedLogDensity(gradient[i], value[i]);
            }
            return result;
        }

        public double[] updateGradientUnWeightedLogDensity(double[] gradient, double[] value, int from, int to) {
            double[] result = value.clone();
            for (int i = from; i < to; ++i) {
                result[i] = updateGradientUnWeightedLogDensity(gradient[i], value[i]);
            }
            return result;
        }

        public double[][] updateHessianLogDensity(double[][] hessian, double[][] transformationHessian, double[] gradient, double[] value, int from, int to) {

            final int dim = to - from;
            double[][] updatedHessian = new double[dim][dim];

            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    if (i == j) updatedHessian[i][j] = updateDiagonalHessianLogDensity(hessian[i][j], gradient[i], value[i]);
                    else updatedHessian[i][j] = updateOffdiagonalHessianLogDensity(hessian[i][j], transformationHessian[i][j], gradient[i], gradient[j], value[i], value[j]);
                }
            }
            return updatedHessian;
        }

        protected abstract double getGradientLogJacobianInverse(double value); // takes transformed value

        public abstract double gradient(double value);

        @Override
        public double[] gradient(double[] values, int from, int to) {
            double[] result = values.clone();
            for (int i = from; i < to; ++i) {
                result[i] = gradient(values[i]);
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

        public boolean isMultivariate() { return false;}

        public abstract boolean isInInteriorDomain(double value);

        public boolean isInInteriorDomain(double[] values, int from, int to) {
            for (double val : values) {
                if (!isInInteriorDomain(val)) return false;
            }
            return true;
        }
    }

    abstract class MultivariableTransform implements Transform {

        abstract public int getDimension();

        abstract public int getInputDimension();

        abstract public int getOutputDimension();

        public double transform(double value) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

        public double inverse(double value) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

        public double updateGradientLogDensity(double gradient, double value) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

        public double updateDiagonalHessianLogDensity(double diagonalHessian, double gradient, double value) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

        public double updateOffdiagonalHessianLogDensity(double offdiagonalHessian, double transformationHessian, double gradientI, double gradientJ, double valueI, double valueJ) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

        public double updateGradientInverseUnWeightedLogDensity(double gradient, double value) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

        public double updateGradientUnWeightedLogDensity(double gradient, double value) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

        protected double getGradientLogJacobianInverse(double value) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

        public double gradientInverse(double value) {
             throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
         }

        public double gradient(double value) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

        public double getLogJacobian(double value) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }

        public boolean isInInteriorDomain(double value) {
            throw new RuntimeException("Transformation not permitted for this type of parameter, exiting ...");
        }
    }

    abstract class MultivariableTransformWithParameter extends MultivariableTransform {
        abstract public Parameter getParameter();
    }

    abstract class MultivariateTransform extends MultivariableTransform {
        // A class for a multivariate transform

        protected int dim;

        protected int inputDimension;

        protected int outputDimension;

        public MultivariateTransform(int dim){
            this(dim, dim);
        }

        public MultivariateTransform(int inputDimension, int outputDimension) {
            this.inputDimension = inputDimension;
            this.outputDimension = outputDimension;
            this.dim = outputDimension;
        }

        @Override
        public int getInputDimension() {
            return inputDimension;
        }

        @Override
        public int getOutputDimension() {
            return outputDimension;
        }

        @Override
        public int getDimension() {
            return dim;
        }

        protected abstract double[] transform(double[] values);

        @Override
        public final double[] transform(double[] values, int from, int to) {
            assert from == 0 && to == values.length && dim == values.length
                    : "The multivariate transform function can only be applied to the whole array of values.";
            return transform(values);
        }

        protected abstract double[] inverse(double[] values);

        @Override
        public final double[] inverse(double[] values, int from, int to) {
            assert from == 0 && to == values.length && dim == values.length
                    : "The multivariate transform function can only be applied to the whole array of values.";
            return inverse(values);
        }

        protected abstract double getLogJacobian(double[] values);

        @Override
        public final double getLogJacobian(double[] values, int from, int to) {
            assert from == 0 && to == values.length && dim == values.length
                    : "The multivariate transform function can only be applied to the whole array of values.";
            return getLogJacobian(values);
        }

        @Override
        public final double[] updateGradientLogDensity(double[] gradient, double[] values, int from, int to) {
            assert from == 0 && to == values.length && dim == values.length
                    : "The multivariate transform function can only be applied to the whole array of values.";
            return updateGradientLogDensity(gradient, values);
        }

        protected double[] updateGradientLogDensity(double[] gradient, double[] value) {
            // values = untransformed (R)
            double[] transformedValues = transform(value, 0, value.length);
            // Transform Inverse
            double[] updatedGradient = updateGradientInverseUnWeightedLogDensity(gradient, transformedValues);
            // gradient of log jacobian of the inverse
            double[] gradientLogJacobianInverse = getGradientLogJacobianInverse(transformedValues);
            // Add gradient log jacobian
            for (int i = 0; i < gradient.length; i++) {
                updatedGradient[i] += gradientLogJacobianInverse[i];
            }
            return updatedGradient;
        }

        public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] value, int from, int to) {
            throw new RuntimeException("Not yet implemented");
        }

        public double[][] updateHessianLogDensity(double[][] hessian, double[][] transformationHessian, double[] gradient, double[] value, int from, int to) {
            throw new RuntimeException("Not yet implemented");
        }

        protected double[] updateGradientInverseUnWeightedLogDensity(double[] gradient, double[] value) {
            // takes transformed values
            // Jacobian of inverse (transpose)
            double[][] jacobianInverse = computeJacobianMatrixInverse(value);
            return updateGradientJacobian(gradient, jacobianInverse);
        }

        public double[] updateGradientInverseUnWeightedLogDensity(double[] gradient, double[] value, int from, int to) {
            assert from == 0 && to == value.length && dim == value.length
                    : "The multivariate transform function can only be applied to the whole array of values.";
            return updateGradientInverseUnWeightedLogDensity(gradient, value);
        }

        double[] updateGradientJacobian(double[] gradient, double[][] jacobianInverse) {
            // Matrix multiplication
            double[] updatedGradient = new double[gradient.length];
            for (int i = 0; i < gradient.length; i++) {
                for (int j = 0; j < gradient.length; j++) {
                    updatedGradient[i] += jacobianInverse[i][j] * gradient[j];
                }
            }
            return updatedGradient;
        }

        @Override
        public double[] updateGradientUnWeightedLogDensity(double[] gradient, double[] value, int from, int to) {
            // takes untransformed value TODO: more efficient way ?
            return updateGradientInverseUnWeightedLogDensity(gradient, transform(value, from, to), from, to);
        }

        abstract protected double[] getGradientLogJacobianInverse(double[] values); // transformed value

        abstract public double[][] computeJacobianMatrixInverse(double[] values); // transformed values

        public boolean isMultivariate() { return true;}

        protected abstract boolean isInInteriorDomain(double[] values);

        @Override
        public final boolean isInInteriorDomain(double[] values, int from, int to) {
            assert from == 0 && to == values.length && dim == values.length
                    : "The multivariate transform function can only be applied to the whole array of values.";
            return isInInteriorDomain(values);
        }
    }

    abstract class MatrixVariateTransform extends MultivariateTransform {

        protected final int rowDimension;
        protected final int columnDimension;

        public MatrixVariateTransform(int inputDimension, int outputRowDimension, int outputColumnDimension) {
            super(inputDimension, outputRowDimension * outputColumnDimension);
            this.rowDimension = outputRowDimension;
            this.columnDimension = outputColumnDimension;
        }


        public int getRowDimension() {
            return rowDimension;
        }

        public int getColumnDimension() {
            return columnDimension;
        }
    }

    class LogTransform extends UnivariableTransform {

        public double transform(double value) {
            return Math.log(value);
        }

        public double inverse(double value) {
            return Math.exp(value);
        }

        public boolean isInInteriorDomain(double value) {
            return value > 0.0 && !Double.isInfinite(value);
        }

        public double gradientInverse(double value) { return Math.exp(value); }

        public double updateGradientLogDensity(double gradient, double value) {
            // gradient == gradient of inverse()
            // value == gradient of inverse() (value is untransformed)
            // 1.0 == gradient of log Jacobian of inverse()
            return gradient * value + 1.0;
        }

        protected double getGradientLogJacobianInverse(double value) {
            return 1.0;
        }

        @Override
        public double updateDiagonalHessianLogDensity(double diagonalHessian, double gradient, double value) {
            // value == inverse()
            // diagonalHessian == hessian of inverse()
            // gradient == gradient of inverse()
            return value * (gradient + value * diagonalHessian);
        }

        @Override
        public double updateOffdiagonalHessianLogDensity(double offdiagonalHessian, double transfomationHessian, double gradientI, double gradientJ, double valueI, double valueJ) {
            return offdiagonalHessian * valueI * valueJ + gradientJ * transfomationHessian;
        }

        @Override
        public double gradient(double value) {
            return value;
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

        public int getDimension() {
            return -1;
        }

        @Override
        public int getInputDimension() {
            return getDimension();
        }

        @Override
        public int getOutputDimension() {
            return getDimension();
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

        public boolean isInInteriorDomain(double value) {
            return value > 0.0 && !Double.isInfinite(value);
        }

        public boolean isInInteriorDomain(double[] values, int from, int to) {
            for (double val : values) {
                if (!isInInteriorDomain(val)) return false;
            }
            return true;
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

        public double[] updateGradientInverseUnWeightedLogDensity(double[] gradient, double[] value, int from, int to) {
            throw new RuntimeException("Not yet implemented");
        }

        public double[] updateGradientUnWeightedLogDensity(double[] gradient, double[] value, int from, int to) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] value, int from, int to) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double[][] updateHessianLogDensity(double[][] hessian, double[][] transformationHessian, double[] gradient, double[] value, int from, int to) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double updateOffdiagonalHessianLogDensity(double offdiagonalHessian, double transformationHessian, double gradientI, double gradientJ, double valueI, double valueJ) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double gradient(double value) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double[] gradient(double[] values, int from, int to) {
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

        public boolean isMultivariate() { return true;}

        public static void main(String[] args) {

            //specify starting values
            double[] startValues = {1.5, 0.6, 0.9};
            System.err.print("Starting values: ");
            double startSum = 0.0;
            for (double startValue : startValues) {
                System.err.print(startValue + " ");
                startSum += startValue;
            }
            System.err.println("\nSum = " + startSum);

            //perform transformation
            double[] transformedValues = LOG_CONSTRAINED_SUM.transform(startValues, 0, startValues.length-1);
            System.err.print("Transformed values: ");
            for (double transformedValue : transformedValues) {
                System.err.print(transformedValue + " ");
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
            for (double transformedValue : transformedValues) {
                System.err.print(transformedValue + " ");
                endSum += transformedValue;
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

        public boolean isInInteriorDomain(double value) {
            return value > 0.0 && value < 1.0;
        }

        public double gradientInverse(double value) {
            return gradient(inverse(value));
        }

        public double updateGradientLogDensity(double gradient, double value) {
            return gradient * value * (1.0 - value) - (2.0 * value - 1.0);
        }

        protected double getGradientLogJacobianInverse(double value) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double updateDiagonalHessianLogDensity(double diagonalHessian, double gradient, double value) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double updateOffdiagonalHessianLogDensity(double offdiagonalHessian, double transformationHessian, double gradientI, double gradientJ, double valueI, double valueJ) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double gradient(double value) {
            return value * (1.0 - value);
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

    class ScaledLogitTransform extends UnivariableTransform {

        public ScaledLogitTransform() {
            upper = 1.0;
            lower = 0.0;
        }

        public ScaledLogitTransform(double upper, double lower) {
            this.upper = upper;
            this.lower = lower;
        }

        public double transform(double value) {
            return Math.log((value - lower) / (upper - value));
        }

        public double inverse(double value) {
            double x = Math.exp(-value);
            return (upper + lower * x) / (1.0 + x);
        }

        public boolean isInInteriorDomain(double value) {
            return value > lower && value < upper;
        }

        public double gradientInverse(double value) {
            throw new RuntimeException("Not yet implemented");
        }

        public double updateGradientLogDensity(double gradient, double value) {
            throw new RuntimeException("Not yet implemented");
        }

        protected double getGradientLogJacobianInverse(double value) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double updateDiagonalHessianLogDensity(double diagonalHessian, double gradient, double value) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double updateOffdiagonalHessianLogDensity(double offdiagonalHessian, double transformationHessian, double gradientI, double gradientJ, double valueI, double valueJ) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double gradient(double value) {
            return (value - lower) * (upper - value) / (upper - lower);
        }

        public String getTransformName() {
            return "logit";
        }

        public double getLogJacobian(double value) {
            return Math.log(upper - lower) - Math.log(upper - value) - Math.log(value - lower);
        }

        private final double upper;
        private final double lower;
    }

    class FisherZTransform extends UnivariableTransform {

        public double transform(double value) {
//            return 0.5 * (Math.log(1.0 + value) - Math.log(1.0 - value));
            return FastMath.atanh(value);
        }

        public double inverse(double value) {
//            return (Math.exp(2 * value) - 1) / (Math.exp(2 * value) + 1);
            return FastMath.tanh(value);  // optional: Math.tanh(value);
        }

        public boolean isInInteriorDomain(double value) {
            return value > -1.0 && value < 1.0;
        }

        public double gradientInverse(double value) {
            return 1.0 - Math.pow(inverse(value), 2);
        }

        public double updateGradientLogDensity(double gradient, double value) {
            // 1 - value^2 : gradient of inverse (value is untransformed)
            // - 2*value : gradient of log jacobian of inverse
            return (1.0 - value * value) * gradient  - 2 * value;
        }

        protected double getGradientLogJacobianInverse(double value) {
            // - 2*value : gradient of log jacobian of inverse (value is transformed)
            return -2 * inverse(value);
        }

        @Override
        public double updateDiagonalHessianLogDensity(double diagonalHessian, double gradient, double value) {
            return (1.0 - value * value) * (diagonalHessian * (1.0 - value * value) - 2.0 * gradient * value - 2.0);
        }

        @Override
        public double updateOffdiagonalHessianLogDensity(double offdiagonalHessian, double transformationHessian, double gradientI, double gradientJ, double valueI, double valueJ) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double gradient(double value) {
            return 1.0 - Math.pow(value, 2);
        }

        public String getTransformName() {
            return "fisherz";
        }

        public double getLogJacobian(double value) {
            return -Math.log1p(-value) - Math.log1p(value);
        }
    }

    class NegateTransform extends UnivariableTransform {

        public double transform(double value) {
            return -value;
        }

        public double inverse(double value) {
            return -value;
        }

        public boolean isInInteriorDomain(double value) {
            return true;
        }

        public double updateGradientLogDensity(double gradient, double value) {
            // -1 == gradient of inverse()
            // 0.0 == gradient of log Jacobian of inverse()
            return -gradient;
        }

        protected double getGradientLogJacobianInverse(double value) {
            return 0.0;
        }

        @Override
        public double updateDiagonalHessianLogDensity(double diagonalHessian, double gradient, double value) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double updateOffdiagonalHessianLogDensity(double offdiagonalHessian, double transformationHessian, double gradientI, double gradientJ, double valueI, double valueJ) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double gradient(double value) {
            return -1.0;
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

        public PowerTransform(double power){
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
        public boolean isInInteriorDomain(double value) {
            if (power == (int) power) {
                if (power >= 0) return true;
                return value != 0;
            }
            return value > 0;
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
        public double updateDiagonalHessianLogDensity(double diagonalHessian, double gradient, double value) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double updateOffdiagonalHessianLogDensity(double offdiagonalHessian, double transformationHessian, double gradientI, double gradientJ, double valueI, double valueJ) {
            throw new RuntimeException("Not yet implemented");
        }

        public double updateGradientInverseUnWeightedLogDensity(double gradient, double value) {
            throw new RuntimeException("not implemented yet");
        }

        @Override
        protected double getGradientLogJacobianInverse(double value) {
            throw new RuntimeException("not implemented yet");
        }

        @Override
        public double gradient(double value) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double getLogJacobian(double value) {
            throw new RuntimeException("not implemented yet");
        }
    }

    class ReciprocalTransform extends UnivariableTransform {

        @Override
        public double updateDiagonalHessianLogDensity(double diagonalHessian, double gradient, double value) {
            return 0;
        }

        @Override
        public double updateOffdiagonalHessianLogDensity(double offdiagonalHessian, double transformationHessian, double gradientI, double gradientJ, double valueI, double valueJ) {
            return 0;
        }

        @Override
        public String getTransformName() {
            return "reciprocal transform";
        }

        @Override
        public double transform(double value) {
            return 1.0 / value;
        }

        @Override
        public double inverse(double value) {
            return 1.0 / value;
        }

        @Override
        public double gradientInverse(double value) {
            throw new RuntimeException("not yet implemented");
        }

        @Override
        protected double getGradientLogJacobianInverse(double value) {
            throw new RuntimeException("not yet implemented");
        }

        @Override
        public double gradient(double value) {
            throw new RuntimeException("not yet implemented");
        }

        @Override
        public double getLogJacobian(double value) {
            throw new RuntimeException("not yet implemented");
        }

        @Override
        public boolean isInInteriorDomain(double value) {
            throw new RuntimeException("not yet implemented");
        }
    }

    class InverseSumTransform extends UnivariableTransform {
        private double sum;

        InverseSumTransform() {
            this.sum = 1;
        }

        InverseSumTransform(double sum) {
            this.sum = sum;
        }

        @Override
        public double updateDiagonalHessianLogDensity(double diagonalHessian, double gradient, double value) {
            return 0;
        }

        @Override
        public double updateOffdiagonalHessianLogDensity(double offdiagonalHessian, double transformationHessian,
                                                         double gradientI, double gradientJ, double valueI,
                                                         double valueJ) {
            return 0;
        }

        @Override
        public String getTransformName() {
            return "inversedSum transform";
        }

        @Override
        public double transform(double value) {
            return value / (value * sum - 1);
        }

        @Override
        public double inverse(double value) {
            return transform(value);
        }

        @Override
        public double gradientInverse(double value) {
            return 0;
        }

        @Override
        protected double getGradientLogJacobianInverse(double value) {
            return 0;
        }

        @Override
        public double gradient(double value) {
            return 0;
        }

        @Override
        public double getLogJacobian(double value) {
            return 0;
        }

        @Override
        public boolean isInInteriorDomain(double value) {
            return false;
        }
    }

    class NoTransform extends UnivariableTransform {

        public double transform(double value) {
            return value;
        }

        public double inverse(double value) {
            return value;
        }

        public boolean isInInteriorDomain(double value) {
            return true;
        }

        public double updateGradientLogDensity(double gradient, double value) {
            return gradient;
        }

        protected double getGradientLogJacobianInverse(double value) {
            return 0.0;
        }

        @Override
        public double updateDiagonalHessianLogDensity(double diagonalHessian, double gradient, double value) {
            return diagonalHessian;
        }

        @Override
        public double updateOffdiagonalHessianLogDensity(double offdiagonalHessian, double transformationHessian, double gradientI, double gradientJ, double valueI, double valueJ) {
            return offdiagonalHessian;
        }

        @Override
        public double gradient(double value) {
            return 1.0;
        }

        public double gradientInverse(double value) { return 1.0; }

        public String getTransformName() {
            return "none";
        }

        public double getLogJacobian(double value) {
            return 0.0;
        }
    }

    class NoTransformMultivariable extends MultivariableTransform {

        public int getDimension() {
            return -1;
        }

        @Override
        public int getInputDimension() {
            return getDimension();
        }

        @Override
        public int getOutputDimension() {
            return getDimension();
        }

        @Override
        public String getTransformName() {
            return "NoTransformMultivariate";
        }

        @Override
        public double[] transform(double[] values, int from, int to) {
            return subArray(values, from, to);
        }

        private double[] subArray(double[] values, int from, int to) {
            int length = to - from;
            if (length == values.length) return values;
            double[] result = new double[length];
            System.arraycopy(values, to, result, 0, length);
            return result;
        }

        @Override
        public double[] inverse(double[] values, int from, int to) {
            return subArray(values, from, to);
        }

        public boolean isInInteriorDomain(double[] values, int from, int to) {
            return true;
        }

        @Override
        public double[] updateGradientLogDensity(double[] gradient, double[] value, int from, int to) {
            return subArray(gradient, from, to);
        }

        @Override
        public double[] updateGradientInverseUnWeightedLogDensity(double[] gradient, double[] value, int from, int to) {
            return subArray(gradient, from, to);
        }

        @Override
        public double[] updateGradientUnWeightedLogDensity(double[] gradient, double[] value, int from, int to) {
            return subArray(gradient, from, to);
        }

        public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] value, int from, int to) {
            throw new RuntimeException("not implemented yet");
        }

        @Override
        public double[][] updateHessianLogDensity(double[][] hessian, double[][] transformationHessian, double[] gradient, double[] value, int from, int to) {
            throw new RuntimeException("not implemented yet");
        }

        @Override
        public double[] inverse(double[] values, int from, int to, double sum) {
            throw new RuntimeException("Not implemented.");
        }

        @Override
        public double[] gradient(double[] values, int from, int to) {
            return arrayValue(1.0, from, to);
        }

        private double[] arrayValue(double value, int from, int to) {
            int length = to - from;
            double[] result = new double[length];
            for (int i = 0; i < length; i++) {
                result[i] = value;
            }
            return result;
        }

        @Override
        public double[] gradientInverse(double[] values, int from, int to) {
            return arrayValue(1.0, from, to);
        }

        @Override
        public double getLogJacobian(double[] values, int from, int to) {
            return 0.0;
        }

        public boolean isMultivariate() { return false;}
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
        public boolean isInInteriorDomain(double value) {
            return inner.isInInteriorDomain(value);
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
        public double updateDiagonalHessianLogDensity(double diagonalHessian, double gradient, double value) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double updateOffdiagonalHessianLogDensity(double offdiagonalHessian, double transformationHessian, double gradientI, double gradientJ, double valueI, double valueJ) {
            throw new RuntimeException("Not yet implemented");
        }

        protected double getGradientLogJacobianInverse(double value) {
            throw new RuntimeException("not implemented yet");
        }

        @Override
        public double gradient(double value) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double getLogJacobian(double value) {
            return inner.getLogJacobian(value) + outer.getLogJacobian(inner.transform(value));
        }

        private final UnivariableTransform outer;
        private final UnivariableTransform inner;
    }

    class ComposeMultivariable extends MultivariableTransform {

        public ComposeMultivariable(MultivariableTransform outer, MultivariableTransform inner) {
            assert outer.getDimension() == inner.getDimension() : "In ComposeMultivariable, transforms should have the same dimension.";
            this.outer = outer;
            this.inner = inner;
        }

        public int getDimension() {
            return outer.getDimension();
        }

        @Override
        public int getInputDimension() {
            return inner.getInputDimension();
        }

        @Override
        public int getOutputDimension() {
            return outer.getOutputDimension();
        }

        @Override
        public String getTransformName() {
            return "compose." + outer.getTransformName() + "." + inner.getTransformName();
        }

        public Transform getInnerTransform() {
            return inner;
        }

        @Override
        public double[] transform(double[] values, int from, int to) {
            return outer.transform(inner.transform(values, from, to), from, to);
        }

        @Override
        public double[] inverse(double[] values, int from, int to) {
            return inner.inverse(outer.inverse(values, from, to), from, to);
        }

        @Override
        public double[] inverse(double[] values, int from, int to, double sum) {
            throw new RuntimeException("Not yet implemented.");
        }

        @Override
        public boolean isInInteriorDomain(double[] values, int from, int to) {
            return inner.isInInteriorDomain(values, from, to);
        }

        @Override
        public double[] updateGradientLogDensity(double[] gradient, double[] value, int from, int to) {
            return outer.updateGradientLogDensity(
                    inner.updateGradientLogDensity(gradient, value, from, to),
                    inner.transform(value, from, to),
                    from, to);
        }

        @Override
        public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] value, int from, int to) {

            return outer.updateDiagonalHessianLogDensity(
                    inner.updateDiagonalHessianLogDensity(diagonalHessian, gradient, value,from, to),
                    inner.updateGradientLogDensity(gradient, value, from, to),
                    inner.transform(value, from, to),
                    from, to);
        }

        @Override
        public double[][] updateHessianLogDensity(double[][] hessian, double[][] transformationHessian, double[] gradient, double[] value, int from, int to) {
            throw new RuntimeException("not implemented yet");
        }

        @Override
        public double[] updateGradientInverseUnWeightedLogDensity(double gradient[], double[] value, int from, int to) {
            return outer.updateGradientInverseUnWeightedLogDensity(
                    inner.updateGradientInverseUnWeightedLogDensity(gradient, outer.inverse(value, from, to), from, to),
                    value, from, to);
        }

        @Override
        public double[] updateGradientUnWeightedLogDensity(double[] gradient, double[] value, int from, int to) {
            return outer.updateGradientUnWeightedLogDensity(
                    inner.updateGradientUnWeightedLogDensity(gradient, value, from, to),
                    inner.transform(value, from, to), from, to);
        }

        @Override
        public double[] gradient(double[] values, int from, int to) {
            throw new RuntimeException("Not yet implemented.");
        }

        @Override
        public double[] gradientInverse(double[] values, int from, int to) {
            throw new RuntimeException("Not yet implemented.");
        }

        @Override
        public double getLogJacobian(double[] values, int from, int to) {
            return inner.getLogJacobian(values, from, to)
                    + outer.getLogJacobian(inner.transform(values, from, to), from, to);
        }

        public boolean isMultivariate() { return outer.isMultivariate() || inner.isMultivariate();}

        private final MultivariableTransform outer;
        private final MultivariableTransform inner;
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
        public boolean isInInteriorDomain(double value) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double updateGradientLogDensity(double gradient, double value) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double updateDiagonalHessianLogDensity(double diagonalHessian, double gradient, double value) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double updateOffdiagonalHessianLogDensity(double offdiagonalHessian, double transformationHessian, double gradientI, double gradientJ, double valueI, double valueJ) {
            throw new RuntimeException("Not yet implemented");
        }

        protected double getGradientLogJacobianInverse(double value) {
            throw new RuntimeException("not implemented yet");
        }

        @Override
        public double gradient(double value) {
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
            return -inner.getLogJacobian(inner.inverse(value));
        }

        private final UnivariableTransform inner;
    }

    class InverseMultivariate extends MultivariateTransform {

        public InverseMultivariate(MultivariateTransform inner) {
            super(inner.getDimension());
            this.inner = inner;
        }

        @Override
        public String getTransformName() {
            return "inverse." + inner.getTransformName();
        }

        @Override
        protected double[] transform(double[] values) {
            return inner.inverse(values); // Purposefully switched
        }

        @Override
        protected double[] inverse(double[] values) {
            return inner.transform(values); // Purposefully switched
        }

        @Override
        public boolean isInInteriorDomain(double[] values) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] value, int from, int to) {
            throw new RuntimeException("not implemented yet.");
        }

        @Override
        public double[][] updateHessianLogDensity(double[][] hessian, double[][] transformationHessian, double[] gradient, double[] value, int from, int to) {
            throw new RuntimeException("not implemented yet");
        }

        @Override
        public double[] updateGradientInverseUnWeightedLogDensity(double[] gradient, double[] value, int from, int to) {
            throw new RuntimeException("not implemented yet");
        }

        public double[] inverse(double[] values, int from, int to, double sum) {
            throw new RuntimeException("Not relevant.");
        }

        @Override
        public double[] gradient(double[] values, int from, int to) {
            return inner.gradientInverse(values, from, to);
        }

        @Override
        public double[] gradientInverse(double[] values, int from, int to) {
            return inner.gradient(values, from, to);
        }

        @Override
        protected double getLogJacobian(double[] values) {
            return -inner.getLogJacobian(inner.inverse(values));
        }

//        @Override
//        public double[] updateGradientLogDensity(double[] gradient, double[] value, int from, int to) {
//            throw new RuntimeException("Not yet implemented.");
//        }

//        @Override
//        public double[] updateGradientInverse(double[] gradient, double[] value, int from, int to) {
//            throw new RuntimeException("not implemented yet");
//        }

        @Override
        protected double[] updateGradientLogDensity(double[] gradient, double[] transformedValues) {
            // transformedValues = transformed
            // gradient of log jacobian of the inverse
            double[] gradientLogJacobianInverse = inner.getGradientLogJacobianInverse(transformedValues);
            // Add gradient log jacobian
            double[] updatedGradient = new double[gradient.length];
            for (int i = 0; i < gradient.length; i++) {
                updatedGradient[i] = gradient[i] - gradientLogJacobianInverse[i];
            }
            // Jacobian
            double[][] jacobian = computeJacobianMatrix(transformedValues);
            // Matrix Multiplication
            return updateGradientJacobian(updatedGradient, jacobian);
        }

        private double[][] computeJacobianMatrix(double[] transformedValues) {
            Matrix jacobianInverse = new Matrix(inner.computeJacobianMatrixInverse(transformedValues));
            return jacobianInverse.inverse().transpose().toComponents();
        }

        @Override
        public double[][] computeJacobianMatrixInverse(double[] values) {
            // values : untransformed
            Matrix jacobianInverse = new Matrix(inner.computeJacobianMatrixInverse(inner.transform(values)));
            return jacobianInverse.inverse().transpose().toComponents();
        }

        @Override
        protected double[] getGradientLogJacobianInverse(double[] transformedValues) {
            double[] gradient = inner.getGradientLogJacobianInverse(transformedValues);
            for (int i = 0; i < gradient.length; i++) {
                gradient[i] = - gradient[i];
            }
            return gradient;
        }

        private final MultivariateTransform inner;

        @Override
        public int getInputDimension() {
            return getDimension();
        }

        @Override
        public int getOutputDimension() {
            return getDimension();
        }
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

          public Array(Transform transform, int dim, Parameter parameter) {
              List<Transform> repArray = new ArrayList<Transform>();
              for (int i = 0; i < dim; i++) {
                  repArray.add(transform);
              }

              this.parameter = parameter;
              this.array = repArray;
          }

          public int getDimension() {
              return array.size();
          }

        @Override
        public int getInputDimension() {
            return getDimension();
        }

        @Override
        public int getOutputDimension() {
            return getDimension();
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
        public boolean isInInteriorDomain(double[] values, int from, int to) {
            for (int i = from; i < to; ++i) {
                if (!array.get(i).isInInteriorDomain(values[i])) return false;
            }
            return true;
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
        public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] values, int from, int to) {
            final double[] result = values.clone();

            for (int i = from; i < to; ++i) {
                result[i] = array.get(i).updateDiagonalHessianLogDensity(diagonalHessian[i], gradient[i], values[i]);
            }
            return result;
        }

        public double[] updateGradientInverseUnWeightedLogDensity(double[] gradient, double[] values, int from, int to) {

            final double[] result = values.clone();

            for (int i = from; i < to; ++i) {
                result[i] = array.get(i).updateGradientInverseUnWeightedLogDensity(gradient[i], values[i]);
            }
            return result;
        }

        public double[] updateGradientUnWeightedLogDensity(double[] gradient, double[] values, int from, int to) {

            final double[] result = values.clone();

            for (int i = from; i < to; ++i) {
                result[i] = array.get(i).updateGradientUnWeightedLogDensity(gradient[i], values[i]);
            }
            return result;
        }

        @Override
        public double[][] updateHessianLogDensity(double[][] hessian, double[][] transformationHessian, double[] gradient, double[] value, int from, int to) {

            final int dim = to - from;
            double[][] updatedHessian = new double[dim][dim];

            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    if (i == j) updatedHessian[i][j] = array.get(i).updateDiagonalHessianLogDensity(hessian[i][j], gradient[i], value[i]);
                    else {
                        assert(array.get(i).getClass().equals(array.get(j).getClass()));  // TODO: more generic implementation
                        updatedHessian[i][j] = array.get(i).updateOffdiagonalHessianLogDensity(hessian[i][j], transformationHessian[i][j], gradient[i], gradient[j], value[i], value[j]);
                    }
                }
            }
            return updatedHessian;
        }

        @Override
        public double[] gradient(double[] values, int from, int to) {

            final double[] result = values.clone();

            for (int i = from; i < to; ++i) {
                result[i] = array.get(i).gradient(values[i]);
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

        public boolean isMultivariate() { return false;}
    }

    class Collection extends MultivariableTransformWithParameter {

        private final List<ParsedTransform> segments;
        private final Parameter parameter;

        public Collection(List<ParsedTransform> segments, Parameter parameter) {
            this.parameter = parameter;
            this.segments = ensureContiguous(segments);
        }

        public int getDimension() {
            return parameter.getDimension();
        }

        @Override
        public int getInputDimension() {
            return getDimension();
        }

        @Override
        public int getOutputDimension() {
            return getDimension();
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
        public boolean isInInteriorDomain(double[] values, int from, int to) {
            for (ParsedTransform segment : segments) {
                if (from < segment.end && to >= segment.start) {
                    final int begin = Math.max(segment.start, from);
                    final int end = Math.min(segment.end, to);
                    for (int i = begin; i < end; ++i) {
                        if (!segment.transform.isInInteriorDomain(values[i])) return false;
                    }
                }
            }
            return true;
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
        public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] values, int from, int to) {

            final double[] result = values.clone();

            for (ParsedTransform segment : segments) {
                if (from < segment.end && to >= segment.start) {
                    final int begin = Math.max(segment.start, from);
                    final int end = Math.min(segment.end, to);
                    for (int i = begin; i < end; ++i) {
                        result[i] = segment.transform.updateDiagonalHessianLogDensity(diagonalHessian[i], gradient[i], values[i]);
                    }
                }
            }
            return result;
        }

        public double[] updateGradientInverseUnWeightedLogDensity(double[] gradient, double[] values, int from, int to) {

            final double[] result = values.clone();

            for (ParsedTransform segment : segments) {
                if (from < segment.end && to >= segment.start) {
                    final int begin = Math.max(segment.start, from);
                    final int end = Math.min(segment.end, to);
                    for (int i = begin; i < end; ++i) {
                        result[i] = segment.transform.updateGradientInverseUnWeightedLogDensity(gradient[i], values[i]);
                    }
                }
            }
            return result;
        }

        @Override
        public double[] updateGradientUnWeightedLogDensity(double[] gradient, double[] values, int from, int to) {

            final double[] result = values.clone();

            for (ParsedTransform segment : segments) {
                if (from < segment.end && to >= segment.start) {
                    final int begin = Math.max(segment.start, from);
                    final int end = Math.min(segment.end, to);
                    for (int i = begin; i < end; ++i) {
                        result[i] = segment.transform.updateGradientUnWeightedLogDensity(gradient[i], values[i]);
                    }
                }
            }
            return result;
        }

        @Override
        public double[][] updateHessianLogDensity(double[][] hessian, double[][] transformationHessian, double[] gradient, double[] value, int from, int to) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double[] gradient(double[] values, int from, int to) {

            final double[] result = values.clone();

            for (ParsedTransform segment : segments) {
                if (from < segment.end && to >= segment.start) {
                    final int begin = Math.max(segment.start, from);
                    final int end = Math.min(segment.end, to);
                    for (int i = begin; i < end; ++i) {
                        result[i] = segment.transform.gradient(values[i]);
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

        public boolean isMultivariate() { return false;}

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

    class MultivariateArray extends MultivariateTransform {

        private final List<MultivariableTransform> array;

        public MultivariateArray(List<MultivariableTransform> array) {
            super(getInputDimensionArray(array), getOutputDimensionArray(array));
            this.array = array;
        }

        private static int getDimensionArray(List<MultivariableTransform> array) {
            int dim = 0;
            for (MultivariableTransform anArray : array) {
                int dimArray = anArray.getDimension();
                assert dimArray > 0 : "MultivariateArray only allows for transforms with a defined dimension.";
                dim += dimArray;
            }
            return dim;
        }

        private static int getInputDimensionArray(List<MultivariableTransform> array) {
            int inputDimension = 0;
            for (MultivariableTransform anArray : array) {
                inputDimension += anArray.getInputDimension();
            }
            return inputDimension;
        }

        private static int getOutputDimensionArray(List<MultivariableTransform> array) {
            int outputDimension = 0;
            for (MultivariableTransform anArray : array) {
                outputDimension += anArray.getOutputDimension();
            }
            return outputDimension;
        }

        @Override
        protected double[] transform(double[] values) {

            final double[] result = values.clone();

            int offset = 0;
            for (MultivariableTransform anArray : array) {
                int dim = anArray.getDimension();
                double tmp[] = new double[dim];
                System.arraycopy(values, offset, tmp, 0, dim);
                System.arraycopy(anArray.transform(tmp, 0, dim), 0, result, offset, dim);
                offset += dim;
            }
            return result;
        }

        @Override
        protected double[] inverse(double[] values) {

            final double[] result = values.clone();

            int offset = 0;
            for (MultivariableTransform anArray : array) {
                int dim = anArray.getDimension();
                double tmp[] = new double[dim];
                System.arraycopy(values, offset, tmp, 0, dim);
                System.arraycopy(anArray.inverse(tmp, 0, dim), 0, result, offset, dim);
                offset += dim;
            }
            return result;
        }

        @Override
        public double[] inverse(double[] values, int from, int to, double sum) {
            throw new RuntimeException("Not yet implemented.");
        }

        @Override
        public boolean isInInteriorDomain(double[] values) {
            int offset = 0;
            for (MultivariableTransform anArray : array) {
                int dim = anArray.getDimension();
                double tmp[] = new double[dim];
                System.arraycopy(values, offset, tmp, 0, dim);
                if (!anArray.isInInteriorDomain(tmp, 0, dim)) return false;
                offset += dim;
            }
            return true;
        }

        @Override
        public double[] gradientInverse(double[] values, int from, int to) {

            final double[] result = values.clone();

            int offset = 0;
            for (MultivariableTransform anArray : array) {
                int dim = anArray.getDimension();
                double tmp[] = new double[dim];
                System.arraycopy(values, offset, tmp, 0, dim);
                System.arraycopy(anArray.gradientInverse(tmp, from, to), 0, result, offset, dim);
                offset += dim;
            }
            return result;
        }

        @Override
        protected double[] updateGradientLogDensity(double[] gradient, double[] values) {

            final double[] result = values.clone();

            int inputOffset = 0;
            int outputOffset = 0;
            for (MultivariableTransform anArray : array) {
                int inputDimension = anArray.getInputDimension();
                int outputDimension = anArray.getOutputDimension();
                double[] tmpVal = new double[outputDimension];
                System.arraycopy(values, inputOffset, tmpVal, 0, outputDimension);
                double[] tmpGrad = new double[inputDimension];
                System.arraycopy(gradient, inputOffset, tmpGrad, 0, inputDimension);
                System.arraycopy(anArray.updateGradientLogDensity(tmpGrad, tmpVal, 0, outputDimension), 0, result, outputOffset, outputDimension);
                inputOffset += inputDimension;
                outputOffset += outputDimension;
            }
            return result;
        }

        @Override
        public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] values, int from, int to) {
            final double[] result = values.clone();

            int offset = 0;
            for (MultivariableTransform anArray : array) {
                int dim = anArray.getDimension();
                double tmpVal[] = new double[dim];
                System.arraycopy(values, offset, tmpVal, 0, dim);
                double tmpGrad[] = new double[dim];
                System.arraycopy(gradient, offset, tmpGrad, 0, dim);
                double tmpHess[] = new double[dim];
                System.arraycopy(gradient, offset, tmpHess, 0, dim);
                System.arraycopy(anArray.updateDiagonalHessianLogDensity(tmpHess, tmpGrad, tmpVal, from, to), 0, result, offset, dim);
                offset += dim;
            }
            return result;
        }

        @Override
        protected double[] updateGradientInverseUnWeightedLogDensity(double[] gradient, double[] values) {

            final double[] result = values.clone();

            int offset = 0;
            for (MultivariableTransform anArray : array) {
                int dim = anArray.getDimension();
                double tmpVal[] = new double[dim];
                System.arraycopy(values, offset, tmpVal, 0, dim);
                double tmpGrad[] = new double[dim];
                System.arraycopy(gradient, offset, tmpGrad, 0, dim);
                System.arraycopy(anArray.updateGradientInverseUnWeightedLogDensity(tmpGrad, tmpVal, 0, dim), 0, result, offset, dim);
                offset += dim;
            }
            return result;
        }

        @Override
        public double[] gradient(double[] values, int from, int to) {

            final double[] result = values.clone();

            int offset = 0;
            for (MultivariableTransform anArray : array) {
                int dim = anArray.getDimension();
                double tmp[] = new double[dim];
                System.arraycopy(values, offset, tmp, 0, dim);
                System.arraycopy(anArray.gradient(tmp, from, to), 0, result, offset, dim);
                offset += dim;
            }
            return result;
        }

        @Override
        protected double getLogJacobian(double[] values) {

            double sum = 0.0;

            int offset = 0;
            for (MultivariableTransform anArray : array) {
                int dim = anArray.getDimension();
                double tmp[] = new double[dim];
                System.arraycopy(values, offset, tmp, 0, dim);
                sum += anArray.getLogJacobian(tmp, 0, dim);
                offset += dim;
            }
            return sum;
        }

        @Override
        public double[][] updateHessianLogDensity(double[][] hessian, double[][] transformationHessian, double[] gradient, double[] value, int from, int to) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        protected double[] getGradientLogJacobianInverse(double[] values) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double[][] computeJacobianMatrixInverse(double[] values) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public String getTransformName() {
            return "MultivariateArray";
        }

        public boolean isMultivariate() {
            return true;
        }
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

        public static Transform parseTransform(XMLObject xo) {
            final Transform transform = (Transform) xo.getChild(Transform.class);
            final Transform.ParsedTransform parsedTransform
                    = (Transform.ParsedTransform) xo.getChild(Transform.ParsedTransform.class);
            if (transform == null && parsedTransform != null) return parsedTransform.transform;
            return transform;
        }

        public static MultivariableTransform parseMultivariableTransform(Object obj) {
            if (obj instanceof MultivariableTransform) {
                return (MultivariableTransform) obj;
            }
            if (obj instanceof Transform.ParsedTransform) {
                return (MultivariableTransform) ((Transform.ParsedTransform) obj).transform;
            }
            return null;
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
        INVERSE_SUM("inverseSum", new InverseSumTransform()),
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
