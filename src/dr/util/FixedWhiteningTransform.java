/*
 * FixedWhiteningTransform.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 */

package dr.util;

import dr.inference.model.Parameter;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * A fixed dense affine whitening applied after an existing native transform.
 *
 * If the native transform is u = h(x), this transform defines
 *
 *     z = T^{-1} (u - uStar),
 *
 * where T and uStar are fixed.  Its inverse is
 *
 *     x = h^{-1}(uStar + T z).
 *
 * This is intended for fixed-metric local optimization: T can be constructed
 * once from a Fisher metric and then left unchanged while L-BFGS learns the
 * residual curvature.
 *
 * @author Filippo Monti
 */
public final class FixedWhiteningTransform extends Transform.MultivariableTransformWithParameter {

    private final Parameter parameter;
    private final Transform.MultivariableTransform nativeTransform;
    private final double[] reference;
    private final double[][] inverseMap;
    private final double[][] map;
    private final double logAbsDeterminantMap;

    public FixedWhiteningTransform(Parameter parameter,
                                   Transform.MultivariableTransform nativeTransform,
                                   double[] reference,
                                   double[][] inverseMap) {
        if (parameter == null || nativeTransform == null || reference == null || inverseMap == null) {
            throw new IllegalArgumentException("Fixed whitening arguments must be non-null.");
        }
        final int dimension = parameter.getDimension();
        if (nativeTransform.getInputDimension() != dimension ||
                nativeTransform.getOutputDimension() != dimension ||
                reference.length != dimension || inverseMap.length != dimension) {
            throw new IllegalArgumentException("Fixed whitening dimension mismatch.");
        }
        this.map = new double[dimension][dimension];
        for (int i = 0; i < dimension; ++i) {
            if (inverseMap[i] == null || inverseMap[i].length != dimension) {
                throw new IllegalArgumentException("Fixed whitening matrix must be square.");
            }
            for (int j = 0; j < dimension; ++j) {
                if (!Double.isFinite(inverseMap[i][j])) {
                    throw new IllegalArgumentException("Fixed whitening matrix contains a non-finite entry.");
                }
                this.map[i][j] = inverseMap[i][j];
            }
        }
        this.reference = reference.clone();
        for (double value : this.reference) {
            if (!Double.isFinite(value)) {
                throw new IllegalArgumentException("Fixed whitening reference contains a non-finite entry.");
            }
        }
        this.parameter = parameter;
        this.nativeTransform = nativeTransform;

        final DenseMatrix64F mapMatrix = new DenseMatrix64F(this.map);
        final DenseMatrix64F inverseMatrix = new DenseMatrix64F(dimension, dimension);
        if (!CommonOps.invert(mapMatrix, inverseMatrix)) {
            throw new IllegalArgumentException("Fixed whitening matrix must be nonsingular.");
        }
        this.inverseMap = new double[dimension][dimension];
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                this.inverseMap[i][j] = inverseMatrix.get(i, j);
            }
        }
        final double determinant = CommonOps.det(mapMatrix);
        if (!Double.isFinite(determinant) || determinant == 0.0) {
            throw new IllegalArgumentException("Fixed whitening matrix must be finite and nonsingular.");
        }
        this.logAbsDeterminantMap = Math.log(Math.abs(determinant));
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
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

    @Override
    public boolean isMultivariate() {
        return true;
    }

    @Override
    public String getTransformName() {
        return "fixedWhitening." + nativeTransform.getTransformName();
    }

    @Override
    public double[] transform(double[] values, int from, int to) {
        requireWholeVector(values, from, to);
        final double[] nativeValues = nativeTransform.transform(values, from, to);
        final double[] centered = new double[getDimension()];
        for (int i = 0; i < centered.length; ++i) {
            centered[i] = nativeValues[i] - reference[i];
        }
        return multiply(inverseMap, centered);
    }

    @Override
    public double[] inverse(double[] values, int from, int to) {
        requireWholeVector(values, from, to);
        final double[] nativeValues = multiply(map, values);
        for (int i = 0; i < nativeValues.length; ++i) {
            nativeValues[i] += reference[i];
        }
        return nativeTransform.inverse(nativeValues, 0, nativeValues.length);
    }

    @Override
    public double[] inverse(double[] values, int from, int to, double sum) {
        throw new UnsupportedOperationException("A fixed sum is not supported by fixed whitening.");
    }

    @Override
    public boolean isInInteriorDomain(double[] values, int from, int to) {
        requireWholeVector(values, from, to);
        if (!nativeTransform.isInInteriorDomain(values, from, to)) {
            return false;
        }
        final double[] transformed = transform(values, from, to);
        for (double value : transformed) {
            if (!Double.isFinite(value)) {
                return false;
            }
        }
        return true;
    }

    @Override
    public double[] updateGradientUnWeightedLogDensity(double[] gradient,
                                                       double[] values,
                                                       int from,
                                                       int to) {
        requireWholeVector(values, from, to);
        final double[] nativeGradient = nativeTransform.updateGradientUnWeightedLogDensity(
                gradient, values, from, to);
        return multiplyTranspose(map, nativeGradient);
    }

    @Override
    public double[] updateGradientInverseUnWeightedLogDensity(double[] gradient,
                                                              double[] values,
                                                              int from,
                                                              int to) {
        requireWholeVector(values, from, to);
        final double[] nativeValues = multiply(map, values);
        for (int i = 0; i < nativeValues.length; ++i) {
            nativeValues[i] += reference[i];
        }
        final double[] nativeGradient = nativeTransform.updateGradientInverseUnWeightedLogDensity(
                gradient, nativeValues, 0, nativeValues.length);
        return multiplyTranspose(map, nativeGradient);
    }

    @Override
    public double[] updateGradientLogDensity(double[] gradient,
                                             double[] values,
                                             int from,
                                             int to) {
        requireWholeVector(values, from, to);
        // The dense affine map has a constant Jacobian, so it contributes no
        // Jacobian-gradient term.  Only the native transform contributes.
        final double[] nativeGradient = nativeTransform.updateGradientLogDensity(
                gradient, values, from, to);
        return multiplyTranspose(map, nativeGradient);
    }

    @Override
    public double logJacobian(double[] values, int from, int to) {
        requireWholeVector(values, from, to);
        return nativeTransform.logJacobian(values, from, to) - logAbsDeterminantMap;
    }

    @Override
    public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian,
                                                    double[] gradient,
                                                    double[] value,
                                                    int from,
                                                    int to) {
        throw new UnsupportedOperationException("Dense whitening does not expose a diagonal Hessian transform.");
    }

    @Override
    public double[][] updateHessianLogDensity(double[][] hessian,
                                              double[][] transformationHessian,
                                              double[] gradient,
                                              double[] value,
                                              int from,
                                              int to) {
        throw new UnsupportedOperationException("Dense whitening Hessian transformation is not implemented.");
    }

    @Override
    public double[] gradient(double[] values, int from, int to) {
        throw new UnsupportedOperationException("Use the full Jacobian pullback for dense whitening.");
    }

    @Override
    public double[] gradientInverse(double[] values, int from, int to) {
        throw new UnsupportedOperationException("Use the full Jacobian pullback for dense whitening.");
    }

    private void requireWholeVector(double[] values, int from, int to) {
        if (values == null || from != 0 || to != values.length || values.length != getDimension()) {
            throw new IllegalArgumentException("Fixed whitening must be applied to the complete parameter vector.");
        }
    }

    private static double[] multiply(double[][] matrix, double[] vector) {
        final double[] result = new double[vector.length];
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = 0; j < vector.length; ++j) {
                result[i] += matrix[i][j] * vector[j];
            }
        }
        return result;
    }

    private static double[] multiplyTranspose(double[][] matrix, double[] vector) {
        final double[] result = new double[vector.length];
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = 0; j < vector.length; ++j) {
                result[j] += matrix[i][j] * vector[i];
            }
        }
        return result;
    }
}
