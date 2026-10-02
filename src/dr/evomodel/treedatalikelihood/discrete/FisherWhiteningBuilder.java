/*
 * FisherWhiteningBuilder.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 */

package dr.evomodel.treedatalikelihood.discrete;

import org.ejml.data.Complex64F;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.EigenDecomposition;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/** Builds a rank-reduced fixed whitening from a symmetric Fisher metric. */
public final class FisherWhiteningBuilder {

    private static final double RANK_TOLERANCE = 1.0e-10;
    private static final double MINIMUM_PRECISION = 1.0e-12;

    public static final class Result {
        public final double[][] map;
        public final int rank;
        public final int nullity;
        public final double minimumPositiveEigenvalue;
        public final double maximumEigenvalue;

        Result(final double[][] map,
               final int rank,
               final int nullity,
               final double minimumPositiveEigenvalue,
               final double maximumEigenvalue) {
            this.map = map;
            this.rank = rank;
            this.nullity = nullity;
            this.minimumPositiveEigenvalue = minimumPositiveEigenvalue;
            this.maximumEigenvalue = maximumEigenvalue;
        }
    }

    private static final class EigenPair {
        final double value;
        final double[] vector;

        EigenPair(final double value, final double[] vector) {
            this.value = value;
            this.vector = vector;
        }
    }

    private FisherWhiteningBuilder() {
    }

    public static Result build(final double[][] fisher,
                               final double[][] priorMetric,
                               final double eigenFloorFraction,
                               final double nullPrecision,
                               final double trustRadius) {
        validateSquare(fisher, "Fisher metric");
        final int dimension = fisher.length;
        if (priorMetric != null) {
            validateSquare(priorMetric, "prior metric");
            if (priorMetric.length != dimension) {
                throw new IllegalArgumentException("Prior and Fisher metric dimensions differ.");
            }
        }
        if (!(eigenFloorFraction > 0.0) || !Double.isFinite(eigenFloorFraction)) {
            throw new IllegalArgumentException("Fisher eigen-floor fraction must be finite and positive.");
        }
        if (!(nullPrecision > 0.0) || !Double.isFinite(nullPrecision)) {
            throw new IllegalArgumentException("Fisher null precision must be finite and positive.");
        }
        if (!(trustRadius > 0.0) || Double.isNaN(trustRadius)) {
            throw new IllegalArgumentException("Fisher trust radius must be positive.");
        }

        final List<EigenPair> pairs = symmetricEigenpairs(fisher);
        Collections.sort(pairs, new Comparator<EigenPair>() {
            @Override
            public int compare(final EigenPair left, final EigenPair right) {
                return Double.compare(right.value, left.value);
            }
        });
        final double maximum = Math.max(0.0, pairs.get(0).value);
        final double rankThreshold = Math.max(MINIMUM_PRECISION, maximum * RANK_TOLERANCE);
        final List<EigenPair> positive = new ArrayList<EigenPair>();
        final List<EigenPair> nullPairs = new ArrayList<EigenPair>();
        for (final EigenPair pair : pairs) {
            if (pair.value > rankThreshold) {
                positive.add(pair);
            } else {
                nullPairs.add(pair);
            }
        }
        if (positive.isEmpty()) {
            throw new IllegalArgumentException("Fisher metric has no identifiable positive direction.");
        }

        final List<Double> positiveValues = new ArrayList<Double>(positive.size());
        for (final EigenPair pair : positive) {
            positiveValues.add(pair.value);
        }
        Collections.sort(positiveValues);
        final double median = positiveValues.get(positiveValues.size() / 2);
        final double floor = Math.max(MINIMUM_PRECISION, median / eigenFloorFraction);
        final double[][] map = new double[dimension][dimension];
        int column = 0;
        double minimumPositive = Double.POSITIVE_INFINITY;
        for (final EigenPair pair : positive) {
            minimumPositive = Math.min(minimumPositive, pair.value);
            putScaledColumn(map, column++, pair.vector,
                    Math.max(pair.value, floor), trustRadius);
        }

        if (!nullPairs.isEmpty()) {
            final double[][] nullBasis = columns(nullPairs, dimension);
            if (priorMetric == null) {
                for (int j = 0; j < nullPairs.size(); ++j) {
                    putScaledColumn(map, column++, nullPairs.get(j).vector,
                            nullPrecision, trustRadius);
                }
            } else {
                final double[][] projectedPrior = project(priorMetric, nullBasis);
                final List<EigenPair> priorPairs = symmetricEigenpairs(projectedPrior);
                for (final EigenPair priorPair : priorPairs) {
                    final double[] nativeVector = multiply(nullBasis, priorPair.vector);
                    putScaledColumn(map, column++, nativeVector,
                            Math.max(priorPair.value, nullPrecision), trustRadius);
                }
            }
        }
        if (column != dimension) {
            throw new IllegalStateException("Whitening basis is incomplete.");
        }
        return new Result(map, positive.size(), nullPairs.size(), minimumPositive, maximum);
    }

    private static List<EigenPair> symmetricEigenpairs(final double[][] values) {
        final int dimension = values.length;
        final DenseMatrix64F matrix = new DenseMatrix64F(dimension, dimension);
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                final double value = 0.5 * (values[i][j] + values[j][i]);
                if (!Double.isFinite(value)) {
                    throw new IllegalArgumentException("Metric contains a non-finite value.");
                }
                matrix.set(i, j, value);
            }
        }
        final EigenDecomposition<DenseMatrix64F> decomposition =
                DecompositionFactory.eig(dimension, true, true);
        if (!decomposition.decompose(matrix)) {
            throw new IllegalArgumentException("Symmetric metric eigendecomposition failed.");
        }
        final List<EigenPair> result = new ArrayList<EigenPair>(dimension);
        for (int i = 0; i < dimension; ++i) {
            final Complex64F eigenvalue = decomposition.getEigenvalue(i);
            if (Math.abs(eigenvalue.imaginary) > 1.0e-8 * (1.0 + Math.abs(eigenvalue.real))) {
                throw new IllegalArgumentException("Symmetric metric produced a complex eigenvalue.");
            }
            final DenseMatrix64F eigenvector = decomposition.getEigenVector(i);
            if (eigenvector == null) {
                throw new IllegalArgumentException("Metric eigendecomposition omitted an eigenvector.");
            }
            final double[] vector = new double[dimension];
            double normSquared = 0.0;
            for (int j = 0; j < dimension; ++j) {
                vector[j] = eigenvector.get(j, 0);
                normSquared += vector[j] * vector[j];
            }
            final double norm = Math.sqrt(normSquared);
            for (int j = 0; j < dimension; ++j) {
                vector[j] /= norm;
            }
            result.add(new EigenPair(eigenvalue.real, vector));
        }
        return result;
    }

    private static void putScaledColumn(final double[][] map,
                                        final int column,
                                        final double[] vector,
                                        final double precision,
                                        final double trustRadius) {
        double scale = 1.0 / Math.sqrt(Math.max(MINIMUM_PRECISION, precision));
        if (Double.isFinite(trustRadius)) {
            scale = Math.min(scale, trustRadius);
        }
        for (int i = 0; i < vector.length; ++i) {
            map[i][column] = scale * vector[i];
        }
    }

    private static double[][] columns(final List<EigenPair> pairs, final int rows) {
        final double[][] result = new double[rows][pairs.size()];
        for (int j = 0; j < pairs.size(); ++j) {
            for (int i = 0; i < rows; ++i) {
                result[i][j] = pairs.get(j).vector[i];
            }
        }
        return result;
    }

    private static double[][] project(final double[][] metric, final double[][] basis) {
        final int rows = metric.length;
        final int columns = basis[0].length;
        final double[][] result = new double[columns][columns];
        for (int i = 0; i < columns; ++i) {
            for (int j = 0; j < columns; ++j) {
                double value = 0.0;
                for (int r = 0; r < rows; ++r) {
                    for (int c = 0; c < rows; ++c) {
                        value += basis[r][i] * metric[r][c] * basis[c][j];
                    }
                }
                result[i][j] = value;
            }
        }
        return result;
    }

    private static double[] multiply(final double[][] matrix, final double[] vector) {
        final double[] result = new double[matrix.length];
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = 0; j < vector.length; ++j) {
                result[i] += matrix[i][j] * vector[j];
            }
        }
        return result;
    }

    private static void validateSquare(final double[][] matrix, final String name) {
        if (matrix == null || matrix.length == 0) {
            throw new IllegalArgumentException(name + " must be non-empty.");
        }
        for (final double[] row : matrix) {
            if (row == null || row.length != matrix.length) {
                throw new IllegalArgumentException(name + " must be square.");
            }
        }
    }
}
