/*
 * AbstractComputedCompoundMatrix.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.inference.model;

import dr.math.matrixAlgebra.WrappedMatrix;

import java.util.List;

/**
 * Abstract base class for matrix parameters that are computed from underlying parameters
 * but cannot be directly set by users. The matrix entries are computed on-demand from
 * the underlying parameters and cached until the parameters change.
 *
 * This generalizes AbstractTransformedCompoundMatrix to support arbitrary lists of
 * underlying parameters rather than just diagonal + off-diagonal structure.
 *
 * Subclasses must implement:
 * - {@link #computeEntry(int, int)} to define how each matrix entry is computed
 * - {@link #getCachedValue(int, int)} to access cached computed values
 * - {@link #updateCache()} to recompute and cache matrix values
 *
 * The class handles:
 * - Automatic cache invalidation when underlying parameters change
 * - Proper MCMC state management (store/restore/accept)
 * - Event propagation to listeners
 *
 * @author Filippo Monti
 * @author Based on AbstractTransformedCompoundMatrix by Marc Suchard and Paul Bastide
 */
public abstract class AbstractComputedCompoundMatrix extends MatrixParameter {

    protected final int dim;
    private final CompoundParameter allParameters;

    // Cache state management
    private boolean matrixKnown = false;
    private boolean savedMatrixKnown = false;

    /**
     * Construct a computed matrix parameter from a list of underlying parameters.
     *
     * If any parameter in the list is a MatrixParameter or CompoundParameter, it will be
     * automatically flattened to its constituent leaf parameters. This ensures that changes
     * to individual elements are properly detected.
     *
     * @param name The parameter name
     * @param dimension The dimension of the square matrix
     * @param parameters The list of underlying parameters that define this matrix
     */
    protected AbstractComputedCompoundMatrix(String name,
                                             int dimension,
                                             List<Parameter> parameters) {
        super(name);

        this.dim = dimension;
        this.rowDimension = dimension;
        this.columnDimension = dimension;

        this.allParameters = new CompoundParameter(name + ".underlying");
//
//        // Flatten parameters and register them as children
//        List<Parameter> flattenedParams = flattenParameters(parameters);
//
//        for (int i = 0; i < flattenedParams.size(); i++) {
//            Parameter p = flattenedParams.get(i);
//            System.err.println("  [" + i + "] " + p.getClass().getSimpleName() + " " + p.getId());
//        }
//        System.err.println("=================================================");
//
        for (Parameter param : parameters) {
            addParameter(param);
            allParameters.addParameter(param);
        }
    }

//    private static List<Parameter> flattenParameters(List<Parameter> parameters) {
//        List<Parameter> flattened = new java.util.ArrayList<>();
//
//        for (Parameter param : parameters) {
//            if (param instanceof MatrixParameter) {
//                // MatrixParameter is a CompoundParameter of column vectors
//                // Add each column parameter individually
//                MatrixParameter matrixParam = (MatrixParameter) param;
//                for (int i = 0; i < matrixParam.getParameterCount(); i++) {
//                    flattened.add(matrixParam.getParameter(i));
//                }
//            } else if (param instanceof CompoundParameter) {
//                // CompoundParameter may contain multiple sub-parameters
//                // Add each constituent parameter individually
//                CompoundParameter compoundParam = (CompoundParameter) param;
//                for (int i = 0; i < compoundParam.getParameterCount(); i++) {
//                    flattened.add(compoundParam.getParameter(i));
//                }
//            } else {
//                // Simple parameter - add directly
//                flattened.add(param);
//            }
//        }
//
//        return flattened;
//    }

    // ================================================================
    // Abstract methods that subclasses must implement
    // ================================================================

    /**
     * Compute the (row, col) entry of the matrix from the underlying parameters.
     * This method is called when the cache needs to be updated.
     *
     * @param row The row index (0-based)
     * @param col The column index (0-based)
     * @return The computed value at position (row, col)
     */
    protected abstract double computeEntry(int row, int col);

    /**
     * Retrieve the cached value at position (row, col).
     * This method should return the previously computed value without recomputing.
     *
     * @param row The row index (0-based)
     * @param col The column index (0-based)
     * @return The cached value at position (row, col)
     */
    protected abstract double getCachedValue(int row, int col);

    /**
     * Update the cache by recomputing all matrix entries from underlying parameters.
     * This method is called when matrixKnown is false.
     *
     * Typical implementation:
     * <pre>
     * protected void updateCache() {
     *     for (int i = 0; i < dim; i++) {
     *         for (int j = 0; j < dim; j++) {
     *             cachedData[i * dim + j] = computeEntry(i, j);
     *         }
     *     }
     * }
     * </pre>
     */
    protected abstract void updateCache();

    // ================================================================
    // Public API - Matrix access
    // ================================================================

    @Override
    public final double getParameterValue(int row, int col) {
        ensureMatrixUpToDate();
        return getCachedValue(row, col);
    }

    @Override
    public final double getParameterValue(int index) {
        return getParameterValue(index / dim, index % dim);
    }

    @Override
    public final double[][] getParameterAsMatrix() {
        ensureMatrixUpToDate();
        double[][] matrix = new double[dim][dim];
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                matrix[i][j] = getCachedValue(i, j);
            }
        }
        return matrix;
    }

    @Override
    public final void setParameterValue(int index, double value) {
        throw new UnsupportedOperationException(
                "Cannot set values directly on computed matrix parameter '" + getId() + "'. " +
                        "Modify the underlying parameters instead.");
    }

    @Override
    public final void setParameterValue(int row, int column, double value) {
        throw new UnsupportedOperationException(
                "Cannot set values directly on computed matrix parameter '" + getId() + "'. " +
                        "Modify the underlying parameters instead.");
    }

    @Override
    public final void setParameterValueQuietly(int index, double value) {
        throw new UnsupportedOperationException(
                "Cannot set values directly on computed matrix parameter '" + getId() + "'. " +
                        "Modify the underlying parameters instead.");
    }

    @Override
    public final void setParameterValueQuietly(int row, int column, double value) {
        throw new UnsupportedOperationException(
                "Cannot set values directly on computed matrix parameter '" + getId() + "'. " +
                        "Modify the underlying parameters instead.");
    }

    @Override
    public final void setParameterValueNotifyChangedAll(int index, double value) {
        throw new UnsupportedOperationException(
                "Cannot set values directly on computed matrix parameter '" + getId() + "'. " +
                        "Modify the underlying parameters instead.");
    }

    @Override
    public final void setParameterValueNotifyChangedAll(int row, int column, double value) {
        throw new UnsupportedOperationException(
                "Cannot set values directly on computed matrix parameter '" + getId() + "'. " +
                        "Modify the underlying parameters instead.");
    }

    // ================================================================
    // Dimension queries
    // ================================================================

    @Override
    public final int getRowDimension() {
        return dim;
    }

    @Override
    public final int getColumnDimension() {
        return dim;
    }

    @Override
    public final int getDimension() {
        return dim * dim;
    }

    // ================================================================
    // Parameter access
    // ================================================================

    /**
     * Get the compound parameter containing all underlying parameters.
     *
     * @return The compound parameter with all underlying parameters
     */
    public final CompoundParameter getAllParameters() {
        return allParameters;
    }

    /**
     * Get the number of underlying parameters.
     *
     * @return The number of parameters that define this matrix
     */
    public final int getNumberOfUnderlyingParameters() {
        return allParameters.getParameterCount();
    }

    /**
     * Get a specific underlying parameter by index.
     *
     * @param index The parameter index
     * @return The parameter at the given index
     */
    public final Parameter getUnderlyingParameter(int index) {
        return allParameters.getParameter(index);
    }

    // ================================================================
    // Cache management
    // ================================================================

    /**
     * Ensure the cached matrix is up to date.
     * If matrixKnown is false, calls updateCache() and marks the matrix as known.
     */
    private void ensureMatrixUpToDate() {
        if (!matrixKnown) {
            updateCache();
            matrixKnown = true;
        }
    }

    /**
     * Invalidate the matrix cache.
     * This is called automatically when any underlying parameter changes.
     */
    protected final void invalidateCache() {
        matrixKnown = false;
    }

    /**
     * Check if the matrix cache is currently valid.
     *
     * @return true if the cache is up to date, false otherwise
     */
    protected final boolean isCacheValid() {
        return matrixKnown;
    }

    // ================================================================
    // Event handling
    // ================================================================

    @Override
    public void fireParameterChangedEvent() {
        // Invalidate cache when this parameter fires an event
        // (which happens when any child parameter changes via CompoundParameter.variableChangedEvent)
        matrixKnown = false;

        // Propagate event upward to any listeners (e.g., MultivariateElasticModel)
        super.fireParameterChangedEvent();
    }

    // ================================================================
    // MCMC state management
    // ================================================================

    @Override
    protected void storeValues() {
        super.storeValues();
        savedMatrixKnown = matrixKnown;
    }

    @Override
    protected void restoreValues() {
        super.restoreValues();
        matrixKnown = savedMatrixKnown;
    }

    @Override
    protected void acceptValues() {
        super.acceptValues();
        // matrixKnown doesn't need special handling on accept
    }

    // ================================================================
    // Reporting
    // ================================================================

    @Override
    public String getReport() {
        return new WrappedMatrix.ArrayOfArray(getParameterAsMatrix()).toString();
    }

    @Override
    public double[] getAttributeValue() {
        double[] values = new double[dim * dim];
        int index = 0;
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                values[index++] = getParameterValue(i, j);
            }
        }
        return values;
    }

    @Override
    public String getDimensionName(int index) {
        int row = index / dim;
        int col = index % dim;
        return getId() + "[" + row + "," + col + "]";
    }
}