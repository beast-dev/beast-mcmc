/*
 * GraphicalParameterBound.java
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

/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public interface GraphicalParameterBound extends GeneralBoundsProvider {

    Parameter getParameter();

    int[] getConnectedParameterIndices(int index);

    double getFixedLowerBound(int index);

    double getFixedUpperBound(int index);

    /**
     * The orientation of the constraint between two connected coordinates. A tie (x_i == x_j) carries no order, so
     * the reflective engine asks the provider which side of the wall x_i == x_j is the constrained space.
     *
     * @param i a coordinate index
     * @param j a coordinate connected to i (an entry of getConnectedParameterIndices(i))
     * @return true if the constraint is x_i &lt;= x_j, false if it is x_i &gt;= x_j
     * @throws UnsupportedOperationException by default; a provider without connected coordinates never gets the call
     */
    default boolean isBelow(int i, int j) {
        throw new UnsupportedOperationException(getClass().getName() + " does not give the orientation of "
                + "connected coordinates");
    }

    /**
     * Checks that the provider's parameter maps to coordinates the way the provider assumes. The default accepts any
     * parameter; a provider whose constraints depend on that mapping (for example, coordinate i being the height of
     * tree node i + #tips) overrides it.
     *
     * @throws IllegalArgumentException if the provider's parameter does not map to coordinates the way the provider
     *                                  assumes
     */
    default void checkParameterMapping() {
    }

    /**
     * The same check for the parameter that a sampler moves under these bounds, which need not be getParameter()
     * (the reflective HMC operator takes its parameter from its own XML element).
     *
     * @param parameter the parameter whose coordinates are moved under these bounds
     * @throws IllegalArgumentException if parameter does not map to coordinates the way the provider assumes
     */
    default void checkParameterMapping(Parameter parameter) {
    }

    @Override
    default int countTiedConstraints(double[] values) {
        int count = 0;
        for (int i = 0; i < values.length; i++) {
            if (values[i] == getFixedLowerBound(i) || values[i] == getFixedUpperBound(i)) {
                ++count;
            }
            final int[] connectedIndices = getConnectedParameterIndices(i);
            if (connectedIndices != null) {
                for (int j : connectedIndices) {
                    if (j > i && values[i] == values[j]) {
                        ++count;
                    }
                }
            }
        }
        return count;
    }

    class FixedBound implements GraphicalParameterBound {

        private final Parameter parameter;
        private final Bounds<Double> bounds;

        public FixedBound(Parameter parameter) {
            this.parameter = parameter;
            this.bounds = parameter.getBounds();
        }

        @Override
        public Parameter getParameter() {
            return parameter;
        }

        @Override
        public int[] getConnectedParameterIndices(int index) {
            return null;
        }

        @Override
        public double getFixedLowerBound(int index) {
            return bounds.getLowerLimit(index);
        }

        @Override
        public double getFixedUpperBound(int index) {
            return bounds.getUpperLimit(index);
        }

        @Override
        public boolean isWithinBounds(double[] values) {
            for (int i = 0; i < values.length; i++) {
                final double value = values[i];
                if (!(value >= getFixedLowerBound(i) && value <= getFixedUpperBound(i))) {
                    return false;
                }
            }
            return true;
        }
    }
}
