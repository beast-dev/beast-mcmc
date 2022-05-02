/*
 * DiscontinuousLikelihoodAtBoundary.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public interface GraphicalParameterBound {

    Parameter getParameter();

    int[] getConnectedParameterIndices(int index);

    double getFixedLowerBound(int index);

    double getFixedUpperBound(int index);

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
    }
}
