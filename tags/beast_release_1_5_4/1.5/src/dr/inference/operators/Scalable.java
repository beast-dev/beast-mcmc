/*
 * Scalable.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inference.operators;

import dr.inference.model.Bounds;
import dr.inference.model.Parameter;

/**
 * A generic interface for objects capable of scaling.
 * <p/>
 * A default impelementation for any parameter.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: UpDownOperator.java,v 1.25 2005/06/14 10:40:34 rambaut Exp $
 */

public interface Scalable {
    /**
     * @param factor scaling factor
     * @param nDims
     * @return Number of dimentions.
     * @throws OperatorFailedException
     */
    int scale(double factor, int nDims) throws OperatorFailedException;

    /**
     * @return Name for display purposes.
     */
    String getName();

    public class Default implements Scalable {
        private final Parameter parameter;

        public Default(Parameter p) {
            this.parameter = p;
        }

        public int scale(double factor, int nDims) throws OperatorFailedException {
            assert nDims <= 0;
            final int dimension = parameter.getDimension();

            for (int i = 0; i < dimension; ++i) {
                parameter.setParameterValue(i, parameter.getParameterValue(i) * factor);
            }

            final Bounds<Double> bounds = parameter.getBounds();

            for (int i = 0; i < dimension; i++) {
                final double value = parameter.getParameterValue(i);
                if (value < bounds.getLowerLimit(i) || value > bounds.getUpperLimit(i)) {
                    throw new OperatorFailedException("proposed value outside boundaries");
                }
            }
            return dimension;
        }

        public String getName() {
            return parameter.getParameterName();
        }


        /**
         *
         * This method scales values of all dimesions quitely except the last one where is also signals that
         * all values in all dimensions have been scaled. This prevents repeated firing of Events when in fact
         * only ONE parameter is changed.
         *
         */
        public int scaleAllAndNotify(double factor, int nDims) throws OperatorFailedException {

            assert nDims <= 0;
            final int dimension = parameter.getDimension();
            final int dimMinusOne = dimension-1;

            for(int i = 0; i < dimMinusOne; ++i) {
                parameter.setParameterValueQuietly(i, parameter.getParameterValue(i) * factor);
            }

            parameter.setParameterValueNotifyChangedAll(dimMinusOne, parameter.getParameterValue(dimMinusOne) * factor);

            final Bounds<Double> bounds = parameter.getBounds();

            for(int i = 0; i < dimension; i++) {
                final double value = parameter.getParameterValue(i);
                if( value < bounds.getLowerLimit(i) || value > bounds.getUpperLimit(i) ) {
                    throw new OperatorFailedException("proposed value outside boundaries");
                }
            }
            return dimension;
        }
    }
}
