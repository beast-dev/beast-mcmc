/*
 * ContinuousDiffusionIntegrator.java
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

package dr.evomodel.treedatalikelihood.continuous.cdi;

/**
 * @author Marc A. Suchard
 */
public interface ContinuousDiffusionIntegrator {
    int OPERATION_TUPLE_SIZE = 7;
    int NONE = -1;

    void finalize() throws Throwable;

    void setPartialMean(int bufferIndex, final double[] mean);

    void getPartialMean(int bufferIndex, final double[] mean);

    void setPartialPrecision(int bufferIndex, final double[] precision);

    void getPartialPrecision(int bufferIndex, final double[] precision);

    void setDiffusionPrecision(int diffusionIndex, final double[] matrix);

    void updatePartials(final int[] operations, int operationCount);

    InstanceDetails getDetails();


    class Basic implements ContinuousDiffusionIntegrator {

        private int instance = -1;
        private InstanceDetails details = new InstanceDetails();

        private final int numTraits;
        private final int dimMean;
        private final int dimPrecision;
        private final int bufferCount;
        private final int diffusionCount;

        public Basic(
                final int numTraits,
                final int dimMean,
                final int dimPrecision,
                final int bufferCount,
                final int diffusionCount
        ) {

            this.numTraits = numTraits;
            this.dimMean = dimMean;
            this.dimPrecision = dimPrecision;
            this.bufferCount = bufferCount;
            this.diffusionCount = diffusionCount;

            assert(dimPrecision == 1
                    || dimPrecision == dimMean
                    || dimPrecision == dimMean * dimMean);
            assert(bufferCount > 0);
            assert(diffusionCount > 0);

        }

        public void finalize() throws Throwable {
            super.finalize();
        }

        public void setPartialMean(int bufferIndex, final double[] mean) {

        }

        public void getPartialMean(int bufferIndex, final double[] mean) {

        }

        public void setPartialPrecision(int bufferIndex, final double[] precision) {

        }

        public void getPartialPrecision(int bufferIndex, final double[] precision) {

        }

        public void setDiffusionPrecision(int precisionIndex, final double[] matrix) {

        }

        public void updatePartials(final int[] operations, int operationCount) {

        }

        public InstanceDetails getDetails() {
            return details;
        }
    }
}
