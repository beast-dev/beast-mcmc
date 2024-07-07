/*
 * AdaptableVector.java
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

package dr.math;

import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;

/**
 * @author Marc A. Suchard
 * @author Guy Baele
 * @author Xiang Ji
 */
public interface AdaptableVector {

    public ReadableVector getMean();

    public int getUpdateCount();

    public void update(ReadableVector x);

    class Default implements AdaptableVector {
        final protected int dim;
        final private double[] oldMeans;
        final private double[] newMeans;

        protected int updates;

        public Default(int dim) {
            this.dim = dim;
            this.oldMeans = new double[dim];
            this.newMeans = new double[dim];

            updates = 0;
        }

        @Override
        public ReadableVector getMean() {
            return new WrappedVector.Raw(newMeans);
        }

        @Override
        public int getUpdateCount() {
            return updates;
        }

        @Override
        public void update(ReadableVector x) {
            updates++;
            for (int i = 0; i < dim; i++) { //TODO: swap the pointers
                oldMeans[i] = newMeans[i];
                newMeans[i] = ((oldMeans[i] * (updates - 1)) + x.get(i)) / updates;
            }
        }

        public double getOldMeans(int index) {
            return oldMeans[index];
        }

        public double getNewMeans(int index) {
            return newMeans[index];
        }
    }

    class LimitedMemory implements AdaptableVector {

        final private int dim;
        final private double[][] meanQueue;
        private double[] mean;
        final private int queueSize;
        private int updateIndex;
        private int updates;

        public LimitedMemory(int dim, int queueSize) {
            this.dim = dim;
            this.queueSize = queueSize;
            this.meanQueue = new double[queueSize][dim];
            this.mean = new double[dim];
            this.updateIndex = 0;
            this.updates = 0;
        }

        @Override
        public ReadableVector getMean() {
            return new WrappedVector.Raw(mean);
        }

        @Override
        public int getUpdateCount() {
            return updates;
        }

        @Override
        public void update(ReadableVector x) {
            for (int i = 0; i < dim; i++) {
                double increment = (x.get(i) - meanQueue[updateIndex][i]) / ((double) queueSize);
                mean[i] += increment;
                meanQueue[updateIndex][i] = x.get(i);
            }
            updates++;
            updateIndex = (updateIndex + 1) % queueSize;
        }
    }

    class AdaptableVariance extends Default {

        final private double[] meanSquaredValues;
        final private double[] variance;

        public AdaptableVariance(int dim) {
            super(dim);
            this.meanSquaredValues = new double[dim];
            this.variance = new double[dim];
        }

        @Override
        public void update(ReadableVector x) {
            super.update(x);
            for (int i = 0; i < dim; i++) {
                meanSquaredValues[i] = ((updates - 1.0) * meanSquaredValues[i] + x.get(i) * x.get(i)) / updates;
                variance[i] = Math.abs(meanSquaredValues[i] - getNewMeans(i) * getNewMeans(i));
            }
        }

        public double[] getVariance() {
            return variance.clone();
        }
    }

}
