/*
 * RootProcessDelegate.java
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

package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;

/**
 * @author Marc A. Suchard
 * @version $Id$
 */
public interface RootProcessDelegate {

    int getExtraPartialBufferCount();

    int getExtraMatrixBufferCount();

    void calculateRootLogLikelihood(ContinuousDiffusionIntegrator cdi, int rootIndex, final double[] logLike,
                                    boolean incrementOuterProducts);

    void setRootPartial(ContinuousDiffusionIntegrator cdi);

    int getDegreesOfFreedom();

    abstract class Abstract implements RootProcessDelegate {

        protected final ConjugateRootTraitPrior prior;
        private final int priorBufferIndex;
        private final int numTraits;

        protected abstract double getPseudoObservations();

        public Abstract(final ConjugateRootTraitPrior prior, int numTraits,
                        int partialBufferCount, int matrixBufferCount) {
            this.prior = prior;
            this.numTraits = numTraits;

            this.priorBufferIndex = partialBufferCount;
        }

        @Override
        public int getExtraPartialBufferCount() {
            return 2; // TODO Why does 1 not work?
        }

        @Override
        public int getExtraMatrixBufferCount() {
            return 0;
        }

        @Override
        public void calculateRootLogLikelihood(ContinuousDiffusionIntegrator cdi, int rootBufferIndex,
                                               final double[] logLike, boolean incrementOuterProducts) {
            cdi.calculateRootLogLikelihood(rootBufferIndex, priorBufferIndex, logLike, incrementOuterProducts);
        }

        @Override
        public void setRootPartial(ContinuousDiffusionIntegrator cdi) {
            double[] mean = prior.getMean();
            final int dimTrait = mean.length;

            double[] partial = new double[(dimTrait + 1) * numTraits];

            int offset = 0;
            for (int trait = 0; trait < numTraits; ++trait) {
                System.arraycopy(mean, 0, partial, offset, dimTrait);
                offset += dimTrait;
                partial[offset] = getPseudoObservations();
                offset += 1;
            }

            cdi.setPartial(priorBufferIndex, partial);
        }
    }

    class Fixed extends Abstract {

        public Fixed(ConjugateRootTraitPrior prior, int numTraits, int partialBufferCount, int matrixBufferCount) {
            super(prior, numTraits, partialBufferCount, matrixBufferCount);
        }

        @Override
        protected double getPseudoObservations() {
            return Double.POSITIVE_INFINITY;
        }

        @Override
        public int getDegreesOfFreedom() { return 0; }
    }

    class FullyConjugate extends Abstract {

        public FullyConjugate(ConjugateRootTraitPrior prior, int numTraits, int partialBufferCount, int matrixBufferCount) {
            super(prior, numTraits, partialBufferCount, matrixBufferCount);
        }

        @Override
        protected double getPseudoObservations() {
            return prior.getPseudoObservations();
        }

        @Override
        public int getDegreesOfFreedom() { return 1; }
    }
}
