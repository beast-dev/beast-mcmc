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

import dr.evomodel.treedatalikelihood.BufferIndexHelper;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * @author Marc A. Suchard
 * @version $Id$
 */
public interface RootProcessDelegate extends Model {

    int getExtraPartialBufferCount();

    int getExtraMatrixBufferCount();

    void calculateRootLogLikelihood(ContinuousDiffusionIntegrator cdi, int rootIndex, int precisionIndex,
                                    final double[] logLike, boolean incrementOuterProducts, boolean isIntegratedProcess);

    double getPseudoObservations();

    int getPriorBufferIndex();

    abstract class Abstract extends AbstractModel implements RootProcessDelegate {

        protected final ConjugateRootTraitPrior prior;
        private final PrecisionType precisionType;
        private final int numTraits;

        private final int priorBufferIndexOffset;
        private final BufferIndexHelper priorBufferIndex;

        private boolean updatePrior;

        public abstract double getPseudoObservations();

        public Abstract(final ConjugateRootTraitPrior prior,
                        final PrecisionType precisionType, int numTraits,
                        int partialBufferCount, int matrixBufferCount) {

            super("RootProcessDelegate");

            this.prior = prior;
            this.precisionType = precisionType;
            this.numTraits = numTraits;

            this.priorBufferIndexOffset = partialBufferCount;
            priorBufferIndex = new BufferIndexHelper(1, 0);

            if (prior != null) {
                addModel(prior);
            }

            updatePrior = true;
        }

        @Override
        public int getExtraPartialBufferCount() { return 2; }

        @Override
        public int getExtraMatrixBufferCount() {
            return 0;
        }

        @Override
        public int getPriorBufferIndex() { return priorBufferIndexOffset + priorBufferIndex.getOffsetIndex(0); }

        @Override
        public void calculateRootLogLikelihood(ContinuousDiffusionIntegrator cdi, int rootBufferIndex, int precisionIndex,
                                               final double[] logLike, boolean incrementOuterProducts,
                                               boolean isIntegratedProcess) {

            if (updatePrior) {
                setRootPartial(cdi);
                updatePrior = false;
            }

            cdi.calculateRootLogLikelihood(rootBufferIndex, getPriorBufferIndex(), precisionIndex, logLike,
                    incrementOuterProducts, isIntegratedProcess);
        }

        private void setRootPartial(ContinuousDiffusionIntegrator cdi) {
            double[] mean = prior.getMean();
            final int dimTrait = mean.length;

            final int length = precisionType.getPartialsDimension(dimTrait);
            double[] partial = new double[length * numTraits];

            int offset = 0;
            for (int trait = 0; trait < numTraits; ++trait) {
                System.arraycopy(mean, 0, partial, offset, dimTrait);

                final double precision = getPseudoObservations();
                for (int i = 0; i < dimTrait; ++i) {
                    precisionType.fillPrecisionInPartials(partial, offset, i, precision, dimTrait);
                }

                if (precision != 0.0){
                    precisionType.fillEffDimInPartials(partial, offset, dimTrait, dimTrait);
                }

                offset += length;
            }

            priorBufferIndex.flipOffset(0);
            cdi.setPostOrderPartial(getPriorBufferIndex(), partial);
        }

        @Override
        protected void handleModelChangedEvent(Model model, Object object, int index) {
            if (model == prior) {
                updatePrior = true;
                fireModelChanged(object);
            } else {
                throw new IllegalArgumentException("Unknown submodel");
            }
        }

        @Override
        protected void storeState() {
            priorBufferIndex.storeState();
        }

        @Override
        protected void restoreState() {
            priorBufferIndex.restoreState();
        }

        @Override
        protected void acceptState() { }

        @Override
        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            throw new IllegalArgumentException("No subvariables");
        }
    }

    class Fixed extends Abstract {

        public Fixed(ConjugateRootTraitPrior prior, PrecisionType precisionType,
                     int numTraits, int partialBufferCount, int matrixBufferCount) {
            super(prior, precisionType, numTraits, partialBufferCount, matrixBufferCount);
        }

        @Override
        public double getPseudoObservations() {
            return Double.POSITIVE_INFINITY;
        }
    }

    class FullyConjugate extends Abstract {

        FullyConjugate(ConjugateRootTraitPrior prior, PrecisionType precisionType,
                              int numTraits, int partialBufferCount, int matrixBufferCount) {
            super(prior, precisionType, numTraits, partialBufferCount, matrixBufferCount);
        }

        @Override
        public double getPseudoObservations() {
            return prior.getPseudoObservations();
        }
    }
}
