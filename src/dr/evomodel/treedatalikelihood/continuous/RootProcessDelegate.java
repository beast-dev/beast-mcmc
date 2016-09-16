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

import dr.evolution.tree.Tree;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.BufferIndexHelper;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;

/**
 * @author Marc A. Suchard
 * @version $Id$
 */
public interface RootProcessDelegate {

    int getExtraPartialBufferCount();

    double calculateRootLogLikelihood(ContinuousDiffusionIntegrator cdi, int rootIndex, final double[] logLike);

//    void storeState();

//    void restoreState();

    class FullyConjugate implements RootProcessDelegate {

        private final ConjugateRootTraitPrior prior;
        private final int bufferOffset;

        public FullyConjugate(final ConjugateRootTraitPrior prior, int numTraits, int bufferOffset) {
            this.prior = prior;
            this.bufferOffset = bufferOffset;
        }

        @Override
        public int getExtraPartialBufferCount() { return 2; }

        @Override
        public double calculateRootLogLikelihood(ContinuousDiffusionIntegrator cdi, int rootIndex,
                                               final double[] logLike) {
           return 0.0;
        }
    }
}
